mod adaptive;
pub(crate) mod basic;
mod estimate;
mod rule;
mod singularity;

#[cfg(test)]
mod tests;

pub use adaptive::Adaptive;
pub use basic::Basic;
pub use estimate::IntegralEstimate;
pub use rule::Rule;
pub use singularity::AdaptiveSingularity;

use crate::Limits;

/// User selected error bound type.
pub enum ErrorBound {
    Absolute(f64),
    Relative(f64),
    Either { absolute: f64, relative: f64 },
}

impl ErrorBound {
    #[must_use]
    pub fn tolerance(&self, integral_value: f64) -> f64 {
        match *self {
            ErrorBound::Absolute(v) => v,
            ErrorBound::Relative(v) => v * integral_value.abs(),
            ErrorBound::Either { absolute, relative } => {
                f64::max(absolute, relative * integral_value.abs())
            }
        }
    }
}

#[derive(Debug)]
pub struct Error {
    kind: Kind,
    integral: IntegralEstimate,
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub enum Kind {
    MaximumIterationsReached,
    RoundoffErrorDetected,
    BadIntegrandBehaviour { limits: Limits },
    DoesNotConverge,
    DivergentOrSlowlyConverging,
    UninitialisedWorkspace,
    RelativeBoundNegativeOrZero,
    AbsoluteBoundTooSmall,
    InvalidTolerance,
}

impl Error {
    pub(crate) fn new(kind: Kind, integral: IntegralEstimate) -> Self {
        Self { kind, integral }
    }

    pub(crate) fn unevaluated(kind: Kind) -> Self {
        let output = IntegralEstimate::new();
        Error::new(kind, output)
    }

    #[must_use]
    pub fn kind(&self) -> Kind {
        self.kind
    }

    #[must_use]
    pub fn estimate(&self) -> &IntegralEstimate {
        &self.integral
    }

    #[must_use]
    pub fn result(&self) -> f64 {
        self.integral.result()
    }

    #[must_use]
    pub fn error(&self) -> f64 {
        self.integral.error()
    }

    #[must_use]
    pub fn iterations(&self) -> usize {
        self.integral.iterations()
    }

    #[must_use]
    pub fn function_evaluations(&self) -> usize {
        self.integral.function_evaluations()
    }
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let result = self.result();
        let error = self.error();
        let iterations = self.iterations();
        match self.kind() {
            Kind::MaximumIterationsReached => {
                write!(
                    f,
                    "Maximum number of iterations/subdivisions reached. Try increasing max_iterations. If this yields no improvement it is advised to analyse the integrand to determine integration difficulties. If the position of a local difficulty can be determined, one may gain from splitting the total integration interval and calling the integrator on each sub-interval.\nresult:\t{result:.10e}\nerror\t{error:.10e}\niterations:\t{iterations}."
                )
            }

            Kind::RoundoffErrorDetected => {
                write!(
                    f,
                    "Roundoff error detected. This prevents the requested tolerance from being achieved and the returned error may be under-estimated.\nresult:\t{result:.10e}\nerror\t{error:.10e}\niterations:\t{iterations}."
                )
            }

            Kind::BadIntegrandBehaviour { limits } => {
                let lower = limits.lower();
                let upper = limits.upper();
                write!(
                    f,
                    "Extremely bad integrand behaviour. Possible non-integrable singularity, divergence, or discontinuity detected between ({lower},{upper}).\nresult:\t{result:.10e}\nerror\t{error:.10e}\niterations:\t{iterations}."
                )
            }

            Kind::DoesNotConverge => {
                write!(
                    f,
                    "The algorithm does not converge. Roundoff error is detected in the extrapolation table. It is assumed that the requested tolerance cannot be achieved and the returned result is the best which can be obtained.\nresult:\t{result:.10e}\nerror\t{error:.10e}\niterations:\t{iterations}."
                )
            }

            Kind::DivergentOrSlowlyConverging => {
                write!(
                    f,
                    "The integral is probably divergent or slowly convergent. NOTE: divergence can also occur with any other error kind.\nresult:\t{result:.10e}\nerror\t{error:.10e}\niterations:\t{iterations}."
                )
            }

            Kind::UninitialisedWorkspace => {
                write!(f, "The integration Workspace was not properly initialised. This error should not be possible. If this error is returned, contact the crate maintainers.")
            }

            Kind::RelativeBoundNegativeOrZero => {
                write!(
                    f,
                    "Invalid error bound: relative bound must be non-zero and positive."
                )
            }

            Kind::AbsoluteBoundTooSmall => {
                write!(
                    f,
                    "Invalid error bound: absolute bound must be larger than 50.0 * f64::EPSILON."
                )
            }

            Kind::InvalidTolerance => {
                write!(f, "Invalid tolerance: relative bound must be non-zero and positive and absolute bound must be larger than 50.0 * f64::EPSILON.")
            }
        }
    }
}

impl std::error::Error for Error {}

pub(crate) fn rescale_error(error: f64, result_abs: f64, result_asc: f64) -> f64 {
    let mut error = error.abs();

    if result_asc != 0.0 && error != 0.0 {
        let scale = (200.0 * error / result_asc).powf(1.5);

        if scale < 1.0 {
            error = result_asc * scale;
        } else {
            error = result_asc;
        }
    }

    if result_abs > f64::MIN_POSITIVE / (50.0 * f64::EPSILON) {
        let min_error = 50.0 * f64::EPSILON * result_abs;

        if min_error > error {
            error = min_error;
        }
    }

    error
}

#[inline]
pub(crate) fn subinterval_too_small(limits: Limits) -> bool {
    let lower = limits.lower();
    let upper = limits.upper();
    let midpoint = limits.centre();

    let eps = f64::EPSILON;
    let min = f64::MIN_POSITIVE;

    let tmp = (1.0 + 100.0 * eps) * (midpoint.abs() + 1000.0 * min);

    lower.abs() <= tmp && upper.abs() <= tmp
}
