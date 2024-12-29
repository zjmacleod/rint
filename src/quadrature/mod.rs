mod adaptive;
mod basic;
mod estimate;
mod integrator;
mod region;
mod rule;
mod singularity;

#[cfg(test)]
mod tests;

pub(crate) use integrator::Integrator;
pub(crate) use region::Region;

pub use adaptive::Adaptive;
pub use basic::Basic;
pub use singularity::AdaptiveSingularity;

pub use estimate::IntegralEstimate;
pub use rule::Rule;

use crate::Limits;
use crate::ScalarF64;

/// User selected tolerance type.
///
/// The adaptive routines will return the first approximation, `result`, to the integral which has an
/// absolute `error` smaller than the tolerance set by the choice of [`Tolerance`], where
/// * [`Tolerance::Absolute(abserr)`] specifies an absolute error and returns final [`IntegralEstimate`] when `error <= abserr`,
/// * [`Tolerance::Relative(relerr)`] specifies a relative error and returns final [`IntegralEstimate`] when `error <= relerr * abs(result)`,  
/// * [`Tolerance::Either{ abserr, relerr }`] to return a result as soon as _either_ the relative or absolute error bound has been satisfied.
///
/// [`Tolerance::Absolute(abserr)`]: crate::quadrature::Tolerance#variant.Absolute
/// [`Tolerance::Relative(relerr)`]: crate::quadrature::Tolerance#variant.Relative
/// [`Tolerance::Either{ abserr, relerr }`]: crate::quadrature::Tolerance#variant.Either
pub enum Tolerance {
    /// Specify an absolute error and returns final [`IntegralEstimate`] when `error <= abserr`.
    Absolute(f64),
    /// Specify a relative error and return final [`IntegralEstimate`] when `error <= relerr * abs(result)`.
    Relative(f64),
    /// Specify _both_ an absolute and relative error and return final [`IntegralEstimate`] as soon as _either_ the relative or absolute error bound has been satisfied.
    Either { absolute: f64, relative: f64 },
}

impl Tolerance {
    #[must_use]
    pub(crate) fn tolerance<T: ScalarF64>(&self, integral_value: T) -> f64 {
        match *self {
            Tolerance::Absolute(v) => v,
            Tolerance::Relative(v) => v * integral_value.abs(),
            Tolerance::Either { absolute, relative } => {
                f64::max(absolute, relative * integral_value.abs())
            }
        }
    }
}

#[derive(Debug)]
pub struct Error<T: ScalarF64> {
    kind: Kind,
    integral: IntegralEstimate<T>,
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub enum Kind {
    MaximumIterationsReached,
    RoundoffErrorDetected,
    BadIntegrandBehaviour { limits: Limits },
    DoesNotConverge,
    DivergentOrSlowlyConverging,
    UninitialisedWorkspace,
    AbsoluteBoundNegativeOrZero,
    RelativeBoundTooSmall,
    InvalidTolerance,
}

impl<T: ScalarF64> Error<T> {
    pub(crate) fn new(kind: Kind, integral: IntegralEstimate<T>) -> Self {
        Self { kind, integral }
    }

    pub(crate) fn unevaluated(kind: Kind) -> Self {
        let output = IntegralEstimate::new();
        Error::new(kind, output)
    }

    /// Return the error [`Kind`] which was encountered.
    #[must_use]
    pub fn kind(&self) -> Kind {
        self.kind
    }

    /// Return a reference to the best [`IntegralEstimate`] which was calculated before an error
    /// occurred.
    #[must_use]
    pub fn estimate(&self) -> &IntegralEstimate<T> {
        &self.integral
    }

    /// Return the best estimate of the integral value which was calculated before an error occurred.
    #[must_use]
    pub fn result(&self) -> T {
        self.integral.result()
    }

    /// Return the best estimate of the integral error which was calculated before an error occurred.
    #[must_use]
    pub fn error(&self) -> f64 {
        self.integral.error()
    }

    /// Return the number of iterations used by the integrator before an error occurred.
    #[must_use]
    pub fn iterations(&self) -> usize {
        self.integral.iterations()
    }

    /// Return the number of function evaluations used by the integrator before an error occurred.
    #[must_use]
    pub fn function_evaluations(&self) -> usize {
        self.integral.function_evaluations()
    }
}

impl<T: ScalarF64> std::fmt::Display for Error<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let result = self.result();
        let error = self.error();
        let iterations = self.iterations();
        match self.kind() {
            Kind::MaximumIterationsReached => {
                write!(
                    f,
                    "Maximum number of iterations/subdivisions reached. Try increasing max_iterations. If this yields no improvement it is advised to analyse the integrand to determine integration difficulties. If the position of a local difficulty can be determined, one may gain from splitting the total integration interval and calling the integrator on each sub-interval.\nresult:\t{result:?}\nerror\t{error:.10e}\niterations:\t{iterations}."
                )
            }

            Kind::RoundoffErrorDetected => {
                write!(
                    f,
                    "Roundoff error detected. This prevents the requested tolerance from being achieved and the returned error may be under-estimated.\nresult:\t{result:?}\nerror\t{error:.10e}\niterations:\t{iterations}."
                )
            }

            Kind::BadIntegrandBehaviour { limits } => {
                let lower = limits.lower();
                let upper = limits.upper();
                write!(
                    f,
                    "Extremely bad integrand behaviour. Possible non-integrable singularity, divergence, or discontinuity detected between ({lower},{upper}).\nresult:\t{result:?}\nerror\t{error:.10e}\niterations:\t{iterations}."
                )
            }

            Kind::DoesNotConverge => {
                write!(
                    f,
                    "The algorithm does not converge. Roundoff error is detected in the extrapolation table. It is assumed that the requested tolerance cannot be achieved and the returned result is the best which can be obtained.\nresult:\t{result:?}\nerror\t{error:.10e}\niterations:\t{iterations}."
                )
            }

            Kind::DivergentOrSlowlyConverging => {
                write!(
                    f,
                    "The integral is probably divergent or slowly convergent. NOTE: divergence can also occur with any other error kind.\nresult:\t{result:?}\nerror\t{error:.10e}\niterations:\t{iterations}."
                )
            }

            Kind::UninitialisedWorkspace => {
                write!(f, "The integration Workspace was not properly initialised. This error should not be possible. If this error is returned, contact the crate maintainers.")
            }

            Kind::AbsoluteBoundNegativeOrZero => {
                write!(
                    f,
                    "Invalid error bound: relative bound must be non-zero and positive."
                )
            }

            Kind::RelativeBoundTooSmall => {
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

impl<T: ScalarF64> std::error::Error for Error<T> {}

pub(crate) fn rescale_error(error: f64, result_abs: f64, result_asc: f64) -> f64 {
    let mut error = error;

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
