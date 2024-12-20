use crate::quadrature::IntegralEstimate;

mod adaptive;
mod basic;
mod rule;

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub enum Kind {
    WrongDimensionality,
    MinimumPointsLargerThanMaximum(usize, usize),
    MaximumPointsInsufficient(usize, usize),
}

#[derive(Debug)]
pub struct Error {
    kind: Kind,
    integral: IntegralEstimate,
}

impl Error {
    pub(crate) const fn new(kind: Kind, integral: IntegralEstimate) -> Self {
        Self { kind, integral }
    }

    pub(crate) const fn unevaluated(kind: Kind) -> Self {
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
        match self.kind() {
            Kind::WrongDimensionality => {
                write!(
                    f,
                    "Wrong dimensionality. The multidimensional quadrature routine is optimised to handle dimensionality N in the range 2 <= N <= 15. For 1 dimensional integrals, consider using the specialised 1 dimensional integration routines [`Adaptive`] and/or [`AdaptiveSingularity`]."
                )
            }

            Kind::MinimumPointsLargerThanMaximum(min, max) => {
                write!(
                    f,
                    "Minimum points larger than maximum. The minimum number of points provided, {min}, is larger than the maximum number of points provided, {max}. The maximum number of points must exceed the minimum."
                )
            }

            Kind::MaximumPointsInsufficient(max, min_max) => {
                write!(
                    f,
                    "Maximum number of points insufficient. The minimum value for the maximum number of points is set by the dimensionality of the function, N, to be 2^N + 2 * N * (N + 1) + 1), which for the provided function is {min_max}. The provided maximum number of points, {max}, is less than this minimum value."
                )
            }
        }
    }
}

impl std::error::Error for Error {}
