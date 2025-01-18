#![allow(unused)]
use crate::IntegralEstimate;
use crate::ScalarF64;

mod adaptive;
mod basic;
mod generator;
mod geometry;
mod integrator;
mod region;
mod rule;

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub enum Kind {
    WrongDimensionality,
    MinimumPointsLargerThanMaximum(usize, usize),
    MaximumPointsInsufficient(usize, usize),
}

#[derive(Debug)]
pub struct MultiDimensionalError<T: ScalarF64> {
    kind: Kind,
    integral: IntegralEstimate<T>,
}

impl<T: ScalarF64> MultiDimensionalError<T> {
    pub(crate) const fn new(kind: Kind, integral: IntegralEstimate<T>) -> Self {
        Self { kind, integral }
    }

    pub(crate) fn unevaluated(kind: Kind) -> Self {
        let output = IntegralEstimate::new();
        MultiDimensionalError::new(kind, output)
    }

    #[must_use]
    pub const fn kind(&self) -> Kind {
        self.kind
    }

    #[must_use]
    pub const fn estimate(&self) -> &IntegralEstimate<T> {
        &self.integral
    }

    #[must_use]
    pub const fn result(&self) -> T {
        self.integral.result()
    }

    #[must_use]
    pub const fn error(&self) -> f64 {
        self.integral.error()
    }

    #[must_use]
    pub const fn iterations(&self) -> usize {
        self.integral.iterations()
    }

    #[must_use]
    pub const fn function_evaluations(&self) -> usize {
        self.integral.function_evaluations()
    }
}

impl<T: ScalarF64> std::fmt::Display for MultiDimensionalError<T> {
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

impl<T: ScalarF64> std::error::Error for MultiDimensionalError<T> {}

#[inline]
pub(crate) const fn two_pow_n(n: usize) -> usize {
    let mut exp = n;

    // Never need to check this since NDIM > 2
    //if exp == 0 {
    //    return 1;
    //}
    let mut base = 2;
    let mut acc = 1;

    while exp > 1 {
        if (exp & 1) == 1 {
            acc *= base;
        }
        exp /= 2;
        base *= base;
    }

    acc * base
}

#[inline]
pub(crate) const fn two_pow_n_f64(n: usize) -> f64 {
    two_pow_n(n) as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_two_pow_n() {
        assert_eq!(2, two_pow_n(1));
        assert_eq!(4, two_pow_n(2));
        assert_eq!(8, two_pow_n(3));
        assert_eq!(16, two_pow_n(4));
        assert_eq!(32, two_pow_n(5));
        assert_eq!(64, two_pow_n(6));
        assert_eq!(128, two_pow_n(7));
        assert_eq!(256, two_pow_n(8));
        assert_eq!(512, two_pow_n(9));
        assert_eq!(1024, two_pow_n(10));
        assert_eq!(2048, two_pow_n(11));
        assert_eq!(4096, two_pow_n(12));
        assert_eq!(8192, two_pow_n(13));
        assert_eq!(16384, two_pow_n(14));
        assert_eq!(32768, two_pow_n(15));
    }

    #[test]
    fn test_two_pow_n_f64() {
        assert_eq!(2.0, two_pow_n_f64(1));
        assert_eq!(4.0, two_pow_n_f64(2));
        assert_eq!(8.0, two_pow_n_f64(3));
        assert_eq!(16.0, two_pow_n_f64(4));
        assert_eq!(32.0, two_pow_n_f64(5));
        assert_eq!(64.0, two_pow_n_f64(6));
        assert_eq!(128.0, two_pow_n_f64(7));
        assert_eq!(256.0, two_pow_n_f64(8));
        assert_eq!(512.0, two_pow_n_f64(9));
        assert_eq!(1024.0, two_pow_n_f64(10));
        assert_eq!(2048.0, two_pow_n_f64(11));
        assert_eq!(4096.0, two_pow_n_f64(12));
        assert_eq!(8192.0, two_pow_n_f64(13));
        assert_eq!(16384.0, two_pow_n_f64(14));
        assert_eq!(32768.0, two_pow_n_f64(15));
    }
}
