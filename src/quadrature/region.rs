use num_traits::Zero;
use std::cmp::Ordering;

use crate::quadrature::integrator::Integrator;
use crate::quadrature::rule::Rule;
use crate::IntegralEstimate;
use crate::Integrand;
use crate::Limits;
use crate::ScalarF64;

#[derive(Debug)]
pub(crate) struct Region<T> {
    pub(crate) result: T,
    pub(crate) error: f64,
    pub(crate) result_abs: f64,
    pub(crate) result_asc: f64,
    pub(crate) limits: Limits,
}

impl<T: ScalarF64> PartialEq for Region<T> {
    fn eq(&self, other: &Self) -> bool {
        (self.result == other.result)
            && (self.error == other.error)
            && (self.result_abs == other.result_abs)
            && (self.result_asc == other.result_asc)
            && (self.limits == other.limits)
    }
}

impl<T: ScalarF64> Eq for Region<T> {}

impl<T: ScalarF64> PartialOrd for Region<T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<T: ScalarF64> Ord for Region<T> {
    fn cmp(&self, other: &Self) -> Ordering {
        let mut ordering = self.total_cmp_error(other);
        if let Ordering::Equal = ordering {
            ordering = self.total_cmp_interval_length(other);
        }
        ordering
    }
}

impl<T: ScalarF64> Region<T> {
    pub(crate) fn unevaluated() -> Self {
        let result = T::zero();

        Self {
            error: 0.0,
            result,
            result_abs: 0.0,
            result_asc: 0.0,
            limits: Limits::new(0.0, 0.0),
        }
    }

    pub(crate) const fn with_result(mut self, result: T) -> Self {
        self.result = result;
        self
    }

    pub(crate) const fn with_result_asc(mut self, result_asc: f64) -> Self {
        self.result_asc = result_asc;
        self
    }

    pub(crate) const fn with_result_abs(mut self, result_abs: f64) -> Self {
        self.result_abs = result_abs;
        self
    }

    pub(crate) const fn with_error(mut self, error: f64) -> Self {
        self.error = error;
        self
    }

    pub(crate) const fn with_limits(mut self, limits: Limits) -> Self {
        self.limits = limits;
        self
    }

    #[must_use]
    #[inline]
    pub(crate) const fn result(&self) -> T {
        self.result
    }

    #[must_use]
    #[inline]
    pub(crate) const fn result_asc(&self) -> f64 {
        self.result_asc
    }

    #[must_use]
    #[inline]
    pub(crate) const fn result_abs(&self) -> f64 {
        self.result_abs
    }

    #[must_use]
    #[inline]
    pub(crate) const fn error(&self) -> f64 {
        self.error
    }

    #[must_use]
    pub(crate) const fn limits(&self) -> Limits {
        self.limits
    }

    pub(crate) fn total_cmp_error(&self, other: &Self) -> Ordering {
        self.error.total_cmp(&other.error)
    }

    pub(crate) fn total_cmp_interval_length(&self, other: &Self) -> Ordering {
        let inverse_length = 1.0 / (self.limits.width()).abs();
        let other_inverse_length = 1.0 / (other.limits.width()).abs();
        inverse_length.total_cmp(&other_inverse_length)
    }

    #[allow(clippy::needless_borrow)]
    pub(crate) fn bisect<I: Integrand>(&self, function: &I, rule: &Rule) -> [Region<I::Scalar>; 2] {
        let [lower, upper] = self.limits.bisect();
        let lower_integral = Integrator::new(&function, &rule, lower).integrate();

        let upper_integral = Integrator::new(&function, &rule, upper).integrate();

        [lower_integral, upper_integral]
    }

    pub(crate) fn positivity(&self) -> bool {
        self.result.abs() >= (1.0 - 50.0 * f64::EPSILON) * self.result_abs
    }

    pub(crate) fn abs_interval_length(&self) -> f64 {
        self.limits.width().abs()
    }

    pub(crate) fn estimate(&self, iterations: usize, evaluations: usize) -> IntegralEstimate<T> {
        let result = self.result();
        let error = self.error();
        IntegralEstimate::new()
            .with_result(result)
            .with_error(error)
            .with_iterations(iterations)
            .with_evaluations(evaluations)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ordering_of_basic_internal() {
        use std::collections::BinaryHeap;

        let a: Region<f64> = Region {
            error: 2.0,
            result: 1.0,
            result_abs: 1.0,
            result_asc: 1.0,
            limits: Limits::new(0.0, 1.0),
        };
        let b: Region<f64> = Region {
            error: 1.533,
            result: 1.0,
            result_abs: 1.0,
            result_asc: 1.0,
            limits: Limits::new(0.0, 1.0),
        };
        let c: Region<f64> = Region {
            error: 1.533,
            result: 1.0,
            result_abs: 1.0,
            result_asc: 1.0,
            limits: Limits::new(0.5, 1.0),
        };
        let d: Region<f64> = Region {
            error: 1.60,
            result: 1.0,
            result_abs: 1.0,
            result_asc: 1.0,
            limits: Limits::new(0.0, 1.0),
        };

        let mut bh = BinaryHeap::new();
        bh.push(a);
        bh.push(b);
        bh.push(c);
        bh.push(d);
        let vec = bh.into_sorted_vec();
        let check = vec![
            Region {
                error: 1.533,
                result: 1.0,
                result_abs: 1.0,
                result_asc: 1.0,
                limits: Limits::new(0.0, 1.0),
            },
            Region {
                error: 1.533,
                result: 1.0,
                result_abs: 1.0,
                result_asc: 1.0,
                limits: Limits::new(0.5, 1.0),
            },
            Region {
                error: 1.60,
                result: 1.0,
                result_abs: 1.0,
                result_asc: 1.0,
                limits: Limits::new(0.0, 1.0),
            },
            Region {
                error: 2.0,
                result: 1.0,
                result_abs: 1.0,
                result_asc: 1.0,
                limits: Limits::new(0.0, 1.0),
            },
        ];

        assert_eq!(vec, check);
    }
}
