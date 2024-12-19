use std::cmp::Ordering;

use crate::quadrature::integrator::Integrator;
use crate::quadrature::rule::Rule;
use crate::Integrand;
use crate::Limits;

#[derive(Debug, PartialEq)]
pub(crate) struct Region {
    pub(crate) result: f64,
    pub(crate) error: f64,
    pub(crate) result_abs: f64,
    pub(crate) result_asc: f64,
    pub(crate) limits: Limits,
}

impl Eq for Region {}

impl PartialOrd for Region {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Region {
    fn cmp(&self, other: &Self) -> Ordering {
        let mut ordering = self.total_cmp_error(other);
        if let Ordering::Equal = ordering {
            ordering = self.total_cmp_interval_length(other);
        }
        ordering
    }
}

impl Region {
    pub(crate) fn unevaluated() -> Self {
        Self {
            error: 0.0,
            result: 0.0,
            result_abs: 0.0,
            result_asc: 0.0,
            limits: Limits::new(0.0, 0.0),
        }
    }

    pub(crate) fn with_result(mut self, result: f64) -> Self {
        self.result = result;
        self
    }

    pub(crate) fn with_result_asc(mut self, result_asc: f64) -> Self {
        self.result_asc = result_asc;
        self
    }

    pub(crate) fn with_result_abs(mut self, result_abs: f64) -> Self {
        self.result_abs = result_abs;
        self
    }

    pub(crate) fn with_error(mut self, error: f64) -> Self {
        self.error = error;
        self
    }

    pub(crate) fn with_limits(mut self, limits: Limits) -> Self {
        self.limits = limits;
        self
    }

    #[must_use]
    #[inline]
    pub(crate) fn result(&self) -> f64 {
        self.result
    }

    #[must_use]
    #[inline]
    pub(crate) fn result_asc(&self) -> f64 {
        self.result_asc
    }

    #[must_use]
    #[inline]
    pub(crate) fn result_abs(&self) -> f64 {
        self.result_abs
    }

    #[must_use]
    #[inline]
    pub(crate) fn error(&self) -> f64 {
        self.error
    }

    #[must_use]
    pub(crate) fn limits(&self) -> Limits {
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

    pub(crate) fn bisect<I: Integrand>(&self, function: &I, rule: &Rule) -> [Region; 2] {
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ordering_of_basic_internal() {
        use std::collections::BinaryHeap;

        let a = Region {
            error: 2.0,
            result: 1.0,
            result_abs: 1.0,
            result_asc: 1.0,
            limits: Limits::new(0.0, 1.0),
        };
        let b = Region {
            error: 1.533,
            result: 1.0,
            result_abs: 1.0,
            result_asc: 1.0,
            limits: Limits::new(0.0, 1.0),
        };
        let c = Region {
            error: 1.533,
            result: 1.0,
            result_abs: 1.0,
            result_asc: 1.0,
            limits: Limits::new(0.5, 1.0),
        };
        let d = Region {
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
