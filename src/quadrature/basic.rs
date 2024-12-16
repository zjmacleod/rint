//! A basic Gauss-Kronrod integration routine.
//!
//! This is the most basic integration routine provided by this library. It applies an n-point
//! Gauss-Kronrod integration [`Rule`] to a user supplied function, which is a struct implementing
//! the [`Integrand`] trait.
//! The integration routine tries to be as efficient as possible by reusing function
//! evaluations of the lower n-point Gauss integration in the calculation of the higher-order
//! Kronrod extension. The available Gauss-Kronrod integration rules are:
//! * [`GaussKronrod15`] - 7-point Gauss, 15-point Kronrod, 15 total function evaluations
//! * [`GaussKronrod21`] - 10-point Gauss, 21-point Kronrod, 21 total function evaluations
//! * [`GaussKronrod31`] - 15-point Gauss, 31-point Kronrod, 31 total function evaluations
//! * [`GaussKronrod41`] - 20-point Gauss, 41-point Kronrod, 41 total function evaluations
//! * [`GaussKronrod51`] -  25-point Gauss, 51-point Kronrod, 51 total function evaluations
//! * [`GaussKronrod61`] -  30-point Gauss, 61-point Kronrod, 61 total function evaluations
//!
//! [`GaussKronrod15`]: ../../rule/struct.GaussKronrod15.html
//! [`GaussKronrod21`]: ../../rule/struct.GaussKronrod21.html
//! [`GaussKronrod31`]: ../../rule/struct.GaussKronrod31.html
//! [`GaussKronrod41`]: ../../rule/struct.GaussKronrod41.html
//! [`GaussKronrod51`]: ../../rule/struct.GaussKronrod51.html
//! [`GaussKronrod61`]: ../../rule/struct.GaussKronrod61.html
//!
//! The user defines an [`Integrand`] to be integrated using a fixed Gauss-Kronrod
//! integration [`Rule`] between the two integration limits `lower` and `upper`.
//! By convention `upper > lower`, however this is not mandatory.
//!
//!```rust
//! use rint::{Limits, Integrand};
//! use rint::quadrature::Basic;
//! use rint::quadrature::Rule;
//!
//! /* f1(x) = x^alpha * log(1/x) */
//! /* integ(f1,x,0,1) = 1/(alpha + 1)^2 */
//! pub struct Function1 {
//!     pub alpha: f64,
//! }
//!
//! impl Integrand for Function1 {
//!     fn evaluate(&self, x: f64) -> f64 {
//!         let alpha = self.alpha;
//!         x.powf(alpha) * (1.0 / x).ln()
//!     }
//! }
//!
//! let exp_result: f64 = 7.716049357767090777E-02;
//! let exp_abserr: f64 = 2.990224871000550874E-06;
//!
//! let rel_error_bound = 1.0e-15;
//! let alpha = 2.6;
//!
//! let limits = Limits::new(0.0, 1.0);
//!
//! let function = Function1 { alpha };
//! let rule = Rule::GaussKronrod15;
//! let integral = Basic::new(limits, rule, &function);
//!
//! let integral_result = integral.integrate();
//! let result = integral_result.result();
//! let abserr = integral_result.error();
//!
//! assert!((exp_result - result).abs() / exp_result.abs() < rel_error_bound);
//! assert!((exp_abserr - abserr).abs() / exp_abserr.abs() < rel_error_bound);
//!```
use std::cmp::Ordering;

use crate::quadrature::rule::Rule;
use crate::quadrature::IntegralEstimate;
use crate::Integrand;
use crate::Limits;

/// A basic Gauss-Kronrod integration routine.
///
/// This is the most basic integration routine provided by this library. It applies an n-point
/// Gauss-Kronrod integration [`Rule`] to a user supplied function, which is a struct implementing
/// the [`Integrand`] trait.
/// The integration routine tries to be as efficient as possible by reusing function
/// evaluations of the lower n-point Gauss integration in the calculation of the higher-order
/// Kronrod extension.
pub struct Basic<I>
where
    I: Integrand,
{
    limits: Limits,
    rule: Rule,
    function: I,
}

impl<I> Basic<I>
where
    I: Integrand,
{
    /// Create a new [`Basic`].
    ///
    /// The user first defines a `function` which is a `struct` implementing the
    /// [`Integrand`] trait and selects a Gauss-Kronrod quadrature [`Rule`],
    /// `rule`, to integrate the function with between the integration limis,
    /// `limits`.
    pub fn new(limits: Limits, rule: Rule, function: I) -> Self {
        Self {
            limits,
            rule,
            function,
        }
    }

    /// Integrate `function`, returning a [`Basic`] integration result.
    pub(crate) fn integrate_internal(&self) -> Region {
        self.rule.integrate(&self.limits, &self.function)
    }

    pub fn integrate(&self) -> IntegralEstimate {
        let integral = self.integrate_internal();
        IntegralEstimate::new()
            .with_result(integral.result())
            .with_error(integral.error())
            .with_iterations(1)
            .with_function_evaluations(self.rule.evaluations())
    }

    /// Return the integration [`Limits`].
    pub fn limits(&self) -> Limits {
        self.limits
    }
}

#[derive(Debug, PartialEq)]
pub(crate) struct Region {
    pub(crate) error: f64,
    pub(crate) result: f64,
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

    pub(crate) fn bisect<I: Integrand>(&self, function: &I, rule: Rule) -> [Region; 2] {
        let [lower, upper] = self.limits.bisect();
        let lower_integral = Basic::new(lower, rule, function).integrate_internal();

        let upper_integral = Basic::new(upper, rule, function).integrate_internal();

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
