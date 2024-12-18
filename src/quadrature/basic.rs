use std::cmp::Ordering;

use crate::quadrature::rescale_error;
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
/// evaluations of the lower `n`-point Gauss integration in the calculation of the higher-order
/// Kronrod extension. The total number of function evaluations is `n`.
/// The available Gauss-Kronrod integration rules are:
/// * 15-point: 7-point Gauss, 15-point Kronrod ([`Rule::gk15()`])
/// * 21-point: 10-point Gauss, 21-point Kronrod ([`Rule::gk21()`])
/// * 31-point: 15-point Gauss, 31-point Kronrod ([`Rule::gk31()`])
/// * 41-point: 20-point Gauss, 41-point Kronrod ([`Rule::gk41()`])
/// * 51-point: 25-point Gauss, 51-point Kronrod ([`Rule::gk51()`])
/// * 61-point: 30-point Gauss, 61-point Kronrod ([`Rule::gk61()`])
///
/// The user defines an [`Integrand`] to be integrated using a fixed Gauss-Kronrod
/// integration [`Rule`] between the two integration limits `lower` and `upper`.
/// By convention `upper > lower`, however this is not mandatory.
///
///```rust
/// use rint::{Limits, Integrand};
/// use rint::quadrature::Basic;
/// use rint::quadrature::Rule;
///
/// /* f1(x) = x^alpha * log(1/x) */
/// /* integ(f1,x,0,1) = 1/(alpha + 1)^2 */
/// struct Function1 {
///     alpha: f64,
/// }
///
/// impl Integrand for Function1 {
///     fn evaluate(&self, x: f64) -> f64 {
///         let alpha = self.alpha;
///         x.powf(alpha) * (1.0 / x).ln()
///     }
/// }
///
/// let exp_result: f64 = 7.716049357767090777E-02;
/// let exp_abserr: f64 = 2.990224871000550874E-06;
///
/// let rel_error_bound = 1.0e-15;
/// let alpha = 2.6;
///
/// let limits = Limits::new(0.0, 1.0);
///
/// let function = Function1 { alpha };
/// let rule = Rule::gk15();
/// let integral = Basic::new(&function, rule, limits);
///
/// let integral_result = integral.integrate();
/// let result = integral_result.result();
/// let abserr = integral_result.error();
///
/// println!("{}",(exp_result - result).abs() / exp_result.abs());
/// println!("{}",(exp_abserr - abserr).abs() / exp_abserr.abs());
///
/// assert!((exp_result - result).abs() / exp_result.abs() < rel_error_bound);
/// assert!((exp_abserr - abserr).abs() / exp_abserr.abs() < rel_error_bound);
///```
pub struct Basic<I>
where
    I: Integrand,
{
    function: I,
    rule: Rule,
    limits: Limits,
}

impl<I> Basic<I>
where
    I: Integrand,
{
    /// Create a new [`Basic`] integrator.
    ///
    /// The user first defines a `function` which is a `struct` implementing the
    /// [`Integrand`] trait and selects a Gauss-Kronrod quadrature [`Rule`],
    /// `rule`, to integrate the function with between the integration [`Limits`],
    /// `limits`.
    pub const fn new(function: I, rule: Rule, limits: Limits) -> Self {
        Self {
            limits,
            rule,
            function,
        }
    }

    /// Integrate the given function.
    ///
    /// Applies the user supplied integration rule to obtain an [`IntegralEstimate`], which is the
    /// numerically evaluated estimate of the integral value and error, as well as the number of
    /// function evaluations and integration routine iterations.
    /// Note: for the [`Basic`] integrator the number of iterations is `1`, and the number of
    /// function evaluations for an `n`-point integration rule is `n`.
    pub fn integrate(&self) -> IntegralEstimate {
        let integral = self.integrator().integrate();
        IntegralEstimate::new()
            .with_result(integral.result())
            .with_error(integral.error())
            .with_iterations(1)
            .with_function_evaluations(self.rule.evaluations())
    }

    /// Return the integration [`Limits`].
    pub const fn limits(&self) -> Limits {
        self.limits
    }

    const fn integrator(&self) -> GaussKronrodIntegrator<'_, I> {
        GaussKronrodIntegrator::new(&self.function, &self.rule, self.limits)
    }
}

pub(crate) struct GaussKronrodIntegrator<'a, I>
where
    I: Integrand,
{
    function: &'a I,
    rule: &'a Rule,
    limits: Limits,
}

impl<'a, I> GaussKronrodIntegrator<'a, I>
where
    I: Integrand,
{
    pub(crate) const fn new(function: &'a I, rule: &'a Rule, limits: Limits) -> Self {
        Self {
            limits,
            rule,
            function,
        }
    }

    pub(crate) fn integrate(&self) -> Region {
        let centre = self.limits.centre();
        let half_length = self.limits.half_width();
        let abs_half_length = half_length.abs();
        let f_centre = self.function.evaluate(centre);

        let initial_kronrod = self.rule.kronrod_centre() * f_centre;
        let initial_gauss = if let Some(v) = self.rule.gauss_centre() {
            v * f_centre
        } else {
            0.0
        };
        let initial_abs = initial_kronrod.abs();

        // XXX Can we get rid of this additional allocation?
        // Vec<(kronrod_weight, (rate_plus, rate_minus))>
        let mut function_values: Vec<(f64, (f64, f64))> = Vec::with_capacity(61);

        let (gauss_result, kronrod_shared, abs_shared) = self
            .rule
            .shared_data()
            .iter()
            .map(|data| {
                let point = data.point();
                let gauss = data.gauss();
                let kronrod = data.kronrod();
                let abscissa = half_length * point;
                let rate_plus = self.function.evaluate(centre + abscissa);
                let rate_minus = self.function.evaluate(centre - abscissa);
                let rate = rate_plus + rate_minus;
                let rate_abs = rate_plus.abs() + rate_minus.abs();
                function_values.push((kronrod, (rate_plus, rate_minus)));
                (gauss * rate, kronrod * rate, kronrod * rate_abs)
            })
            .fold((initial_gauss, initial_kronrod, initial_abs), |a, v| {
                (a.0 + v.0, a.1 + v.1, a.2 + v.2)
            });

        let (kronrod_result, abs_result) = self
            .rule
            .extended_data()
            .iter()
            .map(|data| {
                let point = data.point();
                let kronrod = data.kronrod();
                let abscissa = half_length * point;
                let rate_plus = self.function.evaluate(centre + abscissa);
                let rate_minus = self.function.evaluate(centre - abscissa);
                let rate = rate_plus + rate_minus;
                let rate_abs = rate_plus.abs() + rate_minus.abs();
                function_values.push((kronrod, (rate_plus, rate_minus)));
                (kronrod * rate, kronrod * rate_abs)
            })
            .fold((kronrod_shared, abs_shared), |a, v| (a.0 + v.0, a.1 + v.1));

        let mean = kronrod_result * 0.5;

        let initial_asc = self.rule.kronrod_centre() * (f_centre - mean).abs();

        let asc_result = function_values
            .iter()
            .map(|(k, (rp, rm))| k * ((rp - mean).abs() + (rm - mean).abs()))
            .fold(initial_asc, |a, v| a + v);

        let error = (kronrod_result - gauss_result) * half_length;

        let result = kronrod_result * half_length;
        let result_abs = abs_result * abs_half_length;
        let result_asc = asc_result * abs_half_length;
        let error = rescale_error(error, result_abs, result_asc);

        Region::unevaluated()
            .with_error(error)
            .with_result(result)
            .with_result_abs(result_abs)
            .with_result_asc(result_asc)
            .with_limits(self.limits)
    }
}

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
        let lower_integral = GaussKronrodIntegrator::new(&function, &rule, lower).integrate();

        let upper_integral = GaussKronrodIntegrator::new(&function, &rule, upper).integrate();

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
