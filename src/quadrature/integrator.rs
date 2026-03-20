use num_complex::ComplexFloat;
use num_traits::Zero;

use crate::quadrature::region::Region;
use crate::quadrature::rescale_error;
use crate::quadrature::rule::Rule;
use crate::Integrand;
use crate::Limits;

/// The core integrator for the one-dimensional routines.
///
/// This is the core integrator which takes a one-dimensional `function` implementing the
/// [`Integrator`] trait, a Gauss-Kronrod integration [`Rule`], and an integration region defined
/// by the passed [`Limits`]. The [`Basic`] integration routine simply applies this integrator once
/// to the function over the given integration region, while the [`Adaptive`] and
/// [`AdaptiveSingularity`] routines use bisection approaches to concentrate more integration
/// points in regions of higher numerical error, applying the [`Integrator`] to each subsequent
/// bisected region.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub(crate) struct Integrator<'a, I> {
    function: &'a I,
    rule: &'a Rule,
    limits: Limits,
}

impl<'a, I> Integrator<'a, I>
where
    I: Integrand<Point = f64>,
{
    pub(crate) const fn new(function: &'a I, rule: &'a Rule, limits: Limits) -> Self {
        Self {
            function,
            rule,
            limits,
        }
    }

    /// The core integrator implementation for the one-dimensional routines.
    ///
    /// This method is the primary implementation of the Gauss-Kronrod integration rules. A
    /// Gaussian numerical integration rule approximates an integral of a function by performing a
    /// weighted sum of the function evaluated at defined abscissae. The order of an integration
    /// rule, $n$, denotes the number of abscissae, $x_{i}$, at which the function is evaluated and
    /// the number of weights $w_{i}$ for the weighted sum. A Gauss-Kronrod integration rule
    /// combines two rules of different order for efficient estimation of the numerical error.The
    /// rules for an $n$-point Gauss-Kronrod rule contain $m = (n - 1) / 2$ abscissae _shared_ by
    /// the Gaussian and Kronrod rules and an extended set of $n - m$ Kronrod abscissae. The
    /// weighted sum of the full set of $n$ Kronrod function evaluations are used to approximate
    /// the result of the integration, while the weighted sum of the lower order set of $m$
    /// Gaussian points are used to calculate the numerical error in the routine. This approach is
    /// efficient, as only $n$ total function evaluations are required to obtain the result
    /// approximation and error estimate.
    pub(crate) fn integrate(&self) -> Region<I::Scalar> {
        let centre = self.limits.centre();
        let half_length = self.limits.half_width();
        let abs_half_length = half_length.abs();
        let f_centre = self.function.evaluate(centre);

        let initial_kronrod = f_centre * self.rule.kronrod_centre();
        let initial_gauss = if let Some(v) = self.rule.gauss_centre() {
            f_centre * v
        } else {
            <I::Scalar as Zero>::zero()
        };
        let initial_abs = initial_kronrod.abs();

        // XXX Can we get rid of this additional allocation?
        // Vec<(kronrod_weight, (rate_plus, rate_minus))>
        let mut function_values = Vec::with_capacity(61);

        // Calculate shared
        let (gauss_result, kronrod_shared, abs_shared) = self
            .rule
            .shared_data()
            .iter()
            .map(|data| {
                let point = data.point();
                let gauss = data.gauss();
                let kronrod = data.kronrod();
                let x_plus = point.mul_add(half_length, centre);
                let x_minus = point.mul_add(-half_length, centre);
                let rate_plus = self.function.evaluate(x_plus);
                let rate_minus = self.function.evaluate(x_minus);
                let rate = rate_plus + rate_minus;
                let rate_abs = rate_plus.abs() + rate_minus.abs();
                function_values.push((kronrod, (rate_plus, rate_minus)));
                (rate * gauss, rate * kronrod, kronrod * rate_abs)
            })
            .fold((initial_gauss, initial_kronrod, initial_abs), |a, v| {
                (a.0 + v.0, a.1 + v.1, a.2 + v.2)
            });

        // Calculate extended
        let (kronrod_result, abs_result) = self
            .rule
            .extended_data()
            .iter()
            .map(|data| {
                let point = data.point();
                let kronrod = data.kronrod();
                let x_plus = point.mul_add(half_length, centre);
                let x_minus = point.mul_add(-half_length, centre);
                let rate_plus = self.function.evaluate(x_plus);
                let rate_minus = self.function.evaluate(x_minus);
                let rate = rate_plus + rate_minus;
                let rate_abs = rate_plus.abs() + rate_minus.abs();
                function_values.push((kronrod, (rate_plus, rate_minus)));
                (rate * kronrod, kronrod * rate_abs)
            })
            .fold((kronrod_shared, abs_shared), |a, v| (a.0 + v.0, a.1 + v.1));

        let mean = kronrod_result * 0.5;

        let initial_asc = self.rule.kronrod_centre() * (f_centre - mean).abs();

        let asc_result = function_values
            .iter()
            .map(|(k, (rp, rm))| k * ((*rp - mean).abs() + (*rm - mean).abs()))
            .fold(initial_asc, |a, v| a + v);

        let error = (kronrod_result - gauss_result).abs() * abs_half_length;

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
