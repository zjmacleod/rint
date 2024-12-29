use num_complex::ComplexFloat;
use num_traits::Zero;

use crate::quadrature::region::Region;
use crate::quadrature::rescale_error;
use crate::quadrature::rule::Rule;
use crate::Integrand;
use crate::Limits;

pub(crate) struct Integrator<'a, I>
where
    I: Integrand,
{
    function: &'a I,
    rule: &'a Rule,
    limits: Limits,
}

impl<'a, I> Integrator<'a, I>
where
    I: Integrand,
{
    pub(crate) const fn new(function: &'a I, rule: &'a Rule, limits: Limits) -> Self {
        Self {
            function,
            rule,
            limits,
        }
    }

    pub(crate) fn integrate(&self) -> Region<I> {
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
                (rate * gauss, rate * kronrod, kronrod * rate_abs)
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
