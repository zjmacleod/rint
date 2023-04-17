use crate::integration::rescale_error;
use crate::rule::Rule;
use crate::Integrand;

/// The value of a function evaluated with Gauss-Kronrod integration and associated error
/// estimation.
pub struct GaussKronrod {
    result: f64,
    result_abs: f64,
    result_asc: f64,
    error: f64,
}

impl GaussKronrod {
    /// Return the numerically approximated value of the integral.
    #[must_use]
    pub fn result(&self) -> f64 {
        self.result
    }

    /// Return the numerically approximated value of the absolute value of the integral.
    #[must_use]
    pub fn result_abs(&self) -> f64 {
        self.result_abs
    }

    #[must_use]
    pub fn result_asc(&self) -> f64 {
        self.result_asc
    }

    /// Return the numerically approximated error.
    #[must_use]
    pub fn error(&self) -> f64 {
        self.error
    }
}

/// An integral to be evaluated with Gauss-Kronrod quadrature.
///
/// The user defines an [`Integrand`] to be integrated using a Gauss-Kronrod
/// integration [`Rule`] between the two integration limits `lower` and `upper`.
/// By convention `upper > lower`, however this is not mandatory.
pub struct GaussKronrodBasic<I, R>
where
    I: Integrand,
    R: Rule,
{
    lower: f64,
    upper: f64,
    rule: R,
    function: I,
}

impl<I, R> GaussKronrodBasic<I, R>
where
    I: Integrand,
    R: Rule,
{
    /// Create a new [`GaussKronrodBasic`].
    ///
    /// The user first defines a `function` which is a `struct` implementing the
    /// [`Integrand`] trait and selects a Gauss-Kronrod quadrature [`Rule`],
    /// `rule`, to integrate the function with between the integration limis,
    /// `upper` and `lower`.
    pub fn new(lower: f64, upper: f64, rule: R, function: I) -> Self {
        Self {
            lower,
            upper,
            rule,
            function,
        }
    }

    /// Integrate the function and return a [`GaussKronrod`] integration result.
    pub fn integrate(&self) -> GaussKronrod {
        let center = 0.5 * (self.upper + self.lower);
        let half_length = 0.5 * (self.upper - self.lower);
        let abs_half_length = half_length.abs();

        let initial_kronrod = 0.0;
        let initial_gauss = 0.0;
        let initial_abs = 0.0;
        let initial_asc = 0.0;

        let gauss_result = self
            .rule
            .gauss_nodes()
            .into_iter()
            .zip(self.rule.gauss_weights().into_iter())
            .map(|(t, w)| {
                let point = (half_length * t) + center;
                let rate = self.function.evaluate(&point);
                w * rate
            })
            .fold(initial_gauss, |a, v| a + v);

        let mut function_values: Vec<f64> = Vec::with_capacity(61);

        let (kronrod_result, abs_result) = self
            .rule
            .kronrod_nodes()
            .into_iter()
            .zip(self.rule.kronrod_weights().into_iter())
            .map(|(t, w)| {
                let point = (half_length * t) + center;
                let rate = self.function.evaluate(&point);
                let rate_abs = rate.abs();

                function_values.push(rate);
                (w * rate, w * rate_abs)
            })
            .fold((initial_kronrod, initial_abs), |a, v| {
                (a.0 + v.0, a.1 + v.1)
            });

        let mean = kronrod_result * 0.5;

        let asc_result = function_values
            .iter()
            .zip(self.rule.kronrod_weights().into_iter())
            .map(|(f, w)| w * (f - mean).abs())
            .fold(initial_asc, |a, v| a + v);

        let error = (kronrod_result - gauss_result) * half_length;

        let result = kronrod_result * half_length;
        let result_abs = abs_result * abs_half_length;
        let result_asc = asc_result * abs_half_length;
        let error = rescale_error(error, result_abs, result_asc);

        GaussKronrod {
            result,
            result_abs,
            result_asc,
            error,
        }
    }

    /// Return the value of the `upper` integration limit.
    pub fn upper(&self) -> f64 {
        self.upper
    }

    /// Return the value of the `lower` integration limit.
    pub fn lower(&self) -> f64 {
        self.lower
    }
}
