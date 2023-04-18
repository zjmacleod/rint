//! A basic Gauss-Kronrod integration routine.
//!
//! This is the most basic integration routine provided by this library. It applies an n-point
//! Gauss-Kronrod integration [`Rule`] to a user supplied function, which is a struct implementing
//! the [`Integrand`] trait.
//! The integration routine tries to be as efficient as possible by reusing function
//! evaluations of the lower n-point Gauss integration in the calculation of the higher-order
//! Kronrod extension. The available Gauss-Kronrod integration rules are:
//! * [`GaussKronrod15`] - 7-point Gauss, 15-point Kronrod, 8 total function evaluations
//! * [`GaussKronrod21`] - 10-point Gauss, 21-point Kronrod, 11 total function evaluations
//! * [`GaussKronrod31`] - 15-point Gauss, 31-point Kronrod, 16 total function evaluations
//! * [`GaussKronrod41`] - 20-point Gauss, 41-point Kronrod, 21 total function evaluations
//! * [`GaussKronrod51`] -  25-point Gauss, 51-point Kronrod, 26 total function evaluations
//! * [`GaussKronrod61`] -  30-point Gauss, 61-point Kronrod, 31 total function evaluations
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
//! use gauss_kronrod_adaptive_integration::Integrand;
//! use gauss_kronrod_adaptive_integration::integration::GaussKronrodBasic;
//! use gauss_kronrod_adaptive_integration::rule::GaussKronrod15;
//!
//! /* f1(x) = x^alpha * log(1/x) */
//! /* integ(f1,x,0,1) = 1/(alpha + 1)^2 */
//! pub struct Function1 {
//!     pub alpha: f64,
//! }
//!
//! impl Integrand for Function1 {
//!     fn evaluate(&self, x: &f64) -> f64 {
//!         let alpha = self.alpha;
//!         x.powf(alpha) * (1.0 / x).ln()
//!     }
//! }
//!
//! let exp_result: f64 = 7.716049357767090777E-02;
//! let exp_abserr: f64 = 2.990224871000550874E-06;
//! let exp_resabs: f64 = 7.716049357767090777E-02;
//! let exp_resasc: f64 = 4.434273814139995384E-02;
//!
//! let rel_error_bound = 1.0e-15;
//! let alpha = 2.6;
//!
//! let lower = 0.0;
//! let upper = 1.0;
//!
//! let function = Function1 { alpha };
//! let rule = GaussKronrod15;
//! let integral = GaussKronrodBasic::new(lower, upper, rule, function);
//!
//! let integral_result = integral.integrate();
//! let result = integral_result.result();
//! let abserr = integral_result.error();
//! let resabs = integral_result.result_abs();
//! let resasc = integral_result.result_asc();
//!
//! assert!((exp_result - result).abs() / exp_result.abs() < rel_error_bound);
//! assert!((exp_abserr - abserr).abs() / exp_abserr.abs() < rel_error_bound);
//! assert!((exp_resabs - resabs).abs() / exp_resabs.abs() < rel_error_bound);
//! assert!((exp_resasc - resasc).abs() / exp_resasc.abs() < rel_error_bound);
//!```
use crate::integration::rescale_error;
use crate::rule::Rule;
use crate::Integrand;

/// The value of a function evaluated with Gauss-Kronrod integration and associated error
/// estimation.
pub struct Basic {
    result: f64,
    result_abs: f64,
    result_asc: f64,
    error: f64,
}

impl Basic {
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

/// A basic Gauss-Kronrod integration routine.
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

    /// Integrate `function`, returning a [`Basic`] integration result.
    pub fn integrate(&self) -> Basic {
        let centre = 0.5 * (self.upper + self.lower);
        let half_length = 0.5 * (self.upper - self.lower);
        let abs_half_length = half_length.abs();
        let f_centre = self.function.evaluate(&centre);

        let initial_kronrod = self.rule.kronrod_centre() * f_centre;
        let initial_gauss = if let Some(v) = self.rule.gauss_centre() {
            v * f_centre
        } else {
            0.0
        };
        let initial_abs = initial_kronrod.abs();

        // XXX Can we get rid of this additional allocation?
        let mut function_values: Vec<(f64, (f64, f64))> = Vec::with_capacity(61);

        let (gauss_result, kronrod_shared, abs_shared) = self
            .rule
            .shared_nodes()
            .into_iter()
            .zip(
                self.rule
                    .gauss_weights()
                    .into_iter()
                    .zip(self.rule.kronrod_weights().into_iter()),
            )
            .map(|(t, (g, k))| {
                let abscissa = half_length * t;
                let rate_plus = self.function.evaluate(&(centre + abscissa));
                let rate_minus = self.function.evaluate(&(centre - abscissa));
                let rate = rate_plus + rate_minus;
                let rate_abs = rate_plus.abs() + rate_minus.abs();
                function_values.push((k, (rate_plus, rate_minus)));
                (g * rate, k * rate, k * rate_abs)
            })
            .fold((initial_gauss, initial_kronrod, initial_abs), |a, v| {
                (a.0 + v.0, a.1 + v.1, a.2 + v.2)
            });

        let (kronrod_result, abs_result) = self
            .rule
            .extended_nodes()
            .into_iter()
            .zip(self.rule.extended_weights().into_iter())
            .map(|(t, k)| {
                let abscissa = half_length * t;
                let rate_plus = self.function.evaluate(&(centre + abscissa));
                let rate_minus = self.function.evaluate(&(centre - abscissa));
                let rate = rate_plus + rate_minus;
                let rate_abs = rate_plus.abs() + rate_minus.abs();
                function_values.push((k, (rate_plus, rate_minus)));
                (k * rate, k * rate_abs)
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

        Basic {
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
