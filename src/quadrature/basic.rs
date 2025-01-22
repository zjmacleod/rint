use crate::quadrature::{Integrator, Rule};
use crate::IntegralEstimate;
use crate::Integrand;
use crate::Limits;

/// A non-adaptive Gauss-Kronrod integrator.
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
///     type Scalar = f64;
///     fn evaluate(&self, x: f64) -> Self::Scalar {
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
pub struct Basic<I> {
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
            function,
            rule,
            limits,
        }
    }

    /// Integrate the given function.
    ///
    /// Applies the user supplied integration rule to obtain an [`IntegralEstimate`], which is the
    /// numerically evaluated estimate of the integral value and error, as well as the number of
    /// function evaluations and integration routine iterations.
    /// Note: for the [`Basic`] integrator the number of iterations is `1`, and the number of
    /// function evaluations for an `n`-point integration rule is `n`.
    pub fn integrate(&self) -> IntegralEstimate<I::Scalar> {
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

    const fn integrator(&self) -> Integrator<'_, I> {
        Integrator::new(&self.function, &self.rule, self.limits)
    }
}
