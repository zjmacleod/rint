use crate::quadrature::{Integrator, Rule};
use crate::IntegralEstimate;
use crate::Integrand;
use crate::Limits;

/// A non-adaptive Gauss-Kronrod integrator.
///
/// The [`Basic`] integrator is a non-adaptive integrator designed for approximating
/// one-dimensional integrals of the form,
/// $$
/// I = \int_{b}^{a} f(x) dx
/// $$
/// using Gauss-Kronrod integration rules. The function $f(x)$ is encoded in [`Basic`] as something
/// implementing the [`Integrand`] trait. A Gaussian numerical integration rule approximates an
/// integral of a function by performing a weighted sum of the function evaluated at defined
/// points/abscissae. The order of an integration rule, $n$, denotes the number of abscissae,
/// $x_{i}$, at which the function is evaluated and the number of weights $w_{i}$ for the weighted
/// sum, such that the approximation is,
/// $$
/// I \approx \sum_{i = 1}^{n} W_{i} f(X_{i}) = I_{n}
/// $$
/// where the $X_{i}$ and $W_{i}$ are the rescaled abscissae and weights,
/// $$
/// X_{i} = \frac{b + a + (a - b) x_{i}}{2} ~~~~~~~~ W_{i} = \frac{(a - b) w_{i}}{2}
/// $$
/// where the `limits` $a$ and $b$ are encoded in [`Basic`] via [`Limits`] passed to the
/// constructor. A Gauss-Kronrod integration rule combines two rules of different order for
/// efficient estimation of the numerical error. The rules for an $n$-point Gauss-Kronrod rule
/// contain $m = (n - 1) / 2$ abscissae _shared_ by the Gaussian and Kronrod rules and an extended
/// set of $n - m$ Kronrod abscissae. The weighted sum of the full set of $n$ Kronrod function
/// evaluations are used to approximate the result of the integration, while the weighted sum of
/// the lower order set of $m$ Gaussian points are used to calculate the numerical error in the
/// routine,
/// $$
/// E = |I_{n} - I_{m}|
/// $$
/// This approach is efficient, as only $n$ total function evaluations are required to obtain the
/// result approximation and error estimate. See [`Rule`] for the available Gauss-Kronrod rules.
///
/// This routine is *not* adaptive, and runs exactly once on the input function. Thus, it is only
/// suitable for the integration of smooth functions with no problematic regions in the integration
/// region. If higher accuracy is required then the [`Adaptive`] or [`AdaptiveSingularity`]
/// integrators should be used instead.
///
/// [`Adaptive`]: crate::quadrature::Adaptive
/// [`AdaptiveSingularity`]: crate::quadrature::AdaptiveSingularity
///
/// # Example
///
/// Here we present a calculation of the golden ratio $\varphi$ ([`std::f64::consts::GOLDEN_RATIO`])
/// using the integral representation,
/// $$
/// \ln \varphi = \int_{0}^{1/2} \frac{dx}{\sqrt{1 + x^{2}}}
/// $$
///```rust
/// use rint::{Integrand, Limits};
/// use rint::quadrature::{Basic, Rule};
///
/// const PHI: f64 = std::f64::consts::GOLDEN_RATIO;
///
/// struct GoldenRatio;
///
/// impl Integrand for GoldenRatio {
///     type Scalar = f64;
///
///     fn evaluate(&self, x: f64) -> Self::Scalar {
///         1.0 / (1.0 + x.powi(2)).sqrt()
///     }
/// }
///
/// let golden_ratio = GoldenRatio;
/// let limits = Limits::new(0.0,0.5);
/// let rule = Rule::gk15();
/// let integral = Basic::new(&golden_ratio, &rule, limits)
///     .integrate();
///
/// let result = integral.result();
/// let error = integral.error();
/// let abs_actual_error = (PHI.ln() - result).abs();
/// let iters = integral.iterations();
/// assert_eq!(iters, 1);
/// assert!(abs_actual_error < error);
///```
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Basic<'a, I> {
    function: &'a I,
    rule: &'a Rule,
    limits: Limits,
}

impl<'a, I> Basic<'a, I>
where
    I: Integrand,
{
    /// Create a new [`Basic`] integrator.
    ///
    /// The user first defines a `function` which is a `struct` implementing the [`Integrand`]
    /// trait and selects a Gauss-Kronrod quadrature [`Rule`], `rule`, to integrate the function
    /// with between the integration [`Limits`], `limits`.
    pub const fn new(function: &'a I, rule: &'a Rule, limits: Limits) -> Self {
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
    /// function evaluations and integration routine iterations. Note: for the [`Basic`] integrator
    /// the number of iterations is 1, and the number of function evaluations for an $n$-point
    /// integration rule is $n$.
    #[must_use]
    pub fn integrate(&self) -> IntegralEstimate<I::Scalar> {
        let integral = self.integrator().integrate();
        IntegralEstimate::new()
            .with_result(integral.result())
            .with_error(integral.error())
            .with_iterations(1)
            .with_evaluations(self.rule.evaluations())
    }

    /// Return the integration [`Limits`].
    #[must_use]
    pub const fn limits(&self) -> Limits {
        self.limits
    }

    /// Create an [`Integrator`] from the integration data stored in a [`Basic`].
    #[must_use]
    const fn integrator(&self) -> Integrator<'_, I> {
        Integrator::new(self.function, self.rule, self.limits)
    }
}
