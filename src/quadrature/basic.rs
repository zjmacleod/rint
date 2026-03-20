use crate::quadrature::{Integrator, Rule};
use crate::IntegralEstimate;
use crate::Integrand;
use crate::Limits;

/// A non-adaptive Gauss-Kronrod integrator.
///
/// The [`Basic`] integrator applies a Gauss-Kronrod integration [`Rule`] to approximate the
/// integral of a one-dimensional function. It is non-adaptive: it runs exactly once on the input
/// function. Thus, it is only suitable for the integration of smooth functions with no problematic
/// regions in the integration region. If higher accuracy is required then the [`Adaptive`] or
/// [`AdaptiveSingularity`] integrators should be used instead.
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
///     type Point = f64;
///     type Scalar = f64;
///
///     fn evaluate(&self, x: Self::Point) -> Self::Scalar {
///         1.0 / (1.0 + x.powi(2)).sqrt()
///     }
/// }
///
/// # use std::error::Error;
/// # fn main() -> Result<(), Box<dyn Error>> {
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
/// # Ok(())
/// # }
///```
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Basic<'a, I> {
    function: &'a I,
    rule: &'a Rule,
    limits: Limits,
}

impl<'a, I> Basic<'a, I>
where
    I: Integrand<Point = f64>,
{
    /// Create a new [`Basic`] integrator.
    ///
    /// The user first defines a `function` which is something implementing the [`Integrand`]
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
