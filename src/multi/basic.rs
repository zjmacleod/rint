use crate::multi::{Integrator, Rule};
use crate::IntegralEstimate;
use crate::Integrand;
use crate::Limits;
use crate::{InitialisationError, InitialisationErrorKind};

/// A non-adaptive multi-dimensional integrator.
///
/// The [`Basic`] integrator applies a fully-symmetric integration [`Rule`] to approximate the
/// integral of an $N$-dimensional function. It is non-adaptive: it runs exactly once on the input
/// function. Thus, it is only suitable for the integration of smooth functions with no problematic
/// regions in the integration region. If higher accuracy is required then the [`Adaptive`].
///
/// [`Adaptive`]: crate::multi::Adaptive
///
/// # Example
///
/// Here we present a calculation of [Catalan's constant] $G$ using the integral representation:
/// $$
/// G = \int_{0}^{1} \int_{0}^{1} \frac{1}{1 + x^{2} y^{2}} dy dx
/// $$
///
/// [Catalan's constant]: <https://en.wikipedia.org/wiki/Catalan%27s_constant>
///```rust
/// use rint::{Integrand, Limits};
/// use rint::multi::{Basic, Rule13};
///
/// const G: f64 = 0.915_965_594_177_219_015_054_603_514_932_384_110_774;
/// const N: usize = 2;
///
/// struct Catalan;
///
/// impl Integrand for Catalan {
///     type Point = [f64; N];
///     type Scalar = f64;
///
///     fn evaluate(&self, coordinate: &[f64; N]) -> Self::Scalar {
///         let x = coordinate[0];
///         let y = coordinate[1];
///
///         1.0 / (1.0 + x.powi(2) * y.powi(2))
///     }
/// }
///
/// # use std::error::Error;
/// # fn main() -> Result<(), Box<dyn Error>> {
/// let catalan = Catalan;
/// let limits = [Limits::new(0.0,1.0),Limits::new(0.0,1.0)];
/// let rule = Rule13::generate();
/// let integral = Basic::new(&catalan, &rule, limits)?.integrate();
///
/// let result = integral.result();
/// let error = integral.error();
/// let abs_actual_error = (G - result).abs();
/// let iters = integral.iterations();
/// assert_eq!(iters, 1);
/// assert!(abs_actual_error < error);
/// # Ok(())
/// # }
///```
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Basic<'a, I, R, const N: usize> {
    function: &'a I,
    rule: &'a R,
    limits: [Limits; N],
}

impl<'a, I, const N: usize, const FINAL: usize, const TOTAL: usize>
    Basic<'a, I, Rule<N, FINAL, TOTAL>, N>
where
    I: Integrand<Point = [f64; N]>,
{
    /// Create a new [`Basic`] multi-dimensional integrator.
    ///
    /// The user first defines a `function` which is something implementing the [`Integrand`] trait
    /// and selects a fully-symmetric multi-dimensional integration [`Rule`], `rule`, to integrate
    /// the function in the hypercube formed by the [`Limits`], `limits` in each of the $N$
    /// integration directions.
    ///
    /// # Errors
    ///
    /// Will fail if $N < 2$ or $N > 15$. The routines probided in this module are developed
    /// for dimensionalities between $2 \le N \le 15$.
    pub fn new(
        function: &'a I,
        rule: &'a Rule<N, FINAL, TOTAL>,
        limits: [Limits; N],
    ) -> Result<Self, InitialisationError> {
        if N < 2 || N > 15 {
            return Err(InitialisationError::new(
                InitialisationErrorKind::InvalidDimension(N),
            ));
        }
        Ok(Self {
            function,
            rule,
            limits,
        })
    }

    /// Integrate the given function.
    ///
    /// Applies the user supplied integration rule to obtain an [`IntegralEstimate`], which is the
    /// numerically evaluated estimate of the integral value and error, as well as the number of
    /// function evaluations and integration routine iterations. Note: for the [`Basic`] integrator
    /// the number of iterations is 1.
    #[must_use]
    pub fn integrate(&self) -> IntegralEstimate<I::Scalar> {
        let integral = self.integrator().integrate();
        IntegralEstimate::new()
            .with_result(integral.result())
            .with_error(integral.error())
            .with_iterations(1)
            .with_evaluations(self.rule.evaluations())
    }

    /// Create an [`Integrator`] from the integration data stored in a [`Basic`].
    const fn integrator(&self) -> Integrator<'_, I, N, FINAL, TOTAL> {
        Integrator::new(self.function, self.rule, self.limits)
    }
}
