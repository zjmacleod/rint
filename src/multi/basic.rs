use crate::multi::{Integrator, Rule};
use crate::IntegralEstimate;
use crate::Limits;
use crate::MultiDimensionalIntegrand;
use crate::{InitialisationError, InitialisationErrorKind};

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Basic<'a, I, R, const N: usize> {
    function: &'a I,
    rule: &'a R,
    limits: [Limits; N],
}

impl<'a, I, const N: usize, const FINAL: usize, const TOTAL: usize>
    Basic<'a, I, Rule<N, FINAL, TOTAL>, N>
where
    I: MultiDimensionalIntegrand<N>,
{
    /// Create a new [`Basic`] multi-dimensional integrator.
    ///
    /// The user first defines a `function` which is a `struct` implementing the
    /// [`MultiDimensionalIntegrand`] trait and selects a fully-symmetric multi-dimensional
    /// integration [`Rule`], `rule`, to integrate the function in the hypercube formed by the
    /// [`Limits`], `limits` in each of the `N` integration directions.
    ///
    /// # Errors
    /// Will fail if `N < 2` or `N > 15`. The routines probided in this module are developed
    /// for dimensionalities between `2 <= N <= 15`.
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

    #[must_use]
    pub fn integrate(&self) -> IntegralEstimate<I::Scalar> {
        let integral = self.integrator().integrate();
        IntegralEstimate::new()
            .with_result(integral.result())
            .with_error(integral.error())
            .with_iterations(1)
            .with_evaluations(self.rule.evaluations())
    }

    const fn integrator(&self) -> Integrator<'_, I, N, FINAL, TOTAL> {
        Integrator::new(self.function, self.rule, self.limits)
    }
}
