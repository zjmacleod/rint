use crate::multi::{Integrator, Rule};
use crate::IntegralEstimate;
use crate::Limits;
use crate::MultiDimensionalIntegrand;
use crate::{InitialisationError, InitialisationErrorKind};

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Basic<'a, I, R, const NDIM: usize> {
    function: &'a I,
    rule: &'a R,
    limits: [Limits; NDIM],
}

impl<'a, I, const NDIM: usize, const FINAL: usize, const TOTAL: usize>
    Basic<'a, I, Rule<NDIM, FINAL, TOTAL>, NDIM>
where
    I: MultiDimensionalIntegrand<NDIM>,
{
    /// Create a new [`Basic`] multi-dimensional integrator.
    ///
    /// The user first defines a `function` which is a `struct` implementing the
    /// [`MultiDimensionalIntegrand`] trait and selects a fully-symmetric multi-dimensional
    /// integration [`Rule`], `rule`, to integrate the function in the hypercube formed by the
    /// [`Limits`], `limits` in each of the `NDIM` integration directions.
    ///
    /// # Errors
    /// Will fail if `NDIM < 2` or `NDIM > 15`. The routines probided in this module are developed
    /// for dimensionalities between `2 <= NDIM <= 15`.
    pub fn new(
        function: &'a I,
        rule: &'a Rule<NDIM, FINAL, TOTAL>,
        limits: [Limits; NDIM],
    ) -> Result<Self, InitialisationError> {
        if NDIM < 2 || NDIM > 15 {
            return Err(InitialisationError::new(
                InitialisationErrorKind::InvalidDimension(NDIM),
            ));
        };
        Ok(Self {
            function,
            rule,
            limits,
        })
    }

    pub fn integrate(&self) -> IntegralEstimate<I::Scalar> {
        let integral = self.integrator().integrate();
        IntegralEstimate::new()
            .with_result(integral.result())
            .with_error(integral.error())
            .with_iterations(1)
            .with_evaluations(self.rule.evaluations())
    }

    const fn integrator(&self) -> Integrator<'_, I, NDIM, FINAL, TOTAL> {
        Integrator::new(&self.function, &self.rule, self.limits)
    }
}
