use crate::quadrature::rule::Rule;
use crate::quadrature::{Adaptive, Basic};
//use crate::quadrature::{Adaptive, AdaptiveSingularity, Basic};
use crate::{Error, Integrand, Limits, ScalarF64, Tolerance};

/// The value of a function evaluated with Gauss-Kronrod integration and associated error
/// estimation.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct IntegralEstimate<T> {
    result: T,
    error: f64,
    iterations: usize,
    evaluations: usize,
}

/// # Getters
impl<T: ScalarF64> IntegralEstimate<T> {
    /// Return the numerically approximated value of the integral.
    #[must_use]
    pub const fn result(&self) -> T {
        self.result
    }

    /// Return the numerically approximated error.
    #[must_use]
    pub const fn error(&self) -> f64 {
        self.error
    }

    /// Return the number of iterations used in integration.
    #[must_use]
    pub const fn iterations(&self) -> usize {
        self.iterations
    }

    /// Return the number of function evaluations used in the integration.
    #[must_use]
    pub const fn evaluations(&self) -> usize {
        self.evaluations
    }
}

impl<T: ScalarF64> IntegralEstimate<T> {
    pub(crate) fn new() -> Self {
        let result = T::zero();
        Self {
            result,
            error: 0.0,
            iterations: 0,
            evaluations: 0,
        }
    }

    pub(crate) fn with_result(mut self, result: T) -> Self {
        self.result = result;
        self
    }

    pub(crate) fn with_error(mut self, error: f64) -> Self {
        self.error = error;
        self
    }

    pub(crate) fn with_iterations(mut self, iterations: usize) -> Self {
        self.iterations = iterations;
        self
    }

    pub(crate) fn with_evaluations(mut self, evaluations: usize) -> Self {
        self.evaluations = evaluations;
        self
    }
}

/// # Integrators
impl<T: ScalarF64> IntegralEstimate<T> {
    /// Integrate a function using a basic (non-adaptive) Gauss-Kronrod integration rule, see
    /// [`Basic`] for details.
    ///
    /// Note: for the [`Basic`] integrator the number of iterations is `1`, and the number of
    /// function evaluations for an `n`-point integration rule is `n`.
    pub fn basic<I: Integrand<Point = f64>>(
        function: &I,
        rule: &Rule,
        limits: Limits,
    ) -> IntegralEstimate<I::Scalar> {
        Basic::new(function, rule, limits).integrate()
    }

    /// Integrate a function using an adaptive Gauss-Kronrod integration routine, see [`Adaptive`]
    /// for details.
    ///
    /// # Errors
    /// The error enum [`Error`] will return the error variant, which either corresponds to an
    /// [`InitialisationError`] or an [`IntegrationError`].
    ///
    /// [`InitialisationError`]: crate::InitialisationError
    /// [`IntegrationError`]: crate::IntegrationError
    pub fn adaptive<'a, I: Integrand<Point = f64>>(
        function: &'a I,
        rule: &'a Rule,
        limits: Limits,
        tolerance: Tolerance,
        max_iterations: usize,
    ) -> Result<IntegralEstimate<I::Scalar>, Error<I::Scalar>> {
        let res = Adaptive::new(function, rule, limits, tolerance, max_iterations)?.integrate()?;
        Ok(res)
    }

    //    /// Integrate the function with possible integrable singularities and a finite integration
    //    /// interval, see [`AdaptiveSingularity`] for details.
    //    ///
    //    /// # Errors
    //    /// The error enum [`Error`] will return the error variant, which either corresponds to an
    //    /// [`InitialisationError`] or an [`IntegrationError`].
    //    ///
    //    /// [`InitialisationError`]: crate::InitialisationError
    //    /// [`IntegrationError`]: crate::IntegrationError
    //    pub fn adaptive_singularity_finite<I: Integrand>(
    //        function: I,
    //        limits: Limits,
    //        tolerance: Tolerance,
    //        max_iterations: usize,
    //    ) -> Result<IntegralEstimate<I::Scalar>, Error<I::Scalar>> {
    //        let res = AdaptiveSingularity::finite(function, limits, tolerance, max_iterations)?
    //            .integrate()?;
    //        Ok(res)
    //    }
    //
    //    /// Integrate the function with possible integrable singularities and an infinite integration
    //    /// interval, see [`AdaptiveSingularity`] for details.
    //    ///
    //    /// # Errors
    //    /// The error enum [`Error`] will return the error variant, which either corresponds to an
    //    /// [`InitialisationError`] or an [`IntegrationError`].
    //    ///
    //    /// [`InitialisationError`]: crate::InitialisationError
    //    /// [`IntegrationError`]: crate::IntegrationError
    //    pub fn adaptive_singularity_infinite<I: Integrand>(
    //        function: I,
    //        tolerance: Tolerance,
    //        max_iterations: usize,
    //    ) -> Result<IntegralEstimate<I::Scalar>, Error<I::Scalar>> {
    //        let res =
    //            AdaptiveSingularity::infinite(function, tolerance, max_iterations)?.integrate()?;
    //        Ok(res)
    //    }
    //
    //    /// Integrate the function with possible integrable singularities and semi-infinite integration
    //    /// interval (b, Inf), see [`AdaptiveSingularity`] for details.
    //    ///
    //    /// # Errors
    //    /// The error enum [`Error`] will return the error variant, which either corresponds to an
    //    /// [`InitialisationError`] or an [`IntegrationError`].
    //    ///
    //    /// [`InitialisationError`]: crate::InitialisationError
    //    /// [`IntegrationError`]: crate::IntegrationError
    //    pub fn adaptive_singularity_semi_infinite_upper<I: Integrand>(
    //        function: I,
    //        lower: f64,
    //        tolerance: Tolerance,
    //        max_iterations: usize,
    //    ) -> Result<IntegralEstimate<I::Scalar>, Error<I::Scalar>> {
    //        let res =
    //            AdaptiveSingularity::semi_infinite_upper(function, lower, tolerance, max_iterations)?
    //                .integrate()?;
    //        Ok(res)
    //    }
    //
    //    /// Integrate the function with possible integrable singularities and semi-infinite integration
    //    /// interval (-Inf, a), see [`AdaptiveSingularity`] for details.
    //    ///
    //    /// # Errors
    //    /// The error enum [`Error`] will return the error variant, which either corresponds to an
    //    /// [`InitialisationError`] or an [`IntegrationError`].
    //    ///
    //    /// [`InitialisationError`]: crate::InitialisationError
    //    /// [`IntegrationError`]: crate::IntegrationError
    //    pub fn adaptive_singularity_semi_infinite_lower<I: Integrand>(
    //        function: I,
    //        upper: f64,
    //        tolerance: Tolerance,
    //        max_iterations: usize,
    //    ) -> Result<IntegralEstimate<I::Scalar>, Error<I::Scalar>> {
    //        let res =
    //            AdaptiveSingularity::semi_infinite_lower(function, upper, tolerance, max_iterations)?
    //                .integrate()?;
    //        Ok(res)
    //    }
}
