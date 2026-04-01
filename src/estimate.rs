use crate::ScalarF64;

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
