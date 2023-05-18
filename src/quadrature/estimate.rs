use crate::quadrature::basic::BasicInternal;

/// The value of a function evaluated with Gauss-Kronrod integration and associated error
/// estimation.
/// The value of a function evaluated with Gauss-Kronrod integration and associated error
/// estimation.
#[derive(Debug)]
pub struct IntegralEstimate {
    result: f64,
    error: f64,
    iterations: usize,
    function_evaluations: usize,
}

impl IntegralEstimate {
    pub(crate) fn new() -> Self {
        Self {
            result: 0.0,
            error: 0.0,
            iterations: 0,
            function_evaluations: 0,
        }
    }

    pub(crate) fn with_result(mut self, result: f64) -> Self {
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

    pub(crate) fn with_function_evaluations(mut self, function_evaluations: usize) -> Self {
        self.function_evaluations = function_evaluations;
        self
    }

    /// Return the numerically approximated value of the integral.
    #[must_use]
    pub fn result(&self) -> f64 {
        self.result
    }

    /// Return the numerically approximated error.
    #[must_use]
    pub fn error(&self) -> f64 {
        self.error
    }

    /// Return the number of iterations used in the adaptive integration routine.
    #[must_use]
    pub fn iterations(&self) -> usize {
        self.iterations
    }

    /// Return the number of function evaluations used in the adaptive integration routine.
    #[must_use]
    pub fn function_evaluations(&self) -> usize {
        self.function_evaluations
    }

    pub(crate) fn from_basic(
        basic: &BasicInternal,
        iterations: usize,
        function_evaluations: usize,
    ) -> Self {
        let result = basic.result();
        let error = basic.error();
        Self::new()
            .with_result(result)
            .with_error(error)
            .with_iterations(iterations)
            .with_function_evaluations(function_evaluations)
    }
}
