use crate::quadrature::region::Region;
use crate::quadrature::rule::Rule;
use crate::quadrature::{Adaptive, AdaptiveSingularity, Basic, Error, ErrorBound};
use crate::{Integrand, Limits};

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
}

/// Integrators
impl IntegralEstimate {
    /// Integrate a function using a basic (non-adaptive) Gauss-Kronrod integration rule.
    pub fn basic<I: Integrand>(limits: Limits, rule: Rule, function: I) -> Self {
        Basic::new(function, rule, limits).integrate()
    }

    /// Integrate a function using an adaptive Gauss-Kronrod integration routine.
    pub fn adaptive<I: Integrand>(
        limits: Limits,
        error_bound: ErrorBound,
        rule: Rule,
        function: I,
        max_iterations: usize,
    ) -> Result<Self, Error> {
        Adaptive::new(function, rule, limits, error_bound, max_iterations)?.integrate()
    }

    /// Integrate a function with possible singularities in the integration region using an adaptive
    /// Gauss-Kronrod integration routine.
    pub fn adaptive_singularity<I: Integrand>(
        limits: Limits,
        error_bound: ErrorBound,
        function: I,
        max_iterations: usize,
    ) -> Result<Self, Error> {
        AdaptiveSingularity::general(function, limits, error_bound, max_iterations)?.integrate()
    }

    pub fn infinite<I: Integrand>(
        error_bound: ErrorBound,
        function: I,
        max_iterations: usize,
    ) -> Result<Self, Error> {
        AdaptiveSingularity::infinite(function, error_bound, max_iterations)?.integrate()
    }

    pub fn semi_infinite_positive<I: Integrand>(
        lower: f64,
        error_bound: ErrorBound,
        function: I,
        max_iterations: usize,
    ) -> Result<Self, Error> {
        AdaptiveSingularity::semi_infinite_positive(function, lower, error_bound, max_iterations)?
            .integrate()
    }

    pub fn semi_infinite_negative<I: Integrand>(
        upper: f64,
        error_bound: ErrorBound,
        function: I,
        max_iterations: usize,
    ) -> Result<Self, Error> {
        AdaptiveSingularity::semi_infinite_negative(function, upper, error_bound, max_iterations)?
            .integrate()
    }
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

    pub(crate) fn from_region(
        basic: &Region,
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
