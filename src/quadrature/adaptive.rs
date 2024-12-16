use std::collections::binary_heap::BinaryHeap;

use crate::quadrature::basic::Region;
use crate::quadrature::rule::Rule;
use crate::quadrature::subinterval_too_small;
use crate::quadrature::{Basic, Error, ErrorBound, IntegralEstimate, Kind};
use crate::Integrand;
use crate::Limits;

/// An integral to be evaluated with an adaptive Gauss-Kronrod quadrature.
///
/// The user constructs a `function` implementing [`Integrand`], provides
/// integration [`Limits`], and provides an `error_bound`, which can be
/// [`ErrorBound::Absolute`] to work to a specified absolute error,
/// [`ErrorBound::Relative`] to work to a specified relative error,
/// or [`ErrorBound::Either`] to return a result as soon as _either_ the relative
/// or absolute error bound has been satisfied.
pub struct Adaptive<I>
where
    I: Integrand,
{
    limits: Limits,
    error_bound: ErrorBound,
    rule: Rule,
    function: I,
    max_iterations: usize,
}

impl<I> Adaptive<I>
where
    I: Integrand,
{
    /// Create a new [`Adaptive`].
    ///
    /// The user defines a `function` which is a `struct` implementing the
    /// [`Integrand`] trait, and integration [`Limits`].
    ///
    /// # Errors
    /// Function will return an error if the user provided `ErrorBound` does not satisfy the
    /// following constraints:
    ///     - `ErrorBound::Absolute(v) where v > 0.0`,
    ///     - `ErrorBound::Relative(v) where v > 50.0 * f64::EPSILON`,
    ///     - `ErrorBound::Either { absolute, relative } where absolute > 0.0 and relative > 50.0 * f64::EPSILON`.
    pub fn new(
        limits: Limits,
        error_bound: ErrorBound,
        rule: Rule,
        function: I,
        max_iterations: usize,
    ) -> Result<Self, Error> {
        match error_bound {
            ErrorBound::Absolute(v) => {
                if v <= 0.0 {
                    let kind = Kind::RelativeBoundNegativeOrZero;
                    return Err(Error::unevaluated(kind));
                }
            }
            ErrorBound::Relative(v) => {
                if v < 50.0 * f64::EPSILON {
                    let kind = Kind::AbsoluteBoundTooSmall;
                    return Err(Error::unevaluated(kind));
                }
            }
            ErrorBound::Either { absolute, relative } => {
                if absolute <= 0.0 && relative < 50.0 * f64::EPSILON {
                    let kind = Kind::InvalidTolerance;
                    return Err(Error::unevaluated(kind));
                }
            }
        }
        Ok(Self {
            limits,
            error_bound,
            rule,
            function,
            max_iterations,
        })
    }

    fn roundoff(result_abs: f64) -> f64 {
        50.0 * f64::EPSILON * result_abs
    }

    pub(crate) fn check_initial_integration(
        &self,
        initial: &Region,
    ) -> Result<Option<IntegralEstimate>, Error> {
        let tolerance = self.error_bound.tolerance(initial.result());
        let roundoff = Self::roundoff(initial.result_abs());

        if initial.error() <= roundoff && initial.error() > tolerance {
            let output = IntegralEstimate::from_region(initial, 1, self.rule.evaluations());
            let kind = Kind::RoundoffErrorDetected;

            Err(Error::new(kind, output))
        } else if (initial.error() <= tolerance
            && initial.error().to_bits() != initial.result_asc().to_bits())
            || initial.error() == 0.0
        {
            let output = IntegralEstimate::from_region(initial, 1, self.rule.evaluations());

            Ok(Some(output))
        } else if self.max_iterations == 1 {
            let output = IntegralEstimate::from_region(initial, 1, self.rule.evaluations());
            let kind = Kind::MaximumIterationsReached;

            Err(Error::new(kind, output))
        } else {
            Ok(None)
        }
    }

    /// Integrate the function and return a [`IntegralEstimate`] integration result.
    ///
    /// # Errors
    /// Integration can fail if user suplied tolerance cannot be achieved within the maximum number
    /// of iterations.
    pub fn integrate(&self) -> Result<IntegralEstimate, Error> {
        let initial = Basic::new(self.limits, self.rule, &self.function).integrate_internal();

        if let Some(output) = self.check_initial_integration(&initial)? {
            return Ok(output);
        }

        let mut workspace = self.initialise_workspace(initial);

        while workspace.iteration < self.max_iterations {
            let previous = workspace.retrieve_largest_error()?;

            let [lower, upper] = previous.bisect(&self.function, self.rule);

            let (result, error) = workspace.improved_result_error(&previous, &lower, &upper);

            let iteration_tolerance = self.error_bound.tolerance(result);

            workspace.push(lower);
            workspace.push(upper);

            if error <= iteration_tolerance {
                break;
            }

            workspace.check_roundoff()?;
            workspace.check_singularity()?;
        }

        let final_iteration = workspace.iteration;
        let output = workspace.integral_estimate();

        if final_iteration == self.max_iterations {
            let kind = Kind::MaximumIterationsReached;
            Err(Error::new(kind, output))
        } else {
            Ok(output)
        }
    }

    /// Return the integration [`Limits`]
    pub fn limits(&self) -> Limits {
        self.limits
    }

    fn initialise_workspace(&self, initial: Region) -> Workspace {
        let mut heap = BinaryHeap::with_capacity(2 * self.max_iterations + 1);

        let iteration = 1;
        let result = initial.result();
        let error = initial.error();
        let limits = initial.limits();
        let roundoff_count = 0;
        let roundoff_on_high_iteration_count = 0;
        let function_evaluations_per_integration = self.rule.evaluations();

        heap.push(initial);

        Workspace {
            heap,
            iteration,
            result,
            error,
            limits,
            roundoff_count,
            roundoff_on_high_iteration_count,
            function_evaluations_per_integration,
        }
    }
}

struct Workspace {
    heap: BinaryHeap<Region>,
    iteration: usize,
    result: f64,
    error: f64,
    limits: Limits,
    roundoff_count: usize,
    roundoff_on_high_iteration_count: usize,
    function_evaluations_per_integration: usize,
}

impl Workspace {
    fn retrieve_largest_error(&mut self) -> Result<Region, Error> {
        self.iteration += 1;
        if let Some(previous) = self.pop() {
            self.limits = previous.limits();
            Ok(previous)
        } else {
            let kind = Kind::UninitialisedWorkspace;
            Err(Error::unevaluated(kind))
        }
    }

    fn pop(&mut self) -> Option<Region> {
        self.heap.pop()
    }

    fn push(&mut self, integral: Region) {
        self.heap.push(integral);
    }

    fn improved_result_error(
        &mut self,
        previous: &Region,
        lower: &Region,
        upper: &Region,
    ) -> (f64, f64) {
        let prev_result = previous.result();
        let prev_error = previous.error();
        let new_result = lower.result() + upper.result();
        let new_error = lower.error() + upper.error();

        if lower.result_asc().to_bits() != lower.error().to_bits()
            && upper.result_asc().to_bits() != upper.error().to_bits()
        {
            let delta = (prev_result - new_result).abs();

            if delta <= 1e-5 * new_result.abs() && new_error >= 0.99 * prev_error {
                self.roundoff_count += 1;
            }
            if self.iteration >= 10 && new_error >= prev_error {
                self.roundoff_on_high_iteration_count += 1;
            }
        }

        self.result += new_result - prev_result;
        self.error += new_error - prev_error;
        let result = self.result;
        let error = self.error;

        (result, error)
    }

    fn check_roundoff(&self) -> Result<(), Error> {
        if self.roundoff_count >= 6 || self.roundoff_on_high_iteration_count >= 20 {
            let output = self.integral_estimate();
            let kind = Kind::RoundoffErrorDetected;
            return Err(Error::new(kind, output));
        }
        Ok(())
    }

    fn check_singularity(&self) -> Result<(), Error> {
        let limits = self.limits;
        if subinterval_too_small(limits) {
            let output = self.integral_estimate();
            let kind = Kind::BadIntegrandBehaviour { limits };
            Err(Error::new(kind, output))
        } else {
            Ok(())
        }
    }

    fn sum_results(&self) -> f64 {
        self.heap.iter().fold(0.0f64, |a, v| a + v.result())
    }

    fn integral_estimate(&self) -> IntegralEstimate {
        let result = self.sum_results();
        let error = self.error;
        let iterations = self.iteration;
        let function_evaluations = (2 * iterations - 1) * self.function_evaluations_per_integration;
        IntegralEstimate::new()
            .with_result(result)
            .with_error(error)
            .with_iterations(iterations)
            .with_function_evaluations(function_evaluations)
    }
}
