use num_traits::Zero;
use std::collections::binary_heap::BinaryHeap;

use crate::ScalarF64;

use crate::quadrature::{Integrator, Region, Rule};
use crate::{
    InitialisationError, IntegralEstimate, Integrand, IntegrationError, IntegrationErrorKind,
    Limits, Tolerance,
};

/// An adaptive Gauss-Kronrod integrator.
///
/// The [`Adaptive`] integrator applies a Gauss-Kronrod integration [`Rule`] to approximate the
/// integral of a one-dimensional function. Unlike the [`Basic`] routine, the routine implemented
/// by [`Adaptive`] is adaptive. After the initial integration, each iteration of the routine picks
/// the previous integration area which has the largest error estimate and bisects this region,
/// updating the estimate to the integal and the total approximated error. The adaptive routine
/// will return the first approximation, `result`, to the integral which has an absolute `error`
/// smaller than the tolerance `tol` encoded through the [`Tolerance`] enum, where
///
/// * [`Tolerance::Absolute`] specifies absolute tolerance and returns final estimate when
/// `error <= tol`,
/// * [`Tolerance::Relative`] specifies a relative error and returns final estimate when
/// `error <= tol * abs(result)`,  
/// * [`Tolerance::Either`] to return a result as soon as _either_ the relative or
/// absolute error bound has been satisfied.
///
/// The routine will end when _either_ one of the tolerance conditions have been satisfied _or_ an
/// error has occurred, see [`Error`] for more details.
///
/// This routine is primarily aimed towards the evaluation of integrals of relatively smooth
/// functions on finite integration domains that are free from singularities. If the function the
/// user wishes to integrate has integrable singularities and/or is defined on an infinite or
/// semi-infinite integration domain then the [`AdaptiveSingularity`] integrator should be
/// preferred. In fact, it can often be the case that the [`AdaptiveSingularity`] integrator is
/// _generally_ more efficient than the [`Adaptive`] routine though this should be benchmarked on a
/// case-by-case basis.
///
/// [`Basic`]: crate::quadrature::Basic
/// [`AdaptiveSingularity`]: crate::quadrature::AdaptiveSingularity
/// [`Error`]: crate::Error
/// [`Tolerance`]: crate::Tolerance
///
/// # Example
///
/// Here we present a calculation of the integral,
/// $$
/// I = \int_{0}^{1} x^{\alpha} \ln \frac{1}{x} dx = \frac{1}{(1+\alpha)^{2}}
/// $$
/// for different values of $\alpha$.
///
///```rust
/// use rint::{Integrand, Limits, Tolerance};
/// use rint::quadrature::{Adaptive, Basic, Rule};
///
/// use std::f64::consts::*;
///
/// struct Function1 {
///     alpha: f64,
/// }
///
/// impl Integrand for Function1 {
///     type Point = f64;
///     type Scalar = f64;
///
///     fn evaluate(&self, x: Self::Point) -> Self::Scalar {
///         let alpha = self.alpha;
///         x.powf(alpha) * (1.0 / x).ln()
///     }
/// }
///
/// # use std::error::Error;
/// # fn main() -> Result<(), Box<dyn Error>> {
/// const TOL: f64 = 1.0e-12;
///
/// let tolerance = Tolerance::Relative(TOL);
/// let limits = Limits::new(0.0, 1.0);
/// let rule = Rule::gk31();
/// let max_iterations = 1000;
///
/// let alpha_values = [2.6, PI, EULER_GAMMA, LOG10_E, 100.0, PI.powi(3)];
///
/// for alpha in alpha_values {
///     let function = Function1 { alpha };
///
///     let integral = Adaptive::new(&function, &rule, limits, tolerance, max_iterations)?
///                     .integrate()?;
///
///     let target = 1.0 / (1.0 + alpha).powi(2);
///     let result = integral.result();
///     let error = integral.error();
///     let abs_actual_error = (result - target).abs();
///     let tol = TOL * result.abs();
///
///     assert!(abs_actual_error < error);
///     assert!(error < tol);
/// }
/// # Ok(())
/// # }
///```
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Adaptive<'a, I> {
    function: &'a I,
    rule: &'a Rule,
    limits: Limits,
    tolerance: Tolerance,
    max_iterations: usize,
}

impl<'a, I> Adaptive<'a, I>
where
    I: Integrand<Point = f64>,
{
    /// Generate a new [`Adaptive`] integrator.
    ///
    /// Initialise an adaptive Gauss-Kronrod integrator. Arguments:
    /// - `function`: A user supplied function to be integrated which is something implementing the
    /// [`Integrand`] trait.
    /// - `rule`: An n-point Gauss-Kronrod integration [`Rule`]
    /// - `limits`: The interval over which the `function` should be integrated, [`Limits`].
    /// - `tolerance`: The tolerance requested by the user. Can be either an absolute tolerance
    /// or relative tolerance. Determines the exit condition of the integration routine, see
    /// [`Tolerance`].
    /// - `max_iterations`: The maximum number of iterations that the adaptive routine should use
    /// to try to satisfy the requested tolerance.
    ///
    /// # Errors
    /// Function can return an error if it receives bad user input. This is primarily related to
    /// using invalid values for the `tolerance`. The returned error is an [`InitialisationError`].
    /// See [`Tolerance`] and [`InitialisationError`] for more details.
    pub fn new(
        function: &'a I,
        rule: &'a Rule,
        limits: Limits,
        tolerance: Tolerance,
        max_iterations: usize,
    ) -> Result<Self, InitialisationError> {
        tolerance.check()?;

        Ok(Self {
            function,
            rule,
            limits,
            tolerance,
            max_iterations,
        })
    }

    /// Integrate the function and return a [`IntegralEstimate`] integration result.
    ///
    /// Adaptively applies the n-point Gauss-Kronrod integration [`Rule`] to the user supplied
    /// function implementing the [`Integrand`] trait to generate an [`IntegralEstimate`] upon
    /// successful completion.
    ///
    /// # Errors
    /// The error type [`IntegrationError`] will return both the error [`IntegrationErrorKind`] and
    /// the [`IntegralEstimate`] obtained before an error was encountered. The integration routine
    /// has several ways of failing:
    /// - The user supplied [`Tolerance`] could not be satisfied within the maximum number of
    /// iterations, see [`IntegrationErrorKind::MaximumIterationsReached`].
    ///
    /// - A roundoff error was detected. Can occur when the calculated numerical error from an
    /// internal integration is smaller than the estimated roundoff, but larger than the tolerance
    /// requested by the user, or when too many successive iterations do not reasonably improve the
    /// integral value and error estimate, see [`IntegrationErrorKind::RoundoffErrorDetected`].
    ///
    /// - Bisection of the highest error region into two subregions results in subregions with
    /// integraion limits that are too small, see [`IntegrationErrorKind::BadIntegrandBehaviour`].
    ///
    /// - An error is encountered when initialising the integration workspace. This is an internal
    /// error, which should not occur in user code, see
    /// [`IntegrationErrorKind::UninitialisedWorkspace`].
    pub fn integrate(&self) -> Result<IntegralEstimate<I::Scalar>, IntegrationError<I::Scalar>> {
        let initial = Integrator::new(self.function, self.rule, self.limits).integrate();

        if let Some(output) = self.check_initial_integration(&initial)? {
            return Ok(output);
        }

        let mut workspace = self.initialise_workspace(initial);

        while workspace.iteration < self.max_iterations {
            // XXX Increases workspace.iteration,
            // should this be taken out into own function for clarity?
            let previous = workspace.retrieve_largest_error()?;

            let [lower, upper] = previous.bisect(self.function, self.rule);

            let (result, error) = workspace.improved_result_error(&previous, &lower, &upper);

            let iteration_tolerance = self.tolerance.tolerance(&result);

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
            let kind = IntegrationErrorKind::MaximumIterationsReached(self.max_iterations);
            Err(IntegrationError::new(output, kind))
        } else {
            Ok(output)
        }
    }

    /// Return the integration [`Limits`]
    #[must_use]
    pub const fn limits(&self) -> Limits {
        self.limits
    }
}

impl<I: Integrand> Adaptive<'_, I> {
    fn initialise_workspace(&self, initial: Region<I::Scalar>) -> Workspace<I::Scalar> {
        let mut heap = BinaryHeap::with_capacity(2 * self.max_iterations + 1);

        let iteration = 1;
        let result = initial.result();
        let error = initial.error();
        let limits = initial.limits();
        let roundoff_count = 0;
        let roundoff_on_high_iteration_count = 0;
        let evaluations_per_integration = self.rule.evaluations();

        heap.push(initial);

        Workspace {
            heap,
            iteration,
            result,
            error,
            limits,
            roundoff_count,
            roundoff_on_high_iteration_count,
            evaluations_per_integration,
        }
    }

    /// Roundoff value used in diagnostics, matches GSL implementation.
    const fn roundoff(result_abs: f64) -> f64 {
        50.0 * f64::EPSILON * result_abs
    }

    /// Check if the initial integration was problematic such that the integration should be
    /// stopped immediately, whether the integration produced an error estimate which already
    /// satisfies the user supplied tolerance constraints, or if the integration was fine but
    /// should be further calculated using the adaptive algorithm to reduce the error estimate.
    pub(crate) fn check_initial_integration(
        &self,
        initial: &Region<I::Scalar>,
    ) -> Result<Option<IntegralEstimate<I::Scalar>>, IntegrationError<I::Scalar>> {
        let tolerance = self.tolerance.tolerance(&initial.result());
        let roundoff = Self::roundoff(initial.result_abs());

        if initial.error() <= roundoff && initial.error() > tolerance {
            let output = initial.estimate(1, self.rule.evaluations());
            let kind = IntegrationErrorKind::RoundoffErrorDetected;

            Err(IntegrationError::new(output, kind))
        } else if (initial.error() <= tolerance
            && initial.error().to_bits() != initial.result_asc().to_bits())
            || initial.error() == 0.0
        {
            let output = initial.estimate(1, self.rule.evaluations());

            Ok(Some(output))
        } else if self.max_iterations == 1 {
            let output = initial.estimate(1, self.rule.evaluations());
            let kind = IntegrationErrorKind::MaximumIterationsReached(self.max_iterations);

            Err(IntegrationError::new(output, kind))
        } else {
            Ok(None)
        }
    }
}

/// The integration workspace.
///
/// We use a workspace for maintaining the state of the adaptive integration routines. The `heap`
/// is a [`BinaryHeap`] storing [`Region`]s which have been integrated, and upon using
/// [`Workspace::pop`] returns the region with the highest `error` approximation which is the next
/// region to be bisected and re-integrated. Also tracked are the current estimates of `result` and
/// `error` across the entire integration region, the number of `iteration`s, and counts of times
/// that roundoff errors may have occurred, used as a diagnostic for terminating the integration
/// early if problematic roundoff behaviour is observed.
struct Workspace<T> {
    heap: BinaryHeap<Region<T>>,
    iteration: usize,
    result: T,
    error: f64,
    limits: Limits,
    roundoff_count: usize,
    roundoff_on_high_iteration_count: usize,
    evaluations_per_integration: usize,
}

impl<T: ScalarF64> Workspace<T> {
    /// Get the next [`Region`] in the [`BinaryHeap`] queue with the largest `error` estimate.
    fn retrieve_largest_error(&mut self) -> Result<Region<T>, IntegrationError<T>> {
        self.iteration += 1;
        if let Some(previous) = self.pop() {
            self.limits = previous.limits();
            Ok(previous)
        } else {
            let kind = IntegrationErrorKind::UninitialisedWorkspace;
            let output = IntegralEstimate::new();
            Err(IntegrationError::new(output, kind))
        }
    }

    fn pop(&mut self) -> Option<Region<T>> {
        self.heap.pop()
    }

    fn push(&mut self, integral: Region<T>) {
        self.heap.push(integral);
    }

    /// Determine the new estimate of the integral of a region which has been bisected.
    ///
    /// The region `previous` has been bisected into two new regions, `lower` and `upper`. This
    /// function returns the updated value of `previous.result` and `previous.error` using the
    /// newly calculated values from the bisected regions, including tests for roundoff errors to
    /// keep track of problematic or difficult regions in the integration.
    fn improved_result_error(
        &mut self,
        previous: &Region<T>,
        lower: &Region<T>,
        upper: &Region<T>,
    ) -> (T, f64) {
        let prev_result = previous.result();
        let prev_error = previous.error();
        let new_result = lower.result() + upper.result();
        let new_error = lower.error() + upper.error();

        // The comparisons was done here at the bit level as we only want to ignore this if we have
        // exact correspondence between the result_asc and error values. XXX should this just
        // compare up to f64::EPSILON, e.g. if (lower.result_asc() - lower.error()) > f64::EPSILON?
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

    /// Check the roundoff counters, triggering a roundoff error if the counts are too high. The
    /// conditions used here are the same as those found in the GSL numerical integration routines
    /// qag.c.
    fn check_roundoff(&self) -> Result<(), IntegrationError<T>> {
        if self.roundoff_count >= 6 || self.roundoff_on_high_iteration_count >= 20 {
            let output = self.integral_estimate();
            let kind = IntegrationErrorKind::RoundoffErrorDetected;
            return Err(IntegrationError::new(output, kind));
        }
        Ok(())
    }

    /// Check the integration for a non-integrable singularity. This is triggered if the next set
    /// of `limits` to be bisected (so has the highest error) corresponds to a region which is too
    /// small to be bisected again.
    fn check_singularity(&self) -> Result<(), IntegrationError<T>> {
        let limits = self.limits;
        if limits.subinterval_too_small() {
            let output = self.integral_estimate();
            let kind = IntegrationErrorKind::BadIntegrandBehaviour(limits);
            Err(IntegrationError::new(output, kind))
        } else {
            Ok(())
        }
    }

    fn sum_results(&self) -> T {
        self.heap
            .iter()
            .fold(<T as Zero>::zero(), |a, v| a + v.result())
    }

    fn integral_estimate(&self) -> IntegralEstimate<T> {
        let result = self.sum_results();
        let error = self.error;
        let iterations = self.iteration;
        let evaluations = (2 * iterations - 1) * self.evaluations_per_integration;
        IntegralEstimate::new()
            .with_result(result)
            .with_error(error)
            .with_iterations(iterations)
            .with_evaluations(evaluations)
    }
}
