use num_traits::Zero;

use crate::quadrature::rule::Rule;
use crate::quadrature::{Adaptive, AdaptiveSingularity, Basic, Error, Tolerance};
use crate::ScalarF64;
use crate::{Integrand, Limits};

/// The value of a function evaluated with Gauss-Kronrod integration and associated error
/// estimation.
#[derive(Debug)]
pub struct IntegralEstimate<T: ScalarF64> {
    result: T,
    error: f64,
    iterations: usize,
    function_evaluations: usize,
}

/// # Getters
impl<T: ScalarF64> IntegralEstimate<T> {
    /// Return the numerically approximated value of the integral.
    #[must_use]
    pub fn result(&self) -> T {
        self.result
    }

    /// Return the numerically approximated error.
    #[must_use]
    pub fn error(&self) -> f64 {
        self.error
    }

    /// Return the number of iterations used in integration.
    #[must_use]
    pub fn iterations(&self) -> usize {
        self.iterations
    }

    /// Return the number of function evaluations used in the integration.
    #[must_use]
    pub fn function_evaluations(&self) -> usize {
        self.function_evaluations
    }
}

/// # Integrators
impl<T: ScalarF64> IntegralEstimate<T> {
    /// Integrate a function using a basic (non-adaptive) Gauss-Kronrod integration rule, see
    /// [`Basic`] for details.
    ///
    /// Note: for the [`Basic`] integrator the number of iterations is `1`, and the number of
    /// function evaluations for an `n`-point integration rule is `n`.
    pub fn basic<I: Integrand>(
        limits: Limits,
        rule: Rule,
        function: I,
    ) -> IntegralEstimate<I::Scalar> {
        Basic::new(function, rule, limits).integrate()
    }

    /// Integrate a function using an adaptive Gauss-Kronrod integration routine, see [`Adaptive`]
    /// for details.
    ///
    /// # Errors
    /// The error type [`Error`] will return both the error [`Kind`] and the [`IntegralEstimate`]
    /// obtained before an error was encountered.
    /// The integration routine has several ways of failing:
    /// - Function will return an error if the user provided `Tolerance` does not satisfy the
    /// following constraints:
    ///     - `Tolerance::Absolute(v)` where `v > 0.0` (kind = [`Kind::AbsoluteBoundNegativeOrZero`])
    ///     - `Tolerance::Relative(v)` where `v > 50.0 * f64::EPSILON` (kind = [`Kind::RelativeBoundTooSmall`])
    ///     - `Tolerance::Either { absolute, relative }` where `absolute > 0.0 and relative > 50.0 * f64::EPSILON` (kind = [`Kind::InvalidTolerance`])
    /// - The user supplied [`Tolerance`] could not be satisfied within the maximum number of
    /// iterations (kind = [`Kind::MaximumIterationsReached`]).
    /// - A roundoff error was detected. Can occur when the calculated numerical error from an
    /// internal integration is smaller than the estimated roundoff, but larger than the tolerance
    /// requested by the user, or when too many successive iterations do not reasonably improve the
    /// integral value and error estimate (kind = [`Kind::RoundoffErrorDetected`]).
    /// - Bisection of the highest error region into two subregions results in subregions with
    /// integraion limits that are too small (kind = [`Kind::BadIntegrandBehaviour`]).
    /// - An error is encountered when initialising the integration workspace. This is an internal
    /// error, which should not occur downstream (kind = [`Kind::UninitialisedWorkspace`]).
    ///
    /// [`Kind`]: crate::quadrature::Kind
    /// [`Kind::DivergentOrSlowlyConverging`]: crate::quadrature::Kind#variant.DivergentOrSlowlyConverging
    /// [`Kind::BadIntegrandBehaviour`]: crate::quadrature::Kind#variant.BadIntegrandBehaviour
    /// [`Kind::DoesNotConverge`]: crate::quadrature::Kind#variant.DoesNotConverge
    /// [`Kind::UninitialisedWorkspace`]: crate::quadrature::Kind#variant.UninitialisedWorkspace
    /// [`Kind::MaximumIterationsReached`]: crate::quadrature::Kind#variant.MaximumIterationsReached
    /// [`Kind::RoundoffErrorDetected`]: crate::quadrature::Kind#variant.RoundoffErrorDetected
    /// [`Kind::InvalidTolerance`]: crate::quadrature::Kind#variant.InvalidTolerance
    /// [`Kind::RelativeBoundTooSmall`]: crate::quadrature::Kind#variant.RelativeBoundTooSmall
    /// [`Kind::AbsoluteBoundNegativeOrZero`]: crate::quadrature::Kind#variant.AbsoluteBoundNegativeOrZero
    pub fn adaptive<I: Integrand>(
        limits: Limits,
        error_bound: Tolerance,
        rule: Rule,
        function: I,
        max_iterations: usize,
    ) -> Result<IntegralEstimate<I::Scalar>, Error<I::Scalar>> {
        Adaptive::new(function, rule, limits, error_bound, max_iterations)?.integrate()
    }

    /// Integrate the function with possible integrable singularities and a finite integration
    /// interval, see [`AdaptiveSingularity`] for details.
    ///
    /// # Errors
    /// The error type [`Error`] will return both the error [`Kind`] and the [`IntegralEstimate`]
    /// obtained before an error was encountered.
    /// The integration routine has several ways of failing:
    /// - Function will return an error if the user provided `Tolerance` does not satisfy the
    /// following constraints:
    ///     - `Tolerance::Absolute(v)` where `v > 0.0` (kind = [`Kind::AbsoluteBoundNegativeOrZero`])
    ///     - `Tolerance::Relative(v)` where `v > 50.0 * f64::EPSILON` (kind = [`Kind::RelativeBoundTooSmall`])
    ///     - `Tolerance::Either { absolute, relative }` where `absolute > 0.0 and relative > 50.0 * f64::EPSILON` (kind = [`Kind::InvalidTolerance`])
    /// - The user supplied [`Tolerance`] could not be satisfied within the maximum number of
    /// iterations (kind = [`Kind::MaximumIterationsReached`]).
    /// - A roundoff error was detected. Can occur when the calculated numerical error from an
    /// internal integration is smaller than the estimated roundoff, but larger than the tolerance
    /// requested by the user, or when too many successive iterations do not reasonably improve the
    /// integral value and error estimate (kind = [`Kind::RoundoffErrorDetected`]).
    /// - Bisection of the highest error region into two subregions results in subregions with
    /// integraion limits that are too small (kind = [`Kind::BadIntegrandBehaviour`]).
    /// - The integral does not converge. Occurs if at least 5 bisections with extrapolation have
    /// failed to improve the integral value and error estimates (kind = [`Kind::DoesNotConverge`])
    /// - The integral is divergent or slowly convergent. Occurs if the ratio of the extrapolation
    /// improved integral estimate and the non-extrapolated estimate is too large or small, or if
    /// the calculated error is larger than the calculated value of the integral (kind =
    /// [`Kind::DivergentOrSlowlyConverging`])
    /// - An error is encountered when initialising the integration workspace. This is an internal
    /// error, which should not occur downstream (kind = [`Kind::UninitialisedWorkspace`]).
    ///
    /// [`Kind`]: crate::quadrature::Kind
    /// [`Kind::DivergentOrSlowlyConverging`]: crate::quadrature::Kind#variant.DivergentOrSlowlyConverging
    /// [`Kind::BadIntegrandBehaviour`]: crate::quadrature::Kind#variant.BadIntegrandBehaviour
    /// [`Kind::DoesNotConverge`]: crate::quadrature::Kind#variant.DoesNotConverge
    /// [`Kind::UninitialisedWorkspace`]: crate::quadrature::Kind#variant.UninitialisedWorkspace
    /// [`Kind::MaximumIterationsReached`]: crate::quadrature::Kind#variant.MaximumIterationsReached
    /// [`Kind::RoundoffErrorDetected`]: crate::quadrature::Kind#variant.RoundoffErrorDetected
    /// [`Kind::InvalidTolerance`]: crate::quadrature::Kind#variant.InvalidTolerance
    /// [`Kind::RelativeBoundTooSmall`]: crate::quadrature::Kind#variant.RelativeBoundTooSmall
    /// [`Kind::AbsoluteBoundNegativeOrZero`]: crate::quadrature::Kind#variant.AbsoluteBoundNegativeOrZero
    pub fn adaptive_singularity_finite<I: Integrand>(
        limits: Limits,
        error_bound: Tolerance,
        function: I,
        max_iterations: usize,
    ) -> Result<IntegralEstimate<I::Scalar>, Error<I::Scalar>> {
        AdaptiveSingularity::finite(function, limits, error_bound, max_iterations)?.integrate()
    }

    /// Integrate the function with possible integrable singularities and an infinite integration
    /// interval, see [`AdaptiveSingularity`] for details.
    ///
    /// # Errors
    /// The error type [`Error`] will return both the error [`Kind`] and the [`IntegralEstimate`]
    /// obtained before an error was encountered.
    /// The integration routine has several ways of failing:
    /// - Function will return an error if the user provided `Tolerance` does not satisfy the
    /// following constraints:
    ///     - `Tolerance::Absolute(v)` where `v > 0.0` (kind = [`Kind::AbsoluteBoundNegativeOrZero`])
    ///     - `Tolerance::Relative(v)` where `v > 50.0 * f64::EPSILON` (kind = [`Kind::RelativeBoundTooSmall`])
    ///     - `Tolerance::Either { absolute, relative }` where `absolute > 0.0 and relative > 50.0 * f64::EPSILON` (kind = [`Kind::InvalidTolerance`])
    /// - The user supplied [`Tolerance`] could not be satisfied within the maximum number of
    /// iterations (kind = [`Kind::MaximumIterationsReached`]).
    /// - A roundoff error was detected. Can occur when the calculated numerical error from an
    /// internal integration is smaller than the estimated roundoff, but larger than the tolerance
    /// requested by the user, or when too many successive iterations do not reasonably improve the
    /// integral value and error estimate (kind = [`Kind::RoundoffErrorDetected`]).
    /// - Bisection of the highest error region into two subregions results in subregions with
    /// integraion limits that are too small (kind = [`Kind::BadIntegrandBehaviour`]).
    /// - The integral does not converge. Occurs if at least 5 bisections with extrapolation have
    /// failed to improve the integral value and error estimates (kind = [`Kind::DoesNotConverge`])
    /// - The integral is divergent or slowly convergent. Occurs if the ratio of the extrapolation
    /// improved integral estimate and the non-extrapolated estimate is too large or small, or if
    /// the calculated error is larger than the calculated value of the integral (kind =
    /// [`Kind::DivergentOrSlowlyConverging`])
    /// - An error is encountered when initialising the integration workspace. This is an internal
    /// error, which should not occur downstream (kind = [`Kind::UninitialisedWorkspace`]).
    ///
    /// [`Kind`]: crate::quadrature::Kind
    /// [`Kind::DivergentOrSlowlyConverging`]: crate::quadrature::Kind#variant.DivergentOrSlowlyConverging
    /// [`Kind::BadIntegrandBehaviour`]: crate::quadrature::Kind#variant.BadIntegrandBehaviour
    /// [`Kind::DoesNotConverge`]: crate::quadrature::Kind#variant.DoesNotConverge
    /// [`Kind::UninitialisedWorkspace`]: crate::quadrature::Kind#variant.UninitialisedWorkspace
    /// [`Kind::MaximumIterationsReached`]: crate::quadrature::Kind#variant.MaximumIterationsReached
    /// [`Kind::RoundoffErrorDetected`]: crate::quadrature::Kind#variant.RoundoffErrorDetected
    /// [`Kind::InvalidTolerance`]: crate::quadrature::Kind#variant.InvalidTolerance
    /// [`Kind::RelativeBoundTooSmall`]: crate::quadrature::Kind#variant.RelativeBoundTooSmall
    /// [`Kind::AbsoluteBoundNegativeOrZero`]: crate::quadrature::Kind#variant.AbsoluteBoundNegativeOrZero
    pub fn adaptive_singularity_infinite<I: Integrand>(
        error_bound: Tolerance,
        function: I,
        max_iterations: usize,
    ) -> Result<IntegralEstimate<I::Scalar>, Error<I::Scalar>> {
        AdaptiveSingularity::infinite(function, error_bound, max_iterations)?.integrate()
    }

    /// Integrate the function with possible integrable singularities and semi-infinite integration
    /// interval (b, Inf), see [`AdaptiveSingularity`] for details.
    ///
    /// # Errors
    /// The error type [`Error`] will return both the error [`Kind`] and the [`IntegralEstimate`]
    /// obtained before an error was encountered.
    /// The integration routine has several ways of failing:
    /// - Function will return an error if the user provided `Tolerance` does not satisfy the
    /// following constraints:
    ///     - `Tolerance::Absolute(v)` where `v > 0.0` (kind = [`Kind::AbsoluteBoundNegativeOrZero`])
    ///     - `Tolerance::Relative(v)` where `v > 50.0 * f64::EPSILON` (kind = [`Kind::RelativeBoundTooSmall`])
    ///     - `Tolerance::Either { absolute, relative }` where `absolute > 0.0 and relative > 50.0 * f64::EPSILON` (kind = [`Kind::InvalidTolerance`])
    /// - The user supplied [`Tolerance`] could not be satisfied within the maximum number of
    /// iterations (kind = [`Kind::MaximumIterationsReached`]).
    /// - A roundoff error was detected. Can occur when the calculated numerical error from an
    /// internal integration is smaller than the estimated roundoff, but larger than the tolerance
    /// requested by the user, or when too many successive iterations do not reasonably improve the
    /// integral value and error estimate (kind = [`Kind::RoundoffErrorDetected`]).
    /// - Bisection of the highest error region into two subregions results in subregions with
    /// integraion limits that are too small (kind = [`Kind::BadIntegrandBehaviour`]).
    /// - The integral does not converge. Occurs if at least 5 bisections with extrapolation have
    /// failed to improve the integral value and error estimates (kind = [`Kind::DoesNotConverge`])
    /// - The integral is divergent or slowly convergent. Occurs if the ratio of the extrapolation
    /// improved integral estimate and the non-extrapolated estimate is too large or small, or if
    /// the calculated error is larger than the calculated value of the integral (kind =
    /// [`Kind::DivergentOrSlowlyConverging`])
    /// - An error is encountered when initialising the integration workspace. This is an internal
    /// error, which should not occur downstream (kind = [`Kind::UninitialisedWorkspace`]).
    ///
    /// [`Kind`]: crate::quadrature::Kind
    /// [`Kind::DivergentOrSlowlyConverging`]: crate::quadrature::Kind#variant.DivergentOrSlowlyConverging
    /// [`Kind::BadIntegrandBehaviour`]: crate::quadrature::Kind#variant.BadIntegrandBehaviour
    /// [`Kind::DoesNotConverge`]: crate::quadrature::Kind#variant.DoesNotConverge
    /// [`Kind::UninitialisedWorkspace`]: crate::quadrature::Kind#variant.UninitialisedWorkspace
    /// [`Kind::MaximumIterationsReached`]: crate::quadrature::Kind#variant.MaximumIterationsReached
    /// [`Kind::RoundoffErrorDetected`]: crate::quadrature::Kind#variant.RoundoffErrorDetected
    /// [`Kind::InvalidTolerance`]: crate::quadrature::Kind#variant.InvalidTolerance
    /// [`Kind::RelativeBoundTooSmall`]: crate::quadrature::Kind#variant.RelativeBoundTooSmall
    /// [`Kind::AbsoluteBoundNegativeOrZero`]: crate::quadrature::Kind#variant.AbsoluteBoundNegativeOrZero
    pub fn adaptive_singularity_semi_infinite_upper<I: Integrand>(
        lower: f64,
        error_bound: Tolerance,
        function: I,
        max_iterations: usize,
    ) -> Result<IntegralEstimate<I::Scalar>, Error<I::Scalar>> {
        AdaptiveSingularity::semi_infinite_upper(function, lower, error_bound, max_iterations)?
            .integrate()
    }

    /// Integrate the function with possible integrable singularities and semi-infinite integration
    /// interval (-Inf, a), see [`AdaptiveSingularity`] for details.
    ///
    /// # Errors
    /// The error type [`Error`] will return both the error [`Kind`] and the [`IntegralEstimate`]
    /// obtained before an error was encountered.
    /// The integration routine has several ways of failing:
    /// - Function will return an error if the user provided `Tolerance` does not satisfy the
    /// following constraints:
    ///     - `Tolerance::Absolute(v)` where `v > 0.0` (kind = [`Kind::AbsoluteBoundNegativeOrZero`])
    ///     - `Tolerance::Relative(v)` where `v > 50.0 * f64::EPSILON` (kind = [`Kind::RelativeBoundTooSmall`])
    ///     - `Tolerance::Either { absolute, relative }` where `absolute > 0.0 and relative > 50.0 * f64::EPSILON` (kind = [`Kind::InvalidTolerance`])
    /// - The user supplied [`Tolerance`] could not be satisfied within the maximum number of
    /// iterations (kind = [`Kind::MaximumIterationsReached`]).
    /// - A roundoff error was detected. Can occur when the calculated numerical error from an
    /// internal integration is smaller than the estimated roundoff, but larger than the tolerance
    /// requested by the user, or when too many successive iterations do not reasonably improve the
    /// integral value and error estimate (kind = [`Kind::RoundoffErrorDetected`]).
    /// - Bisection of the highest error region into two subregions results in subregions with
    /// integraion limits that are too small (kind = [`Kind::BadIntegrandBehaviour`]).
    /// - The integral does not converge. Occurs if at least 5 bisections with extrapolation have
    /// failed to improve the integral value and error estimates (kind = [`Kind::DoesNotConverge`])
    /// - The integral is divergent or slowly convergent. Occurs if the ratio of the extrapolation
    /// improved integral estimate and the non-extrapolated estimate is too large or small, or if
    /// the calculated error is larger than the calculated value of the integral (kind =
    /// [`Kind::DivergentOrSlowlyConverging`])
    /// - An error is encountered when initialising the integration workspace. This is an internal
    /// error, which should not occur downstream (kind = [`Kind::UninitialisedWorkspace`]).
    ///
    /// [`Kind`]: crate::quadrature::Kind
    /// [`Kind::DivergentOrSlowlyConverging`]: crate::quadrature::Kind#variant.DivergentOrSlowlyConverging
    /// [`Kind::BadIntegrandBehaviour`]: crate::quadrature::Kind#variant.BadIntegrandBehaviour
    /// [`Kind::DoesNotConverge`]: crate::quadrature::Kind#variant.DoesNotConverge
    /// [`Kind::UninitialisedWorkspace`]: crate::quadrature::Kind#variant.UninitialisedWorkspace
    /// [`Kind::MaximumIterationsReached`]: crate::quadrature::Kind#variant.MaximumIterationsReached
    /// [`Kind::RoundoffErrorDetected`]: crate::quadrature::Kind#variant.RoundoffErrorDetected
    /// [`Kind::InvalidTolerance`]: crate::quadrature::Kind#variant.InvalidTolerance
    /// [`Kind::RelativeBoundTooSmall`]: crate::quadrature::Kind#variant.RelativeBoundTooSmall
    /// [`Kind::AbsoluteBoundNegativeOrZero`]: crate::quadrature::Kind#variant.AbsoluteBoundNegativeOrZero
    pub fn adaptive_singularity_semi_infinite_lower<I: Integrand>(
        upper: f64,
        error_bound: Tolerance,
        function: I,
        max_iterations: usize,
    ) -> Result<IntegralEstimate<I::Scalar>, Error<I::Scalar>> {
        AdaptiveSingularity::semi_infinite_lower(function, upper, error_bound, max_iterations)?
            .integrate()
    }
}

impl<T: ScalarF64> IntegralEstimate<T> {
    pub(crate) fn new() -> Self {
        let result = <T as Zero>::zero();
        Self {
            result,
            error: 0.0,
            iterations: 0,
            function_evaluations: 0,
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

    pub(crate) fn with_function_evaluations(mut self, function_evaluations: usize) -> Self {
        self.function_evaluations = function_evaluations;
        self
    }
}
