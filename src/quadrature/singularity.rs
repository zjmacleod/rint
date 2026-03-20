use std::collections::binary_heap::BinaryHeap;

use crate::quadrature::{Integrator, Region, Rule};
use crate::{
    InitialisationError, IntegralEstimate, Integrand, IntegrationError, IntegrationErrorKind,
    Limits, ScalarF64, Tolerance,
};

/// An adaptive Gauss-Kronrod integrator with extrapolation for integrable singularities.
///
/// # Overview
///
/// The [`AdaptiveSingularity`] integrator applies a Gauss-Kronrod integration [`Rule`] to
/// approximate the integral of a one-dimensional function. Unlike the [`Basic`] routine, the
/// routine implemented by [`AdaptiveSingularity`] is adaptive. After the initial integration, each
/// iteration of the routine picks the previous integration area which has the largest error
/// estimate and bisects this region, updating the estimate to the integal and the total
/// approximated error. Adaptive routines concentrate new subintervals around the region of highest
/// error. If this region contains an integrable singularity, then the adaptive routine of
/// [`Adaptive`] may fail to obtain a suitable estimate. To overcome this, [`AdaptiveSingularity`]
/// uses an additional Wynn epsilon-algorithm extrapolation proceedure to manage these integrable
/// singularities and provide a suitable estimate. The combination of adaptive iteration and
/// extrapolation makes this routine suitable as a general purpose integrator, and will often
/// out-perform the [`Adaptive`] routine (though this should be benchmarked on a case-by-case
/// basis). It is valid for integrals defined on both finite _and_ (semi-)infinite integration
/// intervals. A choice of Gauss-Kronrod integration rule is not required for the
/// [`AdaptiveSingularity`] integrator. Instead, for general functions on finite intervals the
/// 21-point rule [`Rule::gk21()`] is used, while for (semi-)infinite intervals (see below) the
/// lower 15-point rule [`Rule::gk15()`] is used.
///
/// The adaptive routine will return the first approximation, `result`, to the integral which has
/// an absolute `error` smaller than the tolerance `tol` encoded through the [`Tolerance`] enum,
/// where
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
/// # Finite, infinite and semi-infinite integration regions
///
/// As well as being able to calculate integrals over integrable singularities within finite
/// integration limits , [`AdaptiveSingularity`] can also be used for integrations over infinite
/// and semi-infinite integration limits. Each of these has a dedicated constructor:
/// - [`AdaptiveSingularity::finite`]: finite interval $x \in (b,a)$
/// - [`AdaptiveSingularity::infinite`]: infinite interval $x \in (-\infty,+\infty)$
/// - [`AdaptiveSingularity::semi_infinite_lower`]: infinite lower limit $x \in (-\infty,a)$
/// - [`AdaptiveSingularity::semi_infinite_upper`]: infinite upper limit, $x \in (b,+\infty)$
///
/// The integrations over (semi-)infinite intervals are managed by transforming to a finite
/// interval. This can introduce integrable singularities even to smooth functions, which is why
/// the additional extrapolation provided by [`AdaptiveSingularity`] should be used for these.
///
///
/// ## Example: finite integration region
///
/// Here we present a calculation of [Catalan's constant] $G$ using the integral representation:
/// $$
/// G = \int_{0}^{1} \frac{\mathrm{arcsinh} x}{\sqrt{1 - x^{2}}} d x
/// $$
/// which has a finite integration region $x \in (0,1)$ but difficult integral behaviour as the
/// integration variable $x \to 1$. We use [`AdaptiveSingularity::finite`] to approximate $G$.
///```rust
/// use rint::{Integrand, Limits, Tolerance};
/// use rint::quadrature::AdaptiveSingularity;
///
/// const G: f64 = 0.915_965_594_177_219_015_054_603_514_932_384_110_774;
/// const TOL: f64 = 1.0e-12;
///
/// struct Catalan;
///
/// impl Integrand for Catalan {
///     type Point = f64;
///     type Scalar = f64;
///
///     fn evaluate(&self, x: &Self::Point) -> Self::Scalar {
///         x.asinh() / (1.0 - x.powi(2)).sqrt()
///     }
/// }
///
/// # use std::error::Error;
/// # fn main() -> Result<(), Box<dyn Error>> {
/// let catalan = Catalan;
/// let limits = Limits::new(0.0,1.0);
/// let tolerance = Tolerance::Relative(TOL);
/// let integral = AdaptiveSingularity::finite(catalan, limits, tolerance, 1000)?
///     .integrate()?;  
///
/// let result = integral.result();
/// let error = integral.error();
/// let abs_actual_error = (G - result).abs();
/// let tol = TOL * result.abs();
/// let iters = integral.iterations();
/// assert_eq!(iters, 10);
/// assert!(abs_actual_error < error);
/// assert!(error < tol);
/// # Ok(())
/// # }
///```
///
///
/// ## Example: infinite integration region
///
/// Here we present a calculation of [Euler's constant] $\gamma$ ([`std::f64::consts::EULER_GAMMA`])
/// using the integral representation:
/// $$
/// \gamma = -\int_{-\infty}^{+\infty} x e^{x - e^{x}} dx
/// $$
/// which has an infinite integration region $x\in(-\infty,+\infty)$. We use [`AdaptiveSingularity::infinite`]
/// to approximate $\gamma$.
///```rust
/// use rint::{Integrand, Tolerance};
/// use rint::quadrature::AdaptiveSingularity;
///
/// const GAMMA: f64 = std::f64::consts::EULER_GAMMA;
/// const TOL: f64 = 1.0e-12;
///
/// struct Gamma;
///
/// impl Integrand for Gamma {
///     type Point = f64;
///     type Scalar = f64;
///
///     fn evaluate(&self, x: &Self::Point) -> Self::Scalar {
///         let exponent = x - x.exp();
///         -x * exponent.exp()
///     }
/// }
///
/// # use std::error::Error;
/// # fn main() -> Result<(), Box<dyn Error>> {
/// let gamma = Gamma;
/// let tolerance = Tolerance::Relative(TOL);
/// let integral = AdaptiveSingularity::infinite(gamma, tolerance, 1000)?
///     .integrate()?;
///
/// let result = integral.result();
/// let error = integral.error();
/// let abs_actual_error = (GAMMA - result).abs();
/// let tol = TOL * result.abs();
/// let iters = integral.iterations();
/// assert_eq!(iters, 11);
/// assert!(abs_actual_error < error);
/// assert!(error < tol);
/// # Ok(())
/// # }
///```
///
///
/// ## Example: semi-infinite integration region
///
/// Here we present a calculation of [Catalan's constant] $G$ using the integral representation:
/// $$
/// G = \frac{\pi}{2} \int_{1}^{+\infty} \frac{(x^{4} - 6 x^{2} + 1) \ln \ln x}{(1+x^{2})^{3}} d x
/// $$
/// which has a semi-infinite integration region $x \in (1,+\infty)$ We use
/// [`AdaptiveSingularity::semi_infinite_upper`] to approximate $G$.
///```rust
/// use rint::{Integrand, Limits, Tolerance};
/// use rint::quadrature::AdaptiveSingularity;
///
/// const G: f64 = 0.915_965_594_177_219_015_054_603_514_932_384_110_774;
/// const TOL: f64 = 1.0e-12;
///
/// struct Catalan;
///
/// impl Integrand for Catalan {
///     type Point = f64;
///     type Scalar = f64;
///
///     fn evaluate(&self, x: &Self::Point) -> Self::Scalar {
///         let FRAC_PI_2 = std::f64::consts::FRAC_PI_2;
///         let polynomial = x.powi(4) - 6.0 * x.powi(2) + 1.0;
///         let denominator = (1.0 + x.powi(2)).powi(3);
///         let lnlnx = x.ln().ln();
///         
///         FRAC_PI_2 * polynomial * lnlnx / denominator
///     }
/// }
///
/// # use std::error::Error;
/// # fn main() -> Result<(), Box<dyn Error>> {
/// let catalan = Catalan;
/// let lower = 1.0;
/// let tolerance = Tolerance::Relative(TOL);
/// let integral = AdaptiveSingularity::semi_infinite_upper(catalan, lower, tolerance, 1000)?
///     .integrate()?;
///
/// let result = integral.result();
/// let error = integral.error();
/// let abs_actual_error = (G - result).abs();
/// let tol = TOL * result.abs();
/// let iters = integral.iterations();
/// assert_eq!(iters, 32);
/// assert!(abs_actual_error < error);
/// assert!(error < tol);
/// # Ok(())
/// # }
///```
///
/// [Catalan's constant]: <https://en.wikipedia.org/wiki/Catalan%27s_constant>
/// [Euler's constant]: <https://en.wikipedia.org/wiki/Euler%27s_constant>
/// [`Basic`]: crate::quadrature::Basic
/// [`Adaptive`]: crate::quadrature::Adaptive
/// [`Tolerance`]: crate::Tolerance
/// [`Error`]: crate::Error
///
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct AdaptiveSingularity<I> {
    function: I,
    rule: Rule,
    limits: Limits,
    tolerance: Tolerance,
    max_iterations: usize,
    evaluations_multiplier: usize,
}

impl<I> AdaptiveSingularity<I>
where
    I: Integrand<Point = f64>,
{
    /// Generate a new [`AdaptiveSingularity`] integrator for functions with possible integrable
    /// singularities on finite intervals.
    ///
    ///
    /// Initialise an adaptive Gauss-Kronrod integrator with extrapolation for approximating
    /// integrals of the form,
    /// $$
    /// \int_{b}^{a} f(x) dx.
    /// $$
    /// with finite integration bounds $x \in (b,a)$.
    ///
    ///  Arguments:
    /// - `function`: A user supplied function to be integrated which is something implementing the
    /// [`Integrand`] trait.
    /// - `limits`: The interval over which the `function` should be integrated, [`Limits`].
    /// - `tolerance`: The tolerance requested by the user. Can be either an absolute tolerance
    /// or relative tolerance. Determines the exit condition of the integration routine, see
    /// [`Tolerance`].
    /// - `max_iterations`: The maximum number of iterations that the routine should use to try to
    /// satisfy the requested tolerance.
    ///
    /// # Errors
    /// Function can return an error if it receives bad user input. This is primarily related to
    /// using invalid values for the `tolerance`. The returned error is an [`InitialisationError`].
    /// See [`Tolerance`] and [`InitialisationError`] for more details.
    ///
    ///
    /// # Example
    ///
    /// Here we present a calculation of [Catalan's constant] $G$ using the integral representation:
    /// $$
    /// G = \int_{0}^{1} \frac{\mathrm{arcsinh} x}{\sqrt{1 - x^{2}}} d x
    /// $$
    /// which has a finite integration region $x \in (0,1)$ but difficult integral behaviour as the
    /// integration variable $x \to 1$. We use [`AdaptiveSingularity::finite`] to approximate $G$.
    ///```rust
    /// use rint::{Integrand, Limits, Tolerance};
    /// use rint::quadrature::AdaptiveSingularity;
    ///
    /// const G: f64 = 0.915_965_594_177_219_015_054_603_514_932_384_110_774;
    /// const TOL: f64 = 1.0e-12;
    ///
    /// struct Catalan;
    ///
    /// impl Integrand for Catalan {
    ///     type Point = f64;
    ///     type Scalar = f64;
    ///
    ///     fn evaluate(&self, x: &Self::Point) -> Self::Scalar {
    ///         x.asinh() / (1.0 - x.powi(2)).sqrt()
    ///     }
    /// }
    ///
    /// # use std::error::Error;
    /// # fn main() -> Result<(), Box<dyn Error>> {
    /// let catalan = Catalan;
    /// let limits = Limits::new(0.0,1.0);
    /// let tolerance = Tolerance::Relative(TOL);
    /// let integral = AdaptiveSingularity::finite(catalan, limits, tolerance, 1000)?
    ///     .integrate()?;
    ///
    /// let result = integral.result();
    /// let error = integral.error();
    /// let abs_actual_error = (G - result).abs();
    /// let tol = TOL * result.abs();
    /// let iters = integral.iterations();
    /// assert_eq!(iters, 10);
    /// assert!(abs_actual_error < error);
    /// assert!(error < tol);
    /// # Ok(())
    /// # }
    ///```
    ///
    pub fn finite(
        function: I,
        limits: Limits,
        tolerance: Tolerance,
        max_iterations: usize,
    ) -> Result<Self, InitialisationError> {
        let rule = Rule::gk21();
        Self::new(function, rule, limits, tolerance, max_iterations)
    }

    /// Integrate the function and return a [`IntegralEstimate`] integration result.
    ///
    /// Applies adaptive Gauss-Kronrod integration with extrapolation via the Wynn epsilon-algorithm
    /// to the user supplied function implementing the [`Integrand`] trait to generate an
    /// [`IntegralEstimate`] upon successful completion.
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
    /// - The integral does not converge. Occurs if at least 5 bisections with extrapolation have
    /// failed to improve the integral value and error estimates, see
    /// [`IntegrationErrorKind::DoesNotConverge`].
    ///
    /// - The integral is divergent or slowly convergent. Occurs if the ratio of the extrapolation
    /// improved integral estimate and the non-extrapolated estimate is too large or small, or if
    /// the calculated error is larger than the calculated value of the integral, see
    /// [`IntegrationErrorKind::DivergentOrSlowlyConverging`].
    ///
    /// - An error is encountered when initialising the integration workspace. This is an internal
    /// error, which should not occur downstream , see
    /// [`IntegrationErrorKind::UninitialisedWorkspace`].
    ///
    pub fn integrate(&self) -> Result<IntegralEstimate<I::Scalar>, IntegrationError<I::Scalar>> {
        let initial = Integrator::new(&self.function, &self.rule, self.limits).integrate();

        if let Some(output) = self.check_initial_integration(&initial)? {
            return Ok(output);
        }

        let mut workspace = self.initialise_workspace(initial);

        while workspace.iteration < self.max_iterations {
            let previous = workspace.retrieve_largest_error()?;

            let current_interval = previous.abs_interval_length();

            let [lower, upper] = previous.bisect(&self.function, &self.rule);

            let previous_error = previous.error();
            let iteration_error = lower.error() + upper.error();

            let (result, error) = workspace.improved_result_error(&previous, &lower, &upper);

            let iteration_tolerance = self.tolerance.tolerance(&result);

            workspace.update(lower, upper);

            if error <= iteration_tolerance {
                return workspace.compute_result();
            }

            if workspace.is_err() {
                break;
            }

            if workspace.iteration >= self.max_iterations - 1 {
                workspace.set_error_kind(IntegrationErrorKind::MaximumIterationsReached(
                    self.max_iterations,
                ));
                break;
            }

            if workspace.iteration == 2 {
                workspace.smallest_interval *= 0.375;
                workspace.large_interval_error = workspace.error;
                workspace.table.tolerance = iteration_tolerance;
                workspace.table.append_table(workspace.result);
                continue;
            }

            workspace.update_large_interval_error(
                current_interval,
                previous_error,
                iteration_error,
            );

            // We want to reduce the error in the [`Region`]s with large integration intervals
            // first. If we are _not_ extrapolating, compare the length of the next [`Region`]
            // with the highest `error` with the smallest interval previously encountered. Iff
            // the next [`Region`] is LARGER, return to the start of the loop. Otherwise, the
            // next iteration has the highest error and smallest interval. Store it in the
            // workspace.store and set extrapolate to true.
            if !workspace.extrapolate {
                // peek the next interval with the largest error. If the interval length is larger
                // than the smallest interval previously encountered then go to the start of the
                // loop.
                if let Some(next) = workspace.peek() {
                    if next.abs_interval_length() > workspace.smallest_interval {
                        continue;
                    }
                }

                // the next interval to be bisected has the largest error and
                // smallest interval. Store it to reduce error in large intervals
                if let Some(smallest) = workspace.pop() {
                    workspace.store(smallest);
                }
                workspace.extrapolate = true;
            }

            if !workspace.table.error_detected
                && workspace.large_interval_error > workspace.table.tolerance
            {
                if workspace.remaining_large_intervals() {
                    continue;
                }

                // If we reach here then _all_ the results are now in the store
                // so we mem::swap the store and the heap and continue on to the
                // extrapolation.
                std::mem::swap(&mut workspace.heap, &mut workspace.store);
            }

            let (ext_result, ext_error) = workspace.table.extrapolate(result);

            workspace.check_convergence();

            if ext_error < workspace.table.error {
                let ext_tolerance = self.tolerance.tolerance(&ext_result);
                workspace.update_table_values(ext_result, ext_error, ext_tolerance);
                if workspace.table.error <= workspace.table.tolerance {
                    break;
                }
            }

            if let Some(IntegrationErrorKind::DoesNotConverge) = workspace.error_kind {
                break;
            }

            workspace.prepare_next_iteration();
        }

        workspace.check_error_and_compute()
    }

    /// Return the integration [`Limits`].
    pub const fn limits(&self) -> Limits {
        self.limits
    }
}

impl<I> AdaptiveSingularity<I>
where
    I: Integrand<Point = f64>,
{
    fn new(
        function: I,
        rule: Rule,
        limits: Limits,
        tolerance: Tolerance,
        max_iterations: usize,
    ) -> Result<Self, InitialisationError> {
        let evaluations_multiplier = 1;
        tolerance.check()?;
        Ok(Self {
            function,
            rule,
            limits,
            tolerance,
            max_iterations,
            evaluations_multiplier,
        })
    }

    fn new_with_evaluations_multiplier(
        function: I,
        rule: Rule,
        limits: Limits,
        tolerance: Tolerance,
        max_iterations: usize,
        evaluations_multiplier: usize,
    ) -> Result<Self, InitialisationError> {
        let mut v = Self::new(function, rule, limits, tolerance, max_iterations)?;
        v.evaluations_multiplier = evaluations_multiplier;
        Ok(v)
    }

    const fn roundoff(result_abs: f64) -> f64 {
        100.0 * f64::EPSILON * result_abs
    }

    fn check_initial_integration(
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

    fn initialise_workspace(&self, initial: Region<I::Scalar>) -> Workspace<I::Scalar> {
        let mut heap = BinaryHeap::with_capacity(2 * self.max_iterations + 1);

        let iteration = 1;
        let result = initial.result();
        let error = initial.error();
        let large_interval_error = error;
        let table = ExtrapolationTable::initialise(&initial);
        let smallest_interval = initial.abs_interval_length();
        let extrapolate = false;
        let store = BinaryHeap::with_capacity(2 * self.max_iterations + 1);
        let error_kind = None;
        let initial_absolute_result = initial.result_abs();
        let positive_integrand = initial.positivity();
        let roundoff_count = 0;
        let roundoff_on_high_iteration_count = 0;
        let evaluations_per_integration = self.rule.evaluations() * self.evaluations_multiplier;

        heap.push(initial);

        Workspace {
            heap,
            iteration,
            result,
            error,
            large_interval_error,
            table,
            smallest_interval,
            extrapolate,
            store,
            error_kind,
            initial_absolute_result,
            positive_integrand,
            roundoff_count,
            roundoff_on_high_iteration_count,
            evaluations_per_integration,
        }
    }
}

/// A function defined over an infinite integration interval $x \in (-\infty,+\infty)$
///
/// This is a wrapper around a user supplied function implementing the [`Integrand`] trait which is
/// to be integrated over the infinite integration interval $x \in (-\infty,+\infty)$, i.e. in
/// integrals of the form,
/// $$
/// \int_{-\infty}^{+\infty} f(x) dx.
/// $$
/// This is achieved by transforming the integrand to be defined on the interval (0,1).
pub struct InfiniteInterval<I> {
    function: I,
}

impl<I: Integrand<Point = f64>> InfiniteInterval<I> {
    fn new(function: I) -> Self {
        Self { function }
    }

    fn transform_evaluate(&self, t: f64) -> I::Scalar {
        let x = (1.0 - t) / t;
        let y = self.function.evaluate(&x) + self.function.evaluate(&(-x));
        (y / t) / t
    }
}

impl<I: Integrand<Point = f64>> Integrand for InfiniteInterval<I> {
    type Point = I::Point;
    type Scalar = I::Scalar;
    fn evaluate(&self, x: &Self::Point) -> Self::Scalar {
        self.transform_evaluate(*x)
    }
}

impl<I> AdaptiveSingularity<InfiniteInterval<I>>
where
    I: Integrand<Point = f64>,
{
    /// Generate a new [`AdaptiveSingularity`] integrator for functions with possible integrable
    /// singularities on infinite intervals.
    ///
    ///
    /// Initialise an adaptive Gauss-Kronrod integrator with extrapolation for approximating
    /// integrals of the form,
    /// $$
    /// \int_{-\infty}^{+\infty} f(x) dx.
    /// $$
    /// with finite integration bounds $x \in (-\infty,+\infty)$.
    ///
    ///  Arguments:
    /// - `function`: A user supplied function to be integrated which is something implementing the
    /// [`Integrand`] trait.
    /// - `tolerance`: The tolerance requested by the user. Can be either an absolute tolerance
    /// or relative tolerance. Determines the exit condition of the integration routine, see
    /// [`Tolerance`].
    /// - `max_iterations`: The maximum number of iterations that the routine should use to try to
    /// satisfy the requested tolerance.
    ///
    /// # Errors
    /// Function can return an error if it receives bad user input. This is primarily related to
    /// using invalid values for the `tolerance`. The returned error is an [`InitialisationError`].
    /// See [`Tolerance`] and [`InitialisationError`] for more details.
    ///
    /// # Example
    ///
    ///
    /// Here we present a calculation of [Euler's constant] $\gamma$ ([`std::f64::consts::EULER_GAMMA`])
    /// using the integral representation:
    /// $$
    /// \gamma = -\int_{-\infty}^{+\infty} x e^{x - e^{x}} dx
    /// $$
    /// which has an infinite integration region $x\in(-\infty,+\infty)$. We use [`AdaptiveSingularity::infinite`]
    /// to approximate $\gamma$.
    ///```rust
    /// use rint::{Integrand, Tolerance};
    /// use rint::quadrature::AdaptiveSingularity;
    ///
    /// const GAMMA: f64 = std::f64::consts::EULER_GAMMA;
    /// const TOL: f64 = 1.0e-12;
    ///
    /// struct Gamma;
    ///
    /// impl Integrand for Gamma {
    ///     type Point = f64;
    ///     type Scalar = f64;
    ///
    ///     fn evaluate(&self, x: &Self::Point) -> Self::Scalar {
    ///         let exponent = x - x.exp();
    ///         -x * exponent.exp()
    ///     }
    /// }
    ///
    /// # use std::error::Error;
    /// # fn main() -> Result<(), Box<dyn Error>> {
    /// let gamma = Gamma;
    /// let tolerance = Tolerance::Relative(TOL);
    /// let integral = AdaptiveSingularity::infinite(gamma, tolerance, 1000)?
    ///     .integrate()?;
    ///
    /// let result = integral.result();
    /// let error = integral.error();
    /// let abs_actual_error = (GAMMA - result).abs();
    /// let tol = TOL * result.abs();
    /// let iters = integral.iterations();
    /// assert_eq!(iters, 11);
    /// assert!(abs_actual_error < error);
    /// assert!(error < tol);
    /// # Ok(())
    /// # }
    ///```
    ///
    pub fn infinite(
        function: I,
        tolerance: Tolerance,
        max_iterations: usize,
    ) -> Result<Self, InitialisationError> {
        let rule = Rule::gk15();
        let transformed = InfiniteInterval::new(function);
        let evaluations_multiplier = 2;
        Self::new_with_evaluations_multiplier(
            transformed,
            rule,
            Limits::new(0.0, 1.0),
            tolerance,
            max_iterations,
            evaluations_multiplier,
        )
    }
}

/// A function defined over a semi-infinite integration interval $x \in (b,+\infty)$
///
/// This is a wrapper around a user supplied function implementing the [`Integrand`] trait which is
/// to be integrated over the infinite integration interval $x \in (b,+\infty)$, i.e. in
/// integrals of the form,
/// $$
/// \int_{b}^{+\infty} f(x) dx.
/// $$
/// This is achieved by transforming the integrand to be defined on the interval (0,1).
pub struct SemiInfiniteIntervalPositive<I> {
    function: I,
    lower: f64,
}

impl<I: Integrand<Point = f64>> SemiInfiniteIntervalPositive<I> {
    fn new(function: I, lower: f64) -> Self {
        Self { function, lower }
    }

    fn transform_evaluate(&self, t: f64) -> I::Scalar {
        let x = self.lower + (1.0 - t) / t;
        let y = self.function.evaluate(&x);
        y / (t.powi(2))
    }
}

impl<I: Integrand<Point = f64>> Integrand for SemiInfiniteIntervalPositive<I> {
    type Point = I::Point;
    type Scalar = I::Scalar;
    fn evaluate(&self, x: &Self::Point) -> Self::Scalar {
        self.transform_evaluate(*x)
    }
}

impl<I> AdaptiveSingularity<SemiInfiniteIntervalPositive<I>>
where
    I: Integrand<Point = f64>,
{
    /// Generate a new [`AdaptiveSingularity`] integrator for functions with possible integrable
    /// singularities on semi-infinite intervals with an upper infinite limit.
    ///
    ///
    /// Initialise an adaptive Gauss-Kronrod integrator with extrapolation for approximating
    /// integrals of the form,
    /// $$
    /// \int_{b}^{+\infty} f(x) dx.
    /// $$
    /// with finite integration bounds $x \in (b,+\infty)$.
    ///
    ///  Arguments:
    /// - `function`: A user supplied function to be integrated which is something implementing the
    /// [`Integrand`] trait.
    /// - `lower`: the lower integration bound.
    /// - `tolerance`: The tolerance requested by the user. Can be either an absolute tolerance
    /// or relative tolerance. Determines the exit condition of the integration routine, see
    /// [`Tolerance`].
    /// - `max_iterations`: The maximum number of iterations that the routine should use to try to
    /// satisfy the requested tolerance.
    ///
    /// # Errors
    /// Function can return an error if it receives bad user input. This is primarily related to
    /// using invalid values for the `tolerance`. The returned error is an [`InitialisationError`].
    /// See [`Tolerance`] and [`InitialisationError`] for more details.
    ///
    /// # Example
    ///
    /// Here we present a calculation of [Catalan's constant] $G$ using the integral representation:
    /// $$
    /// G = \frac{\pi}{2} \int_{1}^{+\infty} \frac{(x^{4} - 6 x^{2} + 1) \ln \ln x}{(1+x^{2})^{3}} d x
    /// $$
    /// which has a semi-infinite integration region $x \in (1,+\infty)$ We use
    /// [`AdaptiveSingularity::semi_infinite_upper`] to approximate $G$.
    ///```rust
    /// use rint::{Integrand, Limits, Tolerance};
    /// use rint::quadrature::AdaptiveSingularity;
    ///
    /// const G: f64 = 0.915_965_594_177_219_015_054_603_514_932_384_110_774;
    /// const TOL: f64 = 1.0e-12;
    ///
    /// struct Catalan;
    ///
    /// impl Integrand for Catalan {
    ///     type Point = f64;
    ///     type Scalar = f64;
    ///
    ///     fn evaluate(&self, x: &Self::Point) -> Self::Scalar {
    ///         let FRAC_PI_2 = std::f64::consts::FRAC_PI_2;
    ///         let polynomial = x.powi(4) - 6.0 * x.powi(2) + 1.0;
    ///         let denominator = (1.0 + x.powi(2)).powi(3);
    ///         let lnlnx = x.ln().ln();
    ///         
    ///         FRAC_PI_2 * polynomial * lnlnx / denominator
    ///     }
    /// }
    ///
    /// # use std::error::Error;
    /// # fn main() -> Result<(), Box<dyn Error>> {
    /// let catalan = Catalan;
    /// let lower = 1.0;
    /// let tolerance = Tolerance::Relative(TOL);
    /// let integral = AdaptiveSingularity::semi_infinite_upper(catalan, lower, tolerance, 1000)?
    ///     .integrate()?;
    ///
    /// let result = integral.result();
    /// let error = integral.error();
    /// let abs_actual_error = (G - result).abs();
    /// let tol = TOL * result.abs();
    /// let iters = integral.iterations();
    /// assert_eq!(iters, 32);
    /// assert!(abs_actual_error < error);
    /// assert!(error < tol);
    /// # Ok(())
    /// # }
    ///```
    pub fn semi_infinite_upper(
        function: I,
        lower: f64,
        tolerance: Tolerance,
        max_iterations: usize,
    ) -> Result<Self, InitialisationError> {
        let rule = Rule::gk15();
        let transformed = SemiInfiniteIntervalPositive::new(function, lower);
        let evaluations_multiplier = 2;
        Self::new_with_evaluations_multiplier(
            transformed,
            rule,
            Limits::new(0.0, 1.0),
            tolerance,
            max_iterations,
            evaluations_multiplier,
        )
    }
}

/// A function defined over a semi-infinite integration interval $x \in (-\infty,b)$
///
/// This is a wrapper around a user supplied function implementing the [`Integrand`] trait which is
/// to be integrated over the infinite integration interval $x \in (-\infty,b)$, i.e. in
/// integrals of the form,
/// $$
/// \int_{-\infty}^{b} f(x) dx.
/// $$
/// This is achieved by transforming the integrand to be defined on the interval (0,1).
pub struct SemiInfiniteIntervalNegative<I> {
    function: I,
    upper: f64,
}

impl<I: Integrand<Point = f64>> SemiInfiniteIntervalNegative<I> {
    fn new(function: I, upper: f64) -> Self {
        Self { function, upper }
    }

    fn transform_evaluate(&self, t: f64) -> I::Scalar {
        let x = self.upper - (1.0 - t) / t;
        let y = self.function.evaluate(&x);
        y / (t.powi(2))
    }
}

impl<I: Integrand<Point = f64>> Integrand for SemiInfiniteIntervalNegative<I> {
    type Point = I::Point;
    type Scalar = I::Scalar;
    fn evaluate(&self, x: &Self::Point) -> Self::Scalar {
        self.transform_evaluate(*x)
    }
}

impl<I> AdaptiveSingularity<SemiInfiniteIntervalNegative<I>>
where
    I: Integrand<Point = f64>,
{
    /// Generate a new [`AdaptiveSingularity`] integrator for functions with possible integrable
    /// singularities on semi-infinite intervals with an upper infinite limit.
    ///
    ///
    /// Initialise an adaptive Gauss-Kronrod integrator with extrapolation for approximating
    /// integrals of the form,
    /// $$
    /// \int_{-\infty}^{a} f(x) dx.
    /// $$
    /// with finite integration bounds $x \in (-\infty,a)$.
    ///
    ///  Arguments:
    /// - `function`: A user supplied function to be integrated which is something implementing the
    /// [`Integrand`] trait.
    /// - `upper`: the upper integration bound.
    /// - `tolerance`: The tolerance requested by the user. Can be either an absolute tolerance
    /// or relative tolerance. Determines the exit condition of the integration routine, see
    /// [`Tolerance`].
    /// - `max_iterations`: The maximum number of iterations that the routine should use to try to
    /// satisfy the requested tolerance.
    ///
    /// # Errors
    /// Function can return an error if it receives bad user input. This is primarily related to
    /// using invalid values for the `tolerance`. The returned error is an [`InitialisationError`].
    /// See [`Tolerance`] and [`InitialisationError`] for more details.
    ///
    /// # Example
    ///
    /// Here we present a calculation of [Catalan's constant] $G$ using the integral representation:
    /// $$
    /// G = \frac{\pi}{2} \int^{-1}_{-\infty} \frac{(x^{4} - 6 x^{2} + 1) \ln \ln (-x)}{(1+x^{2})^{3}} d x
    /// $$
    /// which has a semi-infinite integration region $x \in (1,+\infty)$. We use
    /// [`AdaptiveSingularity::semi_infinite_upper`] to approximate $G$.
    ///```rust
    /// use rint::{Integrand, Limits, Tolerance};
    /// use rint::quadrature::AdaptiveSingularity;
    ///
    /// const G: f64 = 0.915_965_594_177_219_015_054_603_514_932_384_110_774;
    /// const TOL: f64 = 1.0e-12;
    ///
    /// struct Catalan;
    ///
    /// impl Integrand for Catalan {
    ///     type Point = f64;
    ///     type Scalar = f64;
    ///
    ///     fn evaluate(&self, x: &Self::Point) -> Self::Scalar {
    ///         let FRAC_PI_2 = std::f64::consts::FRAC_PI_2;
    ///         let polynomial = x.powi(4) - 6.0 * x.powi(2) + 1.0;
    ///         let denominator = (1.0 + x.powi(2)).powi(3);
    ///         let lnlnx = (-x).ln().ln();
    ///         
    ///         FRAC_PI_2 * polynomial * lnlnx / denominator
    ///     }
    /// }
    /// # use std::error::Error;
    /// # fn main() -> Result<(), Box<dyn Error>> {
    /// let catalan = Catalan;
    /// let upper = -1.0;
    /// let tolerance = Tolerance::Relative(TOL);
    /// let integral = AdaptiveSingularity::semi_infinite_lower(catalan, upper, tolerance, 1000)?
    ///     .integrate()?;
    ///
    /// let result = integral.result();
    /// let error = integral.error();
    /// let abs_actual_error = (G - result).abs();
    /// let tol = TOL * result.abs();
    /// let iters = integral.iterations();
    /// assert_eq!(iters, 32);
    /// assert!(abs_actual_error < error);
    /// assert!(error < tol);
    /// # Ok(())
    /// # }
    ///```
    pub fn semi_infinite_lower(
        function: I,
        upper: f64,
        tolerance: Tolerance,
        max_iterations: usize,
    ) -> Result<Self, InitialisationError> {
        let rule = Rule::gk15();
        let transformed = SemiInfiniteIntervalNegative::new(function, upper);
        let evaluations_multiplier = 2;
        Self::new_with_evaluations_multiplier(
            transformed,
            rule,
            Limits::new(0.0, 1.0),
            tolerance,
            max_iterations,
            evaluations_multiplier,
        )
    }
}

/// The integration workspace.
///
/// We use a workspace for maintaining the state of the adaptive integration routines. The `heap`
/// is a [`BinaryHeap`] storing [`Region`]s which have been integrated, and upon using
/// [`Workspace::pop`] returns the region with the highest `error` approximation which is the next
/// region to be bisected and re-integrated.
/// As well as this, the [`AdaptiveSingularity`] routine needs to keep track of the
/// [`ExtrapolationTable`] and the `store`
/// Also tracked are the current estimates of `result` and
/// `error` across the entire integration region, the number of `iteration`s, and counts of times
/// that roundoff errors may have occurred, used as a diagnostic for terminating the integration
/// early if problematic roundoff behaviour is observed.
struct Workspace<T> {
    heap: BinaryHeap<Region<T>>,
    iteration: usize,
    result: T,
    error: f64,
    large_interval_error: f64,
    table: ExtrapolationTable<T>,
    smallest_interval: f64,
    extrapolate: bool,
    store: BinaryHeap<Region<T>>,
    error_kind: Option<IntegrationErrorKind>,
    initial_absolute_result: f64,
    positive_integrand: bool,
    roundoff_count: usize,
    roundoff_on_high_iteration_count: usize,
    evaluations_per_integration: usize,
}

impl<T: ScalarF64> Workspace<T> {
    /// Get the next [`Region`] in the [`BinaryHeap`] with the largest `error` estimate and update
    /// the workspace iteration number.
    fn retrieve_largest_error(&mut self) -> Result<Region<T>, IntegrationError<T>> {
        self.iteration += 1;
        if let Some(previous) = self.pop() {
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

    fn peek(&self) -> Option<&Region<T>> {
        self.heap.peek()
    }

    fn store(&mut self, integral: Region<T>) {
        self.store.push(integral);
    }

    fn sum_results(&self) -> T {
        let initial = T::zero();
        let store = self.store.iter().fold(initial, |a, v| a + v.result());
        let heap = self.heap.iter().fold(initial, |a, v| a + v.result());
        store + heap
    }

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

        if lower.result_asc().to_bits() != lower.error().to_bits()
            && upper.result_asc().to_bits() != upper.error().to_bits()
        {
            let delta = (prev_result - new_result).abs();

            if delta <= 1e-5 * new_result.abs() && new_error >= 0.99 * prev_error {
                if self.extrapolate {
                    self.table.roundoff_count += 1;
                } else {
                    self.roundoff_count += 1;
                }
            } else if self.iteration > 10 && new_error > prev_error {
                self.roundoff_on_high_iteration_count += 1;
            }
        }

        self.result = (self.result + new_result) - prev_result;
        self.error = (self.error + new_error) - prev_error;

        (self.result, self.error)
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

    fn prepare_next_iteration(&mut self) {
        self.extrapolate = false;
        self.large_interval_error = self.error;
        self.heap.append(&mut self.store);
        self.smallest_interval *= 0.5;
    }

    fn compute_result(self) -> Result<IntegralEstimate<T>, IntegrationError<T>> {
        let output = self.integral_estimate();

        if let Some(kind) = self.error_kind {
            Err(IntegrationError::new(output, kind))
        } else {
            Ok(output)
        }
    }

    fn compute_extrapolated_result(self) -> Result<IntegralEstimate<T>, IntegrationError<T>> {
        let evaluations = (2 * self.iteration - 1) * self.evaluations_per_integration;
        let output = self.table.integral_estimate(self.iteration, evaluations);

        if let Some(kind) = self.error_kind {
            Err(IntegrationError::new(output, kind))
        } else {
            Ok(output)
        }
    }

    fn check_error_and_compute(mut self) -> Result<IntegralEstimate<T>, IntegrationError<T>> {
        if (self.table.error - f64::MAX).abs() < f64::EPSILON {
            return self.compute_result();
        }

        if self.is_err() || self.table.is_err() {
            if self.table.is_err() {
                self.table.error += self.table.correction;
            }

            if !self.is_err() {
                self.set_error_kind(IntegrationErrorKind::DoesNotConverge);
            }

            let zero = T::zero();

            if self.table.result != zero && self.result != zero {
                if self.table.error / self.table.result.abs() > self.error / self.result.abs() {
                    return self.compute_result();
                }
            } else if self.table.error > self.error {
                return self.compute_result();
            } else if self.result == zero {
                return self.compute_extrapolated_result();
            }
        }

        let max_area = f64::max(self.table.result.abs(), self.result.abs());

        if !self.positive_integrand && max_area < 0.01 * self.initial_absolute_result {
            return self.compute_extrapolated_result();
        }

        let ratio = self.table.result.abs() / self.result.abs();

        if !(0.01..=100.0).contains(&ratio) || self.error > self.result.abs() {
            self.set_error_kind(IntegrationErrorKind::DivergentOrSlowlyConverging);
        }

        self.compute_extrapolated_result()
    }

    const fn check_roundoff(&mut self) {
        if self.roundoff_count + self.table.roundoff_count >= 10
            || self.roundoff_on_high_iteration_count >= 20
        {
            self.set_error_kind(IntegrationErrorKind::RoundoffErrorDetected);
        }

        if self.table.roundoff_count >= 5 {
            self.table.error_detected = true;
        }
    }

    const fn check_convergence(&mut self) {
        if self.table.ktmin > 5 && self.table.error < 0.001 * self.error {
            self.set_error_kind(IntegrationErrorKind::DoesNotConverge);
        }
    }

    fn remaining_large_intervals(&mut self) -> bool {
        while let Some(next) = self.peek() {
            if next.abs_interval_length() > self.smallest_interval {
                return true;
            }
            if let Some(next) = self.pop() {
                self.store(next);
            }
        }
        false
    }

    const fn is_err(&self) -> bool {
        self.error_kind.is_some()
    }

    const fn set_error_kind(&mut self, kind: IntegrationErrorKind) {
        self.error_kind = Some(kind);
    }

    fn update(&mut self, lower: Region<T>, upper: Region<T>) {
        let lower_limit = lower.limits().lower();
        let upper_limit = upper.limits().upper();
        let limits = Limits::new(lower_limit, upper_limit);

        self.check_roundoff();

        if limits.subinterval_too_small() {
            self.set_error_kind(IntegrationErrorKind::BadIntegrandBehaviour(limits));
        }

        self.push(lower);
        self.push(upper);
    }

    const fn update_large_interval_error(
        &mut self,
        iteration_interval: f64,
        previous_error: f64,
        iteration_error: f64,
    ) {
        self.large_interval_error += -previous_error;

        let half_iteration_interval = iteration_interval * 0.5;
        if half_iteration_interval > self.smallest_interval {
            self.large_interval_error += iteration_error;
        }
    }

    const fn update_table_values(&mut self, ext_result: T, ext_error: f64, tolerance: f64) {
        self.table.ktmin = 0;
        self.table.error = ext_error;
        self.table.result = ext_result;
        self.table.correction = self.large_interval_error;
        self.table.tolerance = tolerance;
    }
}

struct ExtrapolationTable<T> {
    count: usize,
    results: [T; 52],
    cached: Cached<T>,
    error_detected: bool,
    result: T,
    error: f64,
    tolerance: f64,
    ktmin: usize,
    roundoff_count: usize,
    correction: f64,
}

#[derive(Debug, PartialEq)]
enum Cached<T> {
    Empty,
    One(T),
    Two(T, T),
    Full(T, T, T),
}

impl<T: ScalarF64> ExtrapolationTable<T> {
    fn new() -> Self {
        let zero = T::zero();
        let count = 0;
        let results = [zero; 52];
        let cached: Cached<T> = Cached::Empty;
        let error_detected = false;
        let result = zero;
        let error = f64::MAX;
        let tolerance = f64::MAX;
        let ktmin = 0;
        let roundoff_count = 0;
        let correction = 0.0;

        Self {
            count,
            results,
            cached,
            error_detected,
            result,
            error,
            tolerance,
            ktmin,
            roundoff_count,
            correction,
        }
    }

    fn initialise(initial: &Region<T>) -> Self {
        let mut table = Self::new();
        table.result = initial.result();
        table.append_table(initial.result());
        table
    }

    const fn is_err(&self) -> bool {
        self.error_detected
    }

    const fn cache(&mut self, value: T) -> &mut Self {
        match self.cached {
            Cached::Full(_, a, b) | Cached::Two(a, b) => {
                self.cached = Cached::Full(a, b, value);
            }
            Cached::One(a) => {
                self.cached = Cached::Two(a, value);
            }
            Cached::Empty => {
                self.cached = Cached::One(value);
            }
        }
        self
    }

    const fn append_table(&mut self, result: T) -> &mut Self {
        self.results[self.count] = result;
        self.count += 1;
        self
    }

    // Adapted directly from GSL
    #[must_use]
    fn extrapolate(&mut self, value: T) -> (T, f64) {
        self.append_table(value);
        let n_current = self.count - 1;

        let current = self.results[n_current];

        let mut absolute = f64::MAX;
        let mut relative = 5.0 * f64::EPSILON * current.abs();

        let new_element = n_current / 2;
        let n_original = n_current;
        let mut n_final = n_current;

        let mut result = current;
        let mut abserr = f64::MAX;

        if n_current < 2 {
            result = current;
            abserr = f64::max(absolute, relative);
            return (result, abserr);
        }

        self.results[n_current + 2] = self.results[n_current];
        self.results[n_current] = T::max_value();

        for i in 0..new_element {
            let mut res = self.results[n_current - 2 * i + 2];
            let e0 = self.results[n_current - 2 * i - 2];
            let e1 = self.results[n_current - 2 * i - 1];
            let e2 = res;

            let e1abs = e1.abs();
            let delta2 = e2 - e1;
            let err2 = delta2.abs();
            let tol2 = f64::max(e2.abs(), e1abs) * f64::EPSILON;
            let delta3 = e1 - e0;
            let err3 = delta3.abs();
            let tol3 = f64::max(e1abs, e0.abs()) * f64::EPSILON;

            if err2 <= tol2 && err3 <= tol3 {
                // If e0, e1 and e2 are equal to within machine accuracy,
                // convergence is assumed.
                result = res;
                absolute = err2 + err3;
                relative = 5.0 * f64::EPSILON * res.abs();
                abserr = f64::max(absolute, relative);
                return (result, abserr);
            }

            let e3 = self.results[n_current - 2 * i];
            self.results[n_current - 2 * i] = e1;
            let delta1 = e1 - e3;
            let err1 = delta1.abs();
            let tol1 = f64::max(e1abs, e3.abs()) * f64::EPSILON;

            // If two elements are very close to each other, omit a part
            // of the table by adjusting the value of n.
            if err1 <= tol1 || err2 <= tol2 || err3 <= tol3 {
                n_final = 2 * i;
                break;
            }

            let delta1c = delta1.conj();
            let delta2c = delta2.conj();
            let delta3c = delta3.conj();

            let delta1modsq = delta1 * delta1c;
            let delta2modsq = delta2 * delta2c;
            let delta3modsq = delta3 * delta3c;

            let ss = ((delta1c / delta1modsq) + (delta2c / delta2modsq)) - (delta3c / delta3modsq);

            // Test to detect irregular behaviour in the table, and eventually
            // omit a part of the table by adjusting the value of n.
            if (e1 * ss).abs() <= 0.0001 {
                n_final = 2 * i;
                break;
            }

            let ssc = ss.conj();
            let ssmodsq = ss * ssc;

            // Compute a new element and eventually adjust the value of result.
            res = e1 + (ssc / ssmodsq);
            self.results[n_current - 2 * i] = res;

            let error = err2 + (res - e2).abs() + err3;

            if error <= abserr {
                abserr = error;
                result = res;
            }
        }

        let limexp = 50 - 1;
        if n_final == limexp {
            n_final = 2 * (limexp / 2);
        }

        self.shift_table(n_original, n_final, new_element);
        self.count = n_final + 1;

        if let Cached::Full(a, b, c) = self.cached {
            abserr = (result - c).abs() + (result - b).abs() + (result - a).abs();
        } else {
            abserr = f64::MAX;
        }

        self.cache(result);
        abserr = f64::max(abserr, 5.0 * f64::EPSILON * result.abs());

        self.ktmin += 1;

        (result, abserr)
    }

    fn shift_table(&mut self, n_original: usize, n_final: usize, new_element: usize) -> &mut Self {
        if n_original % 2 == 1 {
            for i in 0..=new_element {
                self.results[1 + i * 2] = self.results[i * 2 + 3];
            }
        } else {
            for i in 0..=new_element {
                self.results[i * 2] = self.results[i * 2 + 2];
            }
        }

        if n_original != n_final {
            for i in 0..=n_final {
                self.results[i] = self.results[n_original - n_final + i];
            }
        }
        self
    }

    fn integral_estimate(&self, iterations: usize, evaluations: usize) -> IntegralEstimate<T> {
        let result = self.result;
        let error = self.error;
        IntegralEstimate::new()
            .with_result(result)
            .with_error(error)
            .with_iterations(iterations)
            .with_evaluations(evaluations)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_swapping_cached() {
        let mut table = ExtrapolationTable::new();
        let value0 = 1.0;
        let value1 = 10.0;
        let value2 = 100.0;
        let value = 999.0;

        assert_eq!(table.cached, Cached::Empty);

        table.cache(value0);

        assert_eq!(table.cached, Cached::One(1.0));

        table.cache(value1);

        assert_eq!(table.cached, Cached::Two(1.0, 10.0));

        table.cache(value2);

        assert_eq!(table.cached, Cached::Full(1.0, 10.0, 100.0));

        table.cache(value);

        assert_eq!(table.cached, Cached::Full(10.0, 100.0, 999.0));
    }

    #[test]
    fn test_intitialising() {
        let error = 0.0;
        let result = 0.0;
        let result_abs = 0.0;
        let result_asc = 0.0;
        let limits = Limits::new(0.0, 0.0);
        let value = Region {
            result,
            error,
            result_abs,
            result_asc,
            limits,
        };

        let table = ExtrapolationTable::initialise(&value);

        assert!(table.count == 1);
        assert!((table.results[0] - value.result()) < 1e-10);
    }
}
