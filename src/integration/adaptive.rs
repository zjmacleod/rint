use std::collections::binary_heap::{BinaryHeap, IntoIter};

use crate::integration::basic::BasicInternal;
use crate::integration::{ErrorBound, GaussKronrodBasic};
use crate::rule::Rule;
use crate::Integrand;

/// The value of a function evaluated with Gauss-Kronrod integration and associated error
/// estimation.
#[derive(Debug)]
pub struct Adaptive {
    result: f64,
    error: f64,
    iterations: usize,
}

impl Adaptive {
    pub(crate) fn new(result: f64, error: f64, iterations: usize) -> Self {
        Self {
            result,
            error,
            iterations,
        }
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

    pub(crate) fn from_basic(basic: &BasicInternal, iterations: usize) -> Self {
        let result = basic.result();
        let error = basic.error();
        Self::new(result, error, iterations)
    }

    pub(crate) fn empty() -> Self {
        Self {
            result: f64::NAN,
            error: f64::NAN,
            iterations: 0,
        }
    }
}

/// An integral to be evaluated with an adaptive Gauss-Kronrod quadrature.
///
/// The user constructs a `function` implementing [`Integrand`], provides `upper`
/// and `lower` integration limits, and provides an `error_bound`, which can be
/// [`ErrorBound::Absolute`] to work to a specified absolute error,
/// [`ErrorBound::Relative`] to work to a specified relative error,
/// or [`ErrorBound::Either`] to return a result as soon as _either_ the relative
/// or absolute error bound has been satisfied.
pub struct GaussKronrodAdaptive<'a, I, R>
where
    I: Integrand,
    R: Rule,
{
    lower: f64,
    upper: f64,
    error_bound: ErrorBound,
    rule: R,
    function: &'a I,
    max_iterations: usize,
}

impl<'a, I, R> GaussKronrodAdaptive<'a, I, R>
where
    I: Integrand,
    R: Rule,
{
    /// Create a new [`GaussKronrodAdaptive`].
    ///
    /// The user defines a `function` which is a `struct` implementing the
    /// [`Integrand`] trait, and integration limis `upper` and `lower`.
    ///
    /// # Errors
    /// Function will return an error if the user provided `ErrorBound` does not satisfy the
    /// following constraints:
    ///     - `ErrorBound::Absolute(v) where v > 0.0`,
    ///     - `ErrorBound::Relative(v) where v > 50.0 * f64::EPSILON`,
    ///     - `ErrorBound::Either { absolute, relative } where absolute > 0.0 and relative > 50.0 * f64::EPSILON`.
    pub fn new(
        lower: f64,
        upper: f64,
        error_bound: ErrorBound,
        rule: R,
        function: &'a I,
        max_iterations: usize,
    ) -> Result<Self, Error> {
        match error_bound {
            ErrorBound::Absolute(v) => {
                if v <= 0.0 {
                    let output = Adaptive::new(0.0, 0.0, 0);
                    let kind = Kind::RelativeBoundNegativeOrZero;
                    return Err(Error::new(kind, output));
                }
            }
            ErrorBound::Relative(v) => {
                if v < 50.0 * f64::EPSILON {
                    let output = Adaptive::new(0.0, 0.0, 0);
                    let kind = Kind::AbsoluteBoundTooSmall;
                    return Err(Error::new(kind, output));
                }
            }
            ErrorBound::Either { absolute, relative } => {
                if absolute <= 0.0 && relative < 50.0 * f64::EPSILON {
                    let output = Adaptive::new(0.0, 0.0, 0);
                    let kind = Kind::InvalidTolerance;
                    return Err(Error::new(kind, output));
                }
            }
        }
        Ok(Self {
            lower,
            upper,
            error_bound,
            rule,
            function,
            max_iterations,
        })
    }

    /// Integrate the function and return a [`Adaptive`] integration result.
    ///
    /// # Errors
    /// Integration can fail if user suplied tolerance cannot be achieved within the maximum number
    /// of iterations.
    pub fn integrate(&self) -> Result<Adaptive, Error> {
        let initial = GaussKronrodBasic::new(self.lower, self.upper, self.rule, self.function)
            .integrate_internal();

        if let Some(output) = self.check_initial_integration(&initial)? {
            return Ok(output);
        }

        let mut workspace = Workspace::new(self.max_iterations, initial);

        while workspace.iteration < self.max_iterations {
            let [previous, lower, upper] =
                workspace.bisect_largest_error(self.function, self.rule)?;

            let prev_area = previous.result();
            let prev_error = previous.error();

            let new_area = lower.result() + upper.result();
            let new_error = lower.error() + upper.error();

            let area = workspace.update_current_area(new_area - prev_area);
            let error = workspace.update_error(new_error - prev_error);

            if lower.result_asc().to_bits() != lower.error().to_bits()
                && upper.result_asc().to_bits() != upper.error().to_bits()
            {
                let delta = (prev_area - new_area).abs();

                if delta <= 1e-5 * new_area.abs() && new_error >= 0.99 * prev_error {
                    workspace.roundoff_type1 += 1;
                }
                if workspace.iteration >= 10 && new_error >= prev_error {
                    workspace.roundoff_type2 += 1;
                }
            }

            let iteration_tolerance = self.error_bound.tolerance(&area);

            if error > iteration_tolerance {
                if workspace.roundoff_type1 >= 6 || workspace.roundoff_type2 >= 20 {
                    let output = workspace.integral_estimate();
                    let kind = Kind::FailedToReachToleranceRoundoff;
                    return Err(Error::new(kind, output));
                }

                let lower_limit = previous.lower();
                let upper_limit = previous.upper();
                let midpoint = (upper_limit + lower_limit) * 0.5;

                if subinterval_too_small(lower_limit, midpoint, upper_limit) {
                    let output = workspace.integral_estimate();
                    let kind = Kind::PossibleSingularity {
                        lower: lower_limit,
                        upper: upper_limit,
                    };
                    return Err(Error::new(kind, output));
                }
            }

            workspace.push(lower);
            workspace.push(upper);

            if error < iteration_tolerance {
                break;
            }
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

    /// Return the value of the `upper` integration limit.
    pub fn upper(&self) -> f64 {
        self.upper
    }

    /// Return the value of the `lower` integration limit.
    pub fn lower(&self) -> f64 {
        self.lower
    }

    fn check_initial_integration(
        &self,
        initial: &BasicInternal,
    ) -> Result<Option<Adaptive>, Error> {
        let tolerance = self.error_bound.tolerance(&initial.result());
        let roundoff = initial.roundoff();

        if initial.error() <= roundoff && initial.error() > tolerance {
            let output = Adaptive::from_basic(initial, 1);
            let kind = Kind::FailedToReachToleranceRoundoff;

            Err(Error::new(kind, output))
        } else if (initial.error() <= tolerance
            && initial.error().to_bits() != initial.result_asc().to_bits())
            || initial.error() == 0.0
        {
            let output = Adaptive::from_basic(initial, 1);

            Ok(Some(output))
        } else if self.max_iterations == 1 {
            let output = Adaptive::from_basic(initial, 1);
            let kind = Kind::MaximumIterationsReached;

            Err(Error::new(kind, output))
        } else {
            Ok(None)
        }
    }
}

struct Workspace {
    heap: BinaryHeap<BasicInternal>,
    iteration: usize,
    roundoff_type1: usize,
    roundoff_type2: usize,
    current_area: f64,
    error: f64,
}

impl Workspace {
    fn new(max_iterations: usize, initial: BasicInternal) -> Self {
        let mut heap = BinaryHeap::with_capacity(2 * max_iterations + 1);
        let current_area = initial.result();
        let error = initial.error();
        heap.push(initial);
        Self {
            heap,
            iteration: 1,
            roundoff_type1: 0,
            roundoff_type2: 0,
            current_area,
            error,
        }
    }

    fn bisect_largest_error<F: Integrand, R: Rule>(
        &mut self,
        function: &F,
        rule: R,
    ) -> Result<[BasicInternal; 3], Error> {
        let Some(previous) = self.pop() else {
                let output = Adaptive::empty();
                let kind = Kind::UninitialisedWorkspace;
                return Err(Error::new(kind, output));
            };

        self.iteration += 1;

        let lower = previous.lower();
        let upper = previous.upper();
        let mid = (upper + lower) * 0.5;

        let lower_integral =
            GaussKronrodBasic::new(lower, mid, rule, function).integrate_internal();

        let upper_integral =
            GaussKronrodBasic::new(mid, upper, rule, function).integrate_internal();

        Ok([previous, lower_integral, upper_integral])
    }

    fn sum_results(self) -> f64 {
        self.into_iter().fold(0.0f64, |a, v| a + v.result())
    }

    fn update_current_area(&mut self, diff: f64) -> f64 {
        self.current_area += diff;
        self.current_area
    }

    fn update_error(&mut self, diff: f64) -> f64 {
        self.error += diff;
        self.error
    }

    fn integral_estimate(self) -> Adaptive {
        let error = self.error;
        let iterations = self.iteration;
        let result = self.sum_results();
        Adaptive::new(result, error, iterations)
    }
}

impl std::ops::Deref for Workspace {
    type Target = BinaryHeap<BasicInternal>;
    fn deref(&self) -> &Self::Target {
        &self.heap
    }
}

impl std::ops::DerefMut for Workspace {
    fn deref_mut(&mut self) -> &mut BinaryHeap<BasicInternal> {
        &mut self.heap
    }
}

impl IntoIterator for Workspace {
    type Item = BasicInternal;
    type IntoIter = IntoIter<BasicInternal>;

    fn into_iter(self) -> Self::IntoIter {
        self.heap.into_iter()
    }
}

#[derive(Debug)]
pub enum Kind {
    RelativeBoundNegativeOrZero,
    AbsoluteBoundTooSmall,
    InvalidTolerance,
    FailedToReachToleranceRoundoff,
    MaximumIterationsReached,
    PossibleSingularity { lower: f64, upper: f64 },
    UninitialisedWorkspace,
}

#[derive(Debug)]
pub struct Error {
    kind: Kind,
    integral: Adaptive,
}

impl Error {
    pub(crate) fn new(kind: Kind, integral: Adaptive) -> Self {
        Self { kind, integral }
    }

    #[must_use]
    pub fn kind(&self) -> &Kind {
        &self.kind
    }

    #[must_use]
    pub fn integral(&self) -> &Adaptive {
        &self.integral
    }

    #[must_use]
    pub fn result(&self) -> f64 {
        self.integral.result()
    }

    #[must_use]
    pub fn error(&self) -> f64 {
        self.integral.error()
    }

    #[must_use]
    pub fn iterations(&self) -> usize {
        self.integral.iterations()
    }
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.kind() {
            Kind::RelativeBoundNegativeOrZero => {
                write!(
                    f,
                    "Invalid error bound: relative bound must be non-zero and positive."
                )
            }
            Kind::AbsoluteBoundTooSmall => {
                write!(
                    f,
                    "Invalid error bound: absolute bound must be larger than 50.0 * f64::EPSILON."
                )
            }
            Kind::InvalidTolerance => {
                write!(f, "Invalid tolerance: relative bound must be non-zero and positive and absolute bound must be larger than 50.0 * f64::EPSILON.")
            }
            Kind::FailedToReachToleranceRoundoff => {
                let result = self.result();
                let error = self.error();
                write!(f, "Failed to reach tolerance due to roundoff error. Result: {result}, error: {error}.")
            }
            Kind::MaximumIterationsReached => {
                let result = self.result();
                let error = self.error();
                write!(
                    f,
                    "Maximum number of iterations reached. Result: {result}, error: {error}."
                )
            }
            Kind::PossibleSingularity { lower, upper } => {
                let result = self.result();
                let error = self.error();
                write!(f, "Possible singularity detected between ({lower},{upper}). Result: {result}, error: {error}.")
            }
            Kind::UninitialisedWorkspace => {
                write!(f, "The integration Workspace was not properly initialised.")
            }
        }
    }
}

impl std::error::Error for Error {}

pub(crate) fn subinterval_too_small(lower: f64, midpoint: f64, upper: f64) -> bool {
    let eps = f64::EPSILON;
    let min = f64::MIN_POSITIVE;

    let tmp = (1.0 + 100.0 * eps) * (midpoint.abs() + 1000.0 * min);

    lower.abs() <= tmp && upper.abs() <= tmp
}
