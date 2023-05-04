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
                    let output = Adaptive::empty();
                    let kind = Kind::RelativeBoundNegativeOrZero;
                    return Err(Error::new(kind, output));
                }
            }
            ErrorBound::Relative(v) => {
                if v < 50.0 * f64::EPSILON {
                    let output = Adaptive::empty();
                    let kind = Kind::AbsoluteBoundTooSmall;
                    return Err(Error::new(kind, output));
                }
            }
            ErrorBound::Either { absolute, relative } => {
                if absolute <= 0.0 && relative < 50.0 * f64::EPSILON {
                    let output = Adaptive::empty();
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
            let previous = workspace.retrieve_largest_error()?;

            let [lower, upper] = previous.bisect(self.function, self.rule);
            workspace.update_limits(lower.lower(), upper.upper());

            let (result, error) = workspace.iteration_result_error(&previous, &lower, &upper);

            let iteration_tolerance = self.error_bound.tolerance(&result);

            if error > iteration_tolerance {
                workspace.check_roundoff()?;
                workspace.check_singularity()?;
            }

            workspace.push(lower);
            workspace.push(upper);

            if error <= iteration_tolerance {
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
    result: f64,
    error: f64,
    lower_limit: f64,
    upper_limit: f64,
}

impl Workspace {
    fn new(max_iterations: usize, initial: BasicInternal) -> Self {
        let mut heap = BinaryHeap::with_capacity(2 * max_iterations + 1);

        let result = initial.result();
        let error = initial.error();
        let lower_limit = initial.lower();
        let upper_limit = initial.upper();

        heap.push(initial);

        Self {
            heap,
            iteration: 1,
            roundoff_type1: 0,
            roundoff_type2: 0,
            result,
            error,
            lower_limit,
            upper_limit,
        }
    }

    fn retrieve_largest_error(&mut self) -> Result<BasicInternal, Error> {
        self.iteration += 1;
        if let Some(previous) = self.pop() {
            Ok(previous)
        } else {
            let output = Adaptive::empty();
            let kind = Kind::UninitialisedWorkspace;
            Err(Error::new(kind, output))
        }
    }

    fn update_limits(&mut self, lower: f64, upper: f64) {
        self.lower_limit = lower;
        self.upper_limit = upper;
    }

    fn update_result(&mut self, diff: f64) -> f64 {
        self.result += diff;
        self.result
    }

    fn update_error(&mut self, diff: f64) -> f64 {
        self.error += diff;
        self.error
    }

    fn iteration_result_error(
        &mut self,
        previous: &BasicInternal,
        lower: &BasicInternal,
        upper: &BasicInternal,
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
                self.roundoff_type1 += 1;
            }
            if self.iteration >= 10 && new_error >= prev_error {
                self.roundoff_type2 += 1;
            }
        }

        let result = self.update_result(new_result - prev_result);
        let error = self.update_error(new_error - prev_error);

        (result, error)
    }

    fn check_roundoff(&self) -> Result<(), Error> {
        if self.roundoff_type1 >= 6 || self.roundoff_type2 >= 20 {
            let output = self.integral_estimate();
            let kind = Kind::FailedToReachToleranceRoundoff;
            Err(Error::new(kind, output))
        } else {
            Ok(())
        }
    }

    fn check_singularity(&self) -> Result<(), Error> {
        let lower_limit = self.lower_limit;
        let upper_limit = self.upper_limit;
        let midpoint = (lower_limit + upper_limit) * 0.5;
        if subinterval_too_small(lower_limit, midpoint, upper_limit) {
            let output = self.integral_estimate();
            let kind = Kind::PossibleSingularity {
                lower: lower_limit,
                upper: upper_limit,
            };
            Err(Error::new(kind, output))
        } else {
            Ok(())
        }
    }

    fn sum_results(&self) -> f64 {
        self.iter().fold(0.0f64, |a, v| a + v.result())
    }

    fn integral_estimate(&self) -> Adaptive {
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
