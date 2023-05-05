use crate::integration::basic::BasicInternal;
use crate::integration::workspace::Workspace;
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

    fn roundoff(result_abs: f64) -> f64 {
        50.0 * f64::EPSILON * result_abs
    }

    pub(crate) fn check_initial_integration(
        &self,
        initial: &BasicInternal,
    ) -> Result<Option<Adaptive>, Error> {
        let tolerance = self.error_bound.tolerance(initial.result());
        let roundoff = Self::roundoff(initial.result_abs());

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

        let mut workspace = Workspace::adaptive(self.max_iterations, initial);

        while workspace.iteration() < self.max_iterations {
            let previous = workspace.retrieve_largest_error()?;

            let [lower, upper] = previous.bisect(self.function, self.rule);
            workspace.update_limits(lower.lower(), upper.upper());

            let (result, error) = workspace.iteration_result_error(&previous, &lower, &upper);

            let iteration_tolerance = self.error_bound.tolerance(result);

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

        let final_iteration = workspace.iteration();
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
