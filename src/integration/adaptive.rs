use std::collections::BinaryHeap;

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
}

/// An integral to be evaluated with an adaptive Gauss-Kronrod quadrature.
///
/// The user constructs a `function` implementing [`Integrand`], provides `upper`
/// and `lower` integration limits, and provides an `error_bound`, which can be
/// [`ErrorBound::Absolute`] to work to a specified absolute error,
/// [`ErrorBound::Relative`] to work to a specified relative error,
/// or [`ErrorBound::Either`] to return a result as soon as _either_ the relative
/// or absolute error bound has been satisfied.
///
/// This routine successively applies the Gauss-Kronrod 10-point, 21-point, 43-point and 87-point
/// rules to a numerical integral of the `function` between `(lower, upper)` until an estimate has
/// been achieved within the desired error bound.
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
                    let partial_result = Adaptive::new(0.0, 0.0, 0);
                    let kind = Kind::RelativeBoundNegativeOrZero;
                    return Err(Error::new(kind, partial_result));
                }
            }
            ErrorBound::Relative(v) => {
                if v < 50.0 * f64::EPSILON {
                    let partial_result = Adaptive::new(0.0, 0.0, 0);
                    let kind = Kind::AbsoluteBoundTooSmall;
                    return Err(Error::new(kind, partial_result));
                }
            }
            ErrorBound::Either { absolute, relative } => {
                if absolute <= 0.0 && relative < 50.0 * f64::EPSILON {
                    let partial_result = Adaptive::new(0.0, 0.0, 0);
                    let kind = Kind::InvalidTolerance;
                    return Err(Error::new(kind, partial_result));
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

    /// Integrate the function and return a [`GaussKronrod`] integration result.
    ///
    /// # Errors
    /// Integration can fail if user suplied tolerance cannot be achieved within the maximum number
    /// of iterations.
    pub fn integrate(&self) -> Result<Adaptive, Error> {
        let integral = GaussKronrodBasic::new(self.lower, self.upper, self.rule, self.function)
            .integrate_internal();

        let mut iterations: usize = 1;

        let tolerance = self.error_bound.tolerance(&integral.result());
        let roundoff = integral.roundoff();

        if integral.error() <= roundoff && integral.error() > tolerance {
            let partial_result = Adaptive::from_basic(&integral, iterations);
            let kind = Kind::FailedToReachToleranceRoundoff;
            return Err(Error::new(kind, partial_result));
        } else if (integral.error() <= tolerance
            && integral.error().to_bits() != integral.result_asc().to_bits())
            || integral.error() == 0.0
        {
            let output = Adaptive::from_basic(&integral, iterations);
            return Ok(output);
        } else if self.max_iterations == 1 {
            let partial_result = Adaptive::from_basic(&integral, iterations);
            let kind = Kind::MaximumIterationsReached;
            return Err(Error::new(kind, partial_result));
        }

        let mut area = integral.result();
        let mut error = integral.error();

        let mut results = BinaryHeap::with_capacity(2 * self.max_iterations + 1);
        results.push(integral);

        let mut roundoff_type1 = 0usize;
        let mut roundoff_type2 = 0usize;

        while iterations < self.max_iterations {
            let Some(previous) = results.pop() else {
                // TODO this should be an error? we should _always_ have something to pop
                break;
            };

            iterations += 1;

            let lower = previous.lower();
            let upper = previous.upper();
            let mid = (upper + lower) * 0.5;

            let lower_integral =
                GaussKronrodBasic::new(lower, mid, self.rule, self.function).integrate_internal();

            let upper_integral =
                GaussKronrodBasic::new(mid, upper, self.rule, self.function).integrate_internal();

            let iteration_area = lower_integral.result() + upper_integral.result();
            let iteration_error = lower_integral.error() + upper_integral.error();

            area += iteration_area - previous.result();
            error += iteration_error - previous.error();

            if lower_integral.result_asc().to_bits() != lower_integral.error().to_bits()
                && upper_integral.result_asc().to_bits() != upper_integral.error().to_bits()
            {
                let delta = previous.result() - iteration_area;

                if delta.abs() <= 1.0e-5 * iteration_area.abs()
                    && iteration_error >= 0.99 * previous.error()
                {
                    roundoff_type1 += 1;
                }
                if iterations >= 10 && iteration_error >= previous.error() {
                    roundoff_type2 += 1;
                }
            }

            let iteration_tolerance = self.error_bound.tolerance(&area);

            if error > iteration_tolerance {
                if roundoff_type1 >= 6 || roundoff_type2 >= 20 {
                    let partial_result = Adaptive::new(area, error, iterations);
                    let kind = Kind::FailedToReachToleranceRoundoff;
                    return Err(Error::new(kind, partial_result));
                }

                if subinterval_too_small(lower, mid, upper) {
                    let partial_result = Adaptive::new(area, error, iterations);
                    let kind = Kind::PossibleSingularity { lower, upper };
                    return Err(Error::new(kind, partial_result));
                }
            }

            results.push(lower_integral);
            results.push(upper_integral);

            if error < iteration_tolerance {
                break;
            }
        }

        let result = results.into_iter().fold(0.0f64, |a, v| a + v.result());
        let output = Adaptive::new(result, error, iterations);

        if iterations == self.max_iterations {
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
}

#[derive(Debug)]
pub struct Error {
    kind: Kind,
    integral: Adaptive,
}

impl Error {
    fn new(kind: Kind, integral: Adaptive) -> Self {
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

fn subinterval_too_small(lower: f64, midpoint: f64, upper: f64) -> bool {
    let eps = f64::EPSILON;
    let min = f64::MIN_POSITIVE;

    let tmp = (1.0 + 100.0 * eps) * (midpoint.abs() + 1000.0 * min);

    lower.abs() <= tmp && upper.abs() <= tmp
}
