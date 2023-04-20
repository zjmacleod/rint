use std::collections::BinaryHeap;

use crate::integration::basic::BasicInternal;
use crate::integration::{ErrorBound, GaussKronrodBasic, IntegrationError};
use crate::rule::Rule;
use crate::Integrand;

/// The value of a function evaluated with Gauss-Kronrod integration and associated error
/// estimation.
pub struct Adaptive {
    result: f64,
    error: f64,
    iterations: usize,
}

impl Adaptive {
    /// Return the numerically approximated value of the integral.
    #[must_use]
    pub fn result(&self) -> f64 {
        self.result
    }

    /// Return the numerically approximated value of the absolute value of the integral.
    #[must_use]
    pub fn result_abs(&self) -> f64 {
        self.result_abs
    }

    #[must_use]
    pub fn result_asc(&self) -> f64 {
        self.result_asc
    }

    /// Return the numerically approximated error.
    #[must_use]
    pub fn error(&self) -> f64 {
        self.error
    }

    pub(crate) fn from_basic(basic: &BasicInternal, iterations: usize) -> Self {
        let result = basic.result();
        let result_abs = basic.result_abs();
        let result_asc = basic.result_asc();
        let error = basic.error();
        Self {
            result,
            result_abs,
            result_asc,
            error,
            iterations,
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
    ) -> Result<Self, IntegrationError<Adaptive>> {
        match error_bound {
            ErrorBound::Absolute(v) => {
                if v <= 0.0 {
                    return Err(IntegrationError::RelativeBoundNegativeOrZero(v));
                }
            }
            ErrorBound::Relative(v) => {
                if v < 50.0 * f64::EPSILON {
                    return Err(IntegrationError::AbsoluteBoundTooSmall(v));
                }
            }
            ErrorBound::Either { absolute, relative } => {
                if absolute <= 0.0 && relative < 50.0 * f64::EPSILON {
                    return Err(IntegrationError::InvalidTolerance);
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
    pub fn integrate(&self) -> Result<Adaptive, IntegrationError<Adaptive>> {
        let initial_integration = GaussKronrodBasic::new(
            self.lower,
            self.upper,
            self.rule,
            self.function,
        )
        .integrate_internal();

        let tolerance = self.error_bound.tolerance(&initial_integration.result());
        let roundoff = initial_integration.roundoff();

        if initial_integration.error() <= roundoff
            && initial_integration.error() > tolerance
        {
            let output = Adaptive::from_basic(&initial_integration, 1);
            return Err(IntegrationError::FailedToReachToleranceRoundoff(output));
        } else if (initial_integration.error() <= tolerance
            && initial_integration.error().to_bits()
                != initial_integration.result_asc().to_bits())
            || initial_integration.error() == 0.0
        {
            let output = Adaptive::from_basic(&initial_integration, 1);
            return Ok(output);
        } else if self.max_iterations == 1 {
            let output = Adaptive::from_basic(&initial_integration, 1);
            return Err(IntegrationError::MaximumSubintervalsReached(output));
        }

        let mut iterations: usize = 1;

        let mut area = initial_integration.result();
        let mut error = initial_integration.error();

        let mut results = BinaryHeap::with_capacity(2 * self.max_iterations + 1);
        results.push(initial_integration);
        let mut roundoff_type1 = 0usize;
        let mut roundoff_type2 = 0usize;

        while iterations <= self.max_iterations {
            let Some(previous) = results.pop() else {
                break;
            };

            let midpoint = (previous.upper() + previous.lower()) * 0.5;

            let lower_integration = GaussKronrodBasic::new(
                previous.lower(),
                midpoint,
                self.rule,
                self.function,
            )
            .integrate_internal();

            let upper_integration = GaussKronrodBasic::new(
                midpoint,
                previous.upper(),
                self.rule,
                self.function,
            )
            .integrate_internal();

            let iteration_area =
                lower_integration.result() + upper_integration.result();
            let iteration_error =
                lower_integration.error() + upper_integration.error();

            area += iteration_area - previous.result();
            error += iteration_error - previous.error();

            if lower_integration.result_asc().to_bits()
                != lower_integration.error().to_bits()
                && upper_integration.result_asc().to_bits()
                    != upper_integration.error().to_bits()
            {
                let delta = previous.result() - iteration_area;

                if delta.abs() <= 1.0e-5 * iteration_area.abs()
                    && iteration_error >= 0.99 * previous.error()
                {
                    // TODO
                    roundoff_type1 += 1;
                }
                if iterations >= 10 && iteration_error >= previous.error() {
                    // TODO
                    roundoff_type2 += 1;
                }
            }

            let iteration_tolerance = self.error_bound.tolerance(&area);

            if error > iteration_tolerance {
                if roundoff_type1 >= 6 || roundoff_type2 >= 20 {
                    return Err(IntegrationError::Roundoff);
                }

                // TODO singularity check need subinterval_too_small
                if false {
                    return Err(IntegrationError::PossibleSingularity);
                }
            }

            results.push(lower_integration);
            results.push(upper_integration);

            iterations += 1;

            if error < iteration_tolerance {
                break;
            }
        }

        let result = results.into_iter().fold(0.0f64, |a, v| a + v.result());

        Ok(Adaptive {
            result,
            error,
            iterations,
        })
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
