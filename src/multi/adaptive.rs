use num_complex::ComplexFloat;
use num_traits::Zero;

use std::collections::BinaryHeap;

use crate::multi::rule::AdaptiveErrorCoeff;
use crate::multi::Integrator;
use crate::multi::Region;
use crate::multi::Rule;
use crate::IntegralEstimate;
use crate::MultiDimensionalIntegrand;
use crate::ScalarF64;
use crate::{InitialisationError, InitialisationErrorKind};
use crate::{IntegrationError, IntegrationErrorKind};
use crate::{Limits, Tolerance};

pub struct Adaptive<I, R, const NDIM: usize> {
    function: I,
    rule: R,
    limits: [Limits; NDIM],
    tolerance: Tolerance,
    max_iterations: usize,
}

impl<I, const NDIM: usize, const FINAL: usize, const TOTAL: usize>
    Adaptive<I, Rule<NDIM, FINAL, TOTAL>, NDIM>
where
    I: MultiDimensionalIntegrand<NDIM>,
{
    /// Generate a new adaptive multi-dimensional integrator.
    ///
    /// # Errors
    ///
    /// Returns an [`InitialisationError`] if:
    /// - The requested dimensionality `NDIM` is invalid. `NDIM` must match the dimensionality of
    /// the selected fully-symmetric multi-dimensional integration [`Rule`] and (in general)
    /// satisfy `2 <= NDIM <= 15`.
    /// - The requested [`Tolerance`] is invalid. The tolerance muse satisfy the following
    /// constraints:
    ///     - `Tolerance::Absolute(v)` where `v > 0.0`
    ///     - `Tolerance::Relative(v)` where `v > 50.0 * f64::EPSILON`
    ///     - `Tolerance::Either { absolute, relative }` where `absolute > 0.0 and relative > 50.0 * f64::EPSILON`
    pub fn new(
        function: I,
        rule: Rule<NDIM, FINAL, TOTAL>,
        limits: [Limits; NDIM],
        tolerance: Tolerance,
        max_iterations: usize,
    ) -> Result<Self, InitialisationError> {
        if NDIM < 2 || NDIM > 15 {
            return Err(InitialisationError::new(
                InitialisationErrorKind::InvalidDimension(NDIM),
            ));
        };

        tolerance.check()?;

        Ok(Self {
            function,
            rule,
            limits,
            tolerance,
            max_iterations,
        })
    }

    /// Integrate the function using the adaptive routine and return a [`IntegralEstimate`].
    ///
    /// # Errors
    ///
    /// Returns an [`IntegrationError`] if:
    /// - TODO
    pub fn integrate(&self) -> Result<IntegralEstimate<I::Scalar>, IntegrationError<I::Scalar>> {
        let initial = Integrator::new(&self.function, &self.rule, self.limits).integrate();

        // TODO check initial integration?
        //
        // if let Some(output) = self.check_initial_integration(&initial)? {
        //      return Ok(output);
        // }

        let mut workspace = self.initialise_workspace(initial);

        while workspace.iteration < self.max_iterations {
            let previous = workspace.retrieve_largest_error()?;

            let [lower, upper] = previous.bisect(&self.function, &self.rule);

            let (result, error) = workspace.improved_result_error(
                &previous,
                &lower,
                &upper,
                &self.rule.adaptive_error_coeff(),
            );

            let iteration_tolerance = self.tolerance.tolerance(&result);

            workspace.push(lower);
            workspace.push(upper);

            if error <= iteration_tolerance {
                break;
            }

            //workspace.check_roundoff()?;
            //workspace.check_singularity()?;
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

    fn initialise_workspace(&self, initial: Region<I::Scalar, NDIM>) -> Workspace<I::Scalar, NDIM> {
        let iteration = 1;
        let result = initial.result();
        let error = initial.error();
        let limits = *initial.limits();
        let evaluations_per_integration = self.rule.evaluations();

        let mut heap = BinaryHeap::with_capacity(2 * self.max_iterations + 1);
        heap.push(initial);

        Workspace {
            heap,
            iteration,
            result,
            error,
            limits,
            evaluations_per_integration,
        }
    }
}

struct Workspace<T: ScalarF64, const NDIM: usize> {
    heap: BinaryHeap<Region<T, NDIM>>,
    iteration: usize,
    result: T,
    error: f64,
    limits: [Limits; NDIM],
    evaluations_per_integration: usize,
}

impl<T: ScalarF64, const NDIM: usize> Workspace<T, NDIM> {
    fn retrieve_largest_error(&mut self) -> Result<Region<T, NDIM>, IntegrationError<T>> {
        self.iteration += 1;
        if let Some(previous) = self.pop() {
            self.limits = *previous.limits();
            Ok(previous)
        } else {
            let kind = IntegrationErrorKind::UninitialisedWorkspace;
            let output = IntegralEstimate::new();
            Err(IntegrationError::new(output, kind))
        }
    }

    fn pop(&mut self) -> Option<Region<T, NDIM>> {
        self.heap.pop()
    }

    fn push(&mut self, integral: Region<T, NDIM>) {
        self.heap.push(integral);
    }

    fn improved_result_error(
        &mut self,
        previous: &Region<T, NDIM>,
        lower: &Region<T, NDIM>,
        upper: &Region<T, NDIM>,
        error_coef: &AdaptiveErrorCoeff,
    ) -> (T, f64) {
        let prev_result = previous.result();
        let prev_error = previous.error();
        let new_result = lower.result() + upper.result();
        let lower_error = lower.error();
        let upper_error = upper.error();

        // TODO check roundoff? see quadrature::Adaptive

        self.result += new_result - prev_result;
        self.error += (lower_error + upper_error) - prev_error;

        let result = self.result;
        let error = self.error;

        (result, error)
    }

    //fn check_roundoff(&self) -> Result<(), IntegrationError<T>> {
    //    if self.roundoff_count >= 6 || self.roundoff_on_high_iteration_count >= 20 {
    //        let output = self.integral_estimate();
    //        let kind = IntegrationErrorKind::RoundoffErrorDetected;
    //        return Err(IntegrationError::new(output, kind));
    //    }
    //    Ok(())
    //}

    //fn check_singularity(&self) -> Result<(), IntegrationError<T>> {
    //    let limits = self.limits;
    //    if subinterval_too_small(limits) {
    //        let output = self.integral_estimate();
    //        let kind = IntegrationErrorKind::BadIntegrandBehaviour(limits);
    //        Err(IntegrationError::new(output, kind))
    //    } else {
    //        Ok(())
    //    }
    //}

    fn sum_results(&self) -> T {
        self.heap.iter().fold(T::zero(), |a, v| a + v.result())
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

#[cfg(test)]
mod testswoo {
    use super::*;
    use crate::multi::{Rule07, Rule09, Rule11, Rule13};
    use crate::Error;

    #[test]
    fn compare_adaptive_7point_with_dcuhre_output_ndim_2() -> Result<(), Error<f64>> {
        const NDIM: usize = 2;

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let exponent = -x * y;
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];
            let tol = Tolerance::Relative(1e-7);
            let max_iterations = 100;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations)?;

            let integral_result = integral.integrate()?;

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.79659959933664204;
            let dcuhre_error = 5.2361990216203975E-008;
            let dcuhre_iter = 8;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-7);
        }
        Ok(())
    }
}
