use std::collections::BinaryHeap;

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

            let (result, error) = workspace.improved_result_error(&previous, &lower, &upper);

            let iteration_tolerance = self.tolerance.tolerance(&result);

            workspace.push(lower);
            workspace.push(upper);

            if error <= iteration_tolerance {
                break;
            }

            //workspace.check_roundoff()?;
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

    fn initialise_workspace(&self, initial: Region<I::Scalar, NDIM>) -> Workspace<I::Scalar, NDIM> {
        let iteration = 1;
        let result = initial.result();
        let error = initial.error();
        let limits = *initial.limits();
        let evaluations_per_integration = self.rule.evaluations();
        let bisection_axis = initial.bisection_axis();

        let mut heap = BinaryHeap::with_capacity(2 * self.max_iterations + 1);
        heap.push(initial);

        Workspace {
            heap,
            iteration,
            result,
            error,
            limits,
            bisection_axis,
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
    bisection_axis: usize,
    evaluations_per_integration: usize,
}

impl<T: ScalarF64, const NDIM: usize> Workspace<T, NDIM> {
    fn retrieve_largest_error(&mut self) -> Result<Region<T, NDIM>, IntegrationError<T>> {
        self.iteration += 1;
        if let Some(previous) = self.pop() {
            self.limits = *previous.limits();
            self.bisection_axis = previous.bisection_axis();
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

    fn check_singularity(&self) -> Result<(), IntegrationError<T>> {
        let limits = self.limits[self.bisection_axis];
        if limits.subinterval_too_small() {
            let output = self.integral_estimate();
            let kind = IntegrationErrorKind::BadIntegrandBehaviour(limits);
            Err(IntegrationError::new(output, kind))
        } else {
            Ok(())
        }
    }

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
mod tests_relative {
    use super::*;
    use crate::multi::{Rule07, Rule09, Rule11, Rule13};

    #[test]
    fn compare_adaptive_7point_with_dcuhre_output_ndim_2() {
        const NDIM: usize = 2;
        const TOL: f64 = 1e-7;

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
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.79659959933664204;
            let dcuhre_error = 5.2361990216203975E-008;
            let dcuhre_iter = 8;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-13);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-7);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let exponent = -(x * x + y * y);
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.55774628528556447;
            let dcuhre_error = 3.6639853981745754E-008;
            let dcuhre_iter = 15;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-7);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let exponent = -((x * x).cos().powi(2) * (y * y).cos().powi(2));
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.51478982537195195;
            let dcuhre_error = 5.1101592647406511E-008;
            let dcuhre_iter = 32;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-7);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    f64::cos((1.0 + x * x).ln() * (1.0 + y * y).ln())
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.99331725144532956;
            let dcuhre_error = 8.1154853602071175E-008;
            let dcuhre_iter = 13;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-7);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    (x * x / (2.0 - x.cos())) + (y * y / (2.0 - y.cos()))
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limits = [Limits::new(-2.0, 3.0), Limits::new(1.0, 10.0)];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 915.16141422826388;
            let dcuhre_error = 9.0498160158280369E-005;
            let dcuhre_iter = 122;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-14);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-6);
        }
    }

    #[test]
    fn compare_adaptive_7point_with_dcuhre_output_ndim_3() {
        const NDIM: usize = 3;
        const TOL: f64 = 1e-7;

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent = -x * y * z;
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.89121279725479763;
            let dcuhre_error = 8.4688451939699204E-008;
            let dcuhre_iter = 29;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-8);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent = -(x * x + y * y + z * z);
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.41653838589891773;
            let dcuhre_error = 4.1309505449997359E-008;
            let dcuhre_iter = 50;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-7);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent =
                        -((x * x).cos().powi(2) * (y * y).cos().powi(2) * (z * z).cos().powi(2));
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.57854055859751030;
            let dcuhre_error = 5.7504381834204602E-008;
            let dcuhre_iter = 161;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-14);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-6);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    f64::cos((1.0 + x * x).ln() * (1.0 + y * y).ln() * (1.0 + z * z).ln())
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.99922121344591008;
            let dcuhre_error = 9.8730425677660813E-008;
            let dcuhre_iter = 19;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-8);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    (x * x / (2.0 - x.cos()))
                        + (y * y / (2.0 - y.cos()))
                        + (z * z / (2.0 - z.cos()))
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limits = [
                Limits::new(-2.0, 3.0),
                Limits::new(1.0, 10.0),
                Limits::new(0.0, -1.0),
            ];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = -926.96331826144592;
            let dcuhre_error = 9.1088004807296453E-005;
            let dcuhre_iter = 323;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-14);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-6);
        }
    }

    #[test]
    fn compare_adaptive_9point_with_dcuhre_output_ndim_3() {
        const NDIM: usize = 3;
        const TOL: f64 = 1e-7;

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent = -x * y * z;
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule09::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.89121279904294926;
            let dcuhre_error = 2.9095568211933360E-008;
            let dcuhre_iter = 11;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-9);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent = -(x * x + y * y + z * z);
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule09::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];
            let tol = Tolerance::Relative(1e-7);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.41653838877485561;
            let dcuhre_error = 3.9221810923491651E-008;
            let dcuhre_iter = 3;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-8);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent =
                        -((x * x).cos().powi(2) * (y * y).cos().powi(2) * (z * z).cos().powi(2));
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule09::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.57854055859929987;
            let dcuhre_error = 4.8429415555578287E-008;
            let dcuhre_iter = 30;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-8);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    f64::cos((1.0 + x * x).ln() * (1.0 + y * y).ln() * (1.0 + z * z).ln())
                }
            }

            let function = Function;
            let rule = Rule09::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.99922122294090299;
            let dcuhre_error = 8.6788242268484642E-008;
            let dcuhre_iter = 4;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-9);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    (x * x / (2.0 - x.cos()))
                        + (y * y / (2.0 - y.cos()))
                        + (z * z / (2.0 - z.cos()))
                }
            }

            let function = Function;
            let rule = Rule09::generate().unwrap();
            let limits = [
                Limits::new(-2.0, 3.0),
                Limits::new(1.0, 10.0),
                Limits::new(0.0, -1.0),
            ];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = -926.96332279388207;
            let dcuhre_error = 8.9393505865191549E-005;
            let dcuhre_iter = 49;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-9);
        }
    }

    #[test]
    fn compare_adaptive_11point_with_dcuhre_output_ndim_3() {
        const NDIM: usize = 3;
        const TOL: f64 = 1e-7;

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent = -x * y * z;
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule11::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.89121279809085707;
            let dcuhre_error = 2.7684230235366546E-009;
            let dcuhre_iter = 2;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-15);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent = -(x * x + y * y + z * z);
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule11::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.41653838589354836;
            let dcuhre_error = 1.5049516251831965E-009;
            let dcuhre_iter = 2;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-13);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent =
                        -((x * x).cos().powi(2) * (y * y).cos().powi(2) * (z * z).cos().powi(2));
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule11::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.57854056257168118;
            let dcuhre_error = 3.3561510766170479E-008;
            let dcuhre_iter = 11;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-11);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    f64::cos((1.0 + x * x).ln() * (1.0 + y * y).ln() * (1.0 + z * z).ln())
                }
            }

            let function = Function;
            let rule = Rule11::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.99922121087502347;
            let dcuhre_error = 3.3840200784563892E-008;
            let dcuhre_iter = 2;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-15);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    (x * x / (2.0 - x.cos()))
                        + (y * y / (2.0 - y.cos()))
                        + (z * z / (2.0 - z.cos()))
                }
            }

            let function = Function;
            let rule = Rule11::generate();
            let limits = [
                Limits::new(-2.0, 3.0),
                Limits::new(1.0, 10.0),
                Limits::new(0.0, -1.0),
            ];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = -926.96326270996929;
            let dcuhre_error = 9.1854004280122431E-005;
            let dcuhre_iter = 42;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-9);
        }
    }

    #[test]
    fn compare_adaptive_13point_with_dcuhre_output_ndim_2() {
        const NDIM: usize = 2;
        const TOL: f64 = 1e-7;

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
            let rule = Rule13::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.79659959929705360;
            let dcuhre_error = 4.2849519163043052E-010;
            let dcuhre_iter = 2;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-15);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let exponent = -(x * x + y * y);
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule13::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.55774628535112658;
            let dcuhre_error = 2.1148653359684426E-011;
            let dcuhre_iter = 2;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-15);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let exponent = -((x * x).cos().powi(2) * (y * y).cos().powi(2));
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule13::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.51478982547960483;
            let dcuhre_error = 7.6107864567909130E-009;
            let dcuhre_iter = 7;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-11);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    f64::cos((1.0 + x * x).ln() * (1.0 + y * y).ln())
                }
            }

            let function = Function;
            let rule = Rule13::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.99331725120714986;
            let dcuhre_error = 1.0234173072687019E-010;
            let dcuhre_iter = 2;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-15);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    (x * x / (2.0 - x.cos())) + (y * y / (2.0 - y.cos()))
                }
            }

            let function = Function;
            let rule = Rule13::generate();
            let limits = [Limits::new(-2.0, 3.0), Limits::new(1.0, 10.0)];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(function, rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 915.16141806432074;
            let dcuhre_error = 8.3877707439014823E-005;
            let dcuhre_iter = 37;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-10);
        }
    }
}
