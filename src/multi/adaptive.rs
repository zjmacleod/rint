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

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Adaptive<'a, I, R, const N: usize> {
    function: &'a I,
    rule: &'a R,
    limits: [Limits; N],
    tolerance: Tolerance,
    max_iterations: usize,
}

impl<'a, I, const N: usize, const FINAL: usize, const TOTAL: usize>
    Adaptive<'a, I, Rule<N, FINAL, TOTAL>, N>
where
    I: MultiDimensionalIntegrand<N>,
{
    /// Generate a new adaptive multi-dimensional integrator.
    ///
    /// # Errors
    ///
    /// Returns an [`InitialisationError`] if:
    /// - The requested dimensionality `N` is invalid. `N` must match the dimensionality of
    /// the selected fully-symmetric multi-dimensional integration [`Rule`] and (in general)
    /// satisfy `2 <= N <= 15`.
    /// - The requested [`Tolerance`] is invalid. The tolerance muse satisfy the following
    /// constraints:
    ///     - `Tolerance::Absolute(v)` where `v > 0.0`
    ///     - `Tolerance::Relative(v)` where `v > 50.0 * f64::EPSILON`
    ///     - `Tolerance::Either { absolute, relative }` where `absolute > 0.0 and relative > 50.0 * f64::EPSILON`
    pub fn new(
        function: &'a I,
        rule: &'a Rule<N, FINAL, TOTAL>,
        limits: [Limits; N],
        tolerance: Tolerance,
        max_iterations: usize,
    ) -> Result<Self, InitialisationError> {
        if N < 2 || N > 15 {
            return Err(InitialisationError::new(
                InitialisationErrorKind::InvalidDimension(N),
            ));
        }

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
        let initial = Integrator::new(&self.function, self.rule, self.limits).integrate();

        // TODO check initial integration?
        //
        // if let Some(output) = self.check_initial_integration(&initial)? {
        //      return Ok(output);
        // }

        let mut workspace = self.initialise_workspace(initial);

        while workspace.iteration < self.max_iterations {
            let previous = workspace.retrieve_largest_error()?;

            let [lower, upper] = previous.bisect(&self.function, self.rule);

            let (result, error) = workspace.improved_result_error(&previous, &lower, &upper);

            let iteration_tolerance = self.tolerance.tolerance(&result);

            workspace.push(lower);
            workspace.push(upper);

            if error <= iteration_tolerance {
                // TODO add individual checks on re and im error
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

    fn initialise_workspace(&self, initial: Region<I::Scalar, N>) -> Workspace<I::Scalar, N> {
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

struct Workspace<T: ScalarF64, const N: usize> {
    heap: BinaryHeap<Region<T, N>>,
    iteration: usize,
    result: T,
    error: f64,
    limits: [Limits; N],
    bisection_axis: usize,
    evaluations_per_integration: usize,
}

impl<T: ScalarF64, const N: usize> Workspace<T, N> {
    fn retrieve_largest_error(&mut self) -> Result<Region<T, N>, IntegrationError<T>> {
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

    fn pop(&mut self) -> Option<Region<T, N>> {
        self.heap.pop()
    }

    fn push(&mut self, integral: Region<T, N>) {
        self.heap.push(integral);
    }

    fn improved_result_error(
        &mut self,
        previous: &Region<T, N>,
        lower: &Region<T, N>,
        upper: &Region<T, N>,
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
    #[allow(clippy::too_many_lines)]
    fn compare_adaptive_7point_with_dcuhre_output_ndim_2() {
        const N: usize = 2;
        const TOL: f64 = 1e-7;

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let exponent = -x * y;
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; N];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.796_599_599_336_642_04;
            let dcuhre_error = 5.236_199_021_620_397_5E-008;
            let dcuhre_iter = 8;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-13);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-7);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let exponent = -(x * x + y * y);
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; N];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.557_746_285_285_564_47;
            let dcuhre_error = 3.663_985_398_174_575_4E-008;
            let dcuhre_iter = 15;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-7);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let exponent = -((x * x).cos().powi(2) * (y * y).cos().powi(2));
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; N];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.514_789_825_371_951_95;
            let dcuhre_error = 5.110_159_264_740_651_1E-008;
            let dcuhre_iter = 32;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-7);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    f64::cos((1.0 + x * x).ln() * (1.0 + y * y).ln())
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; N];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.993_317_251_445_329_56;
            let dcuhre_error = 8.115_485_360_207_117_5E-008;
            let dcuhre_iter = 13;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-7);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
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

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 915.161_414_228_263_88;
            let dcuhre_error = 9.049_816_015_828_036_9E-005;
            let dcuhre_iter = 122;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-14);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-6);
        }
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn compare_adaptive_7point_with_dcuhre_output_ndim_3() {
        const N: usize = 3;
        const TOL: f64 = 1e-7;

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
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
            let limits = [limit; N];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.891_212_797_254_797_63;
            let dcuhre_error = 8.468_845_193_969_920_4E-008;
            let dcuhre_iter = 29;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-8);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
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
            let limits = [limit; N];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.416_538_385_898_917_73;
            let dcuhre_error = 4.130_950_544_999_735_9E-008;
            let dcuhre_iter = 50;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-7);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
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
            let limits = [limit; N];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.578_540_558_597_510_30;
            let dcuhre_error = 5.750_438_183_420_460_2E-008;
            let dcuhre_iter = 161;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-14);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-6);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    f64::cos((1.0 + x * x).ln() * (1.0 + y * y).ln() * (1.0 + z * z).ln())
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; N];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.999_221_213_445_910_08;
            let dcuhre_error = 9.873_042_567_766_081_3E-008;
            let dcuhre_iter = 19;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-8);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
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

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = -926.963_318_261_445_92;
            let dcuhre_error = 9.108_800_480_729_645_3E-005;
            let dcuhre_iter = 323;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-14);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-6);
        }
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn compare_adaptive_9point_with_dcuhre_output_ndim_3() {
        const N: usize = 3;
        const TOL: f64 = 1e-7;

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
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
            let limits = [limit; N];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.891_212_799_042_949_26;
            let dcuhre_error = 2.909_556_821_193_336_0E-008;
            let dcuhre_iter = 11;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-9);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
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
            let limits = [limit; N];
            let tol = Tolerance::Relative(1e-7);
            let max_iterations = 1000;

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.416_538_388_774_855_61;
            let dcuhre_error = 3.922_181_092_349_165_1E-008;
            let dcuhre_iter = 3;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-8);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
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
            let limits = [limit; N];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.578_540_558_599_299_87;
            let dcuhre_error = 4.842_941_555_557_828_7E-008;
            let dcuhre_iter = 30;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-8);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    f64::cos((1.0 + x * x).ln() * (1.0 + y * y).ln() * (1.0 + z * z).ln())
                }
            }

            let function = Function;
            let rule = Rule09::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; N];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.999_221_222_940_902_99;
            let dcuhre_error = 8.678_824_226_848_464_2E-008;
            let dcuhre_iter = 4;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-9);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
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

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = -926.963_322_793_882_07;
            let dcuhre_error = 8.939_350_586_519_154_9E-005;
            let dcuhre_iter = 49;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-9);
        }
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn compare_adaptive_11point_with_dcuhre_output_ndim_3() {
        const N: usize = 3;
        const TOL: f64 = 1e-7;

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
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
            let limits = [limit; N];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.891_212_798_090_857_07;
            let dcuhre_error = 2.768_423_023_536_654_6E-009;
            let dcuhre_iter = 2;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-15);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
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
            let limits = [limit; N];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.416_538_385_893_548_36;
            let dcuhre_error = 1.504_951_625_183_196_5E-009;
            let dcuhre_iter = 2;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-13);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
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
            let limits = [limit; N];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.578_540_562_571_681_18;
            let dcuhre_error = 3.356_151_076_617_047_9E-008;
            let dcuhre_iter = 11;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-11);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    f64::cos((1.0 + x * x).ln() * (1.0 + y * y).ln() * (1.0 + z * z).ln())
                }
            }

            let function = Function;
            let rule = Rule11::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; N];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.999_221_210_875_023_47;
            let dcuhre_error = 3.384_020_078_456_389_2E-008;
            let dcuhre_iter = 2;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-15);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
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

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = -926.963_262_709_969_29;
            let dcuhre_error = 9.185_400_428_012_243_1E-005;
            let dcuhre_iter = 42;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-9);
        }
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn compare_adaptive_13point_with_dcuhre_output_ndim_2() {
        const N: usize = 2;
        const TOL: f64 = 1e-7;

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let exponent = -x * y;
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule13::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; N];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.796_599_599_297_053_60;
            let dcuhre_error = 4.284_951_916_304_305_2E-010;
            let dcuhre_iter = 2;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-15);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let exponent = -(x * x + y * y);
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule13::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; N];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.557_746_285_351_126_58;
            let dcuhre_error = 2.114_865_335_968_442_6E-011;
            let dcuhre_iter = 2;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-15);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let exponent = -((x * x).cos().powi(2) * (y * y).cos().powi(2));
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule13::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; N];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.514_789_825_479_604_83;
            let dcuhre_error = 7.610_786_456_790_913_0E-009;
            let dcuhre_iter = 7;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-11);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    f64::cos((1.0 + x * x).ln() * (1.0 + y * y).ln())
                }
            }

            let function = Function;
            let rule = Rule13::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; N];
            let tol = Tolerance::Relative(TOL);
            let max_iterations = 1000;

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 0.993_317_251_207_149_86;
            let dcuhre_error = 1.023_417_307_268_701_9E-010;
            let dcuhre_iter = 2;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-15);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<N> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
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

            let integral = Adaptive::new(&function, &rule, limits, tol, max_iterations).unwrap();

            let integral_result = integral.integrate().unwrap();

            let result = integral_result.result();
            let error = integral_result.error();
            let iter = integral_result.iterations();
            let dcuhre_result = 915.161_418_064_320_74;
            let dcuhre_error = 8.387_770_743_901_482_3E-005;
            let dcuhre_iter = 37;

            assert_eq!(iter, dcuhre_iter);
            assert!(error < TOL * result.abs());
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-10);
        }
    }
}
