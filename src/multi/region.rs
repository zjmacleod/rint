use crate::multi::basic::MultiDimensionalBasic;
use crate::multi::{Kind, MultiDimensionalError};
use crate::Limits;
use crate::MultiDimensionalIntegrand;
use crate::ScalarF64;
use num_complex::ComplexFloat;
use num_traits::Zero;

pub(crate) struct Region<const NDIM: usize, T: ScalarF64> {
    pub(crate) error: f64,
    pub(crate) result: T,
    pub(crate) result_abs: f64,
    pub(crate) error_abs: f64,
    pub(crate) limits: [Limits; NDIM],
    pub(crate) largest_error_axis: usize,
    pub(crate) function_evaluations: usize,
}

impl<const NDIM: usize, T: ScalarF64> Region<NDIM, T> {
    pub(crate) fn unevaluated() -> Self {
        let zero = <T as Zero>::zero();
        Self {
            error: 0.0,
            result: zero,
            result_abs: 0.0,
            error_abs: 0.0,
            limits: [Limits::new(0.0, 0.0); NDIM],
            largest_error_axis: 0,
            function_evaluations: 0,
        }
    }

    pub(crate) fn with_error(mut self, error: f64) -> Self {
        self.error = error;
        self
    }

    pub(crate) fn with_result(mut self, result: T) -> Self {
        self.result = result;
        self
    }

    pub(crate) fn with_result_abs(mut self, abs: f64) -> Self {
        self.result_abs = abs;
        self
    }

    pub(crate) fn with_error_abs(mut self, error_abs: f64) -> Self {
        self.error_abs = error_abs;
        self
    }

    pub(crate) fn with_limits(mut self, limits: [Limits; NDIM]) -> Self {
        self.limits = limits;
        self
    }

    pub(crate) fn with_largest_error_axis(mut self, largest_error_axis: usize) -> Self {
        self.largest_error_axis = largest_error_axis;
        self
    }

    pub(crate) fn with_function_evaluations(mut self, function_evaluations: usize) -> Self {
        self.function_evaluations = function_evaluations;
        self
    }

    pub fn error(&self) -> f64 {
        self.error
    }

    pub fn result(&self) -> T {
        self.result
    }

    pub fn result_abs(&self) -> f64 {
        self.result_abs
    }

    pub fn error_abs(&self) -> f64 {
        self.error_abs
    }

    pub fn limits(&self) -> &[Limits; NDIM] {
        &self.limits
    }

    pub fn largest_error_axis(&self) -> usize {
        self.largest_error_axis
    }

    pub fn function_evaluations(&self) -> usize {
        self.function_evaluations
    }

    pub(crate) fn bisect<I: MultiDimensionalIntegrand<NDIM>>(
        &self,
        function: &I,
    ) -> [Region<NDIM, I::Scalar>; 2] {
        let axis_to_bisect = self.largest_error_axis;
        let previous_limits = self.limits();

        let [lower, upper] = previous_limits[axis_to_bisect].bisect();

        let mut lower_limits = *previous_limits;
        lower_limits[axis_to_bisect] = lower;

        let mut upper_limits = *previous_limits;
        upper_limits[axis_to_bisect] = upper;

        let lower_integral =
            MultiDimensionalBasic::new_unchecked(lower_limits, function).integrate_internal();
        let upper_integral =
            MultiDimensionalBasic::new_unchecked(upper_limits, function).integrate_internal();

        [lower_integral, upper_integral]
    }
}
