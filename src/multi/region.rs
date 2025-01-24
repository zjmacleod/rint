use num_traits::Zero;
use std::cmp::Ordering;

use crate::multi::Integrator;
use crate::multi::Rule;
use crate::Limits;
use crate::MultiDimensionalIntegrand;
use crate::ScalarF64;

#[derive(Debug)]
pub(crate) struct Region<T, const NDIM: usize> {
    pub(crate) error: f64,
    pub(crate) result: T,
    pub(crate) limits: [Limits; NDIM],
    pub(crate) bisection_axis: usize,
    pub(crate) evaluations: usize,
    pub(crate) volume: f64,
}

impl<T: ScalarF64, const NDIM: usize> PartialEq for Region<T, NDIM> {
    fn eq(&self, other: &Self) -> bool {
        (self.result == other.result)
            && (self.error == other.error)
            && (self.bisection_axis == other.bisection_axis)
            && (self.evaluations == other.evaluations)
            && (self.limits == other.limits)
            && (self.volume == other.volume)
    }
}

impl<T: ScalarF64, const NDIM: usize> Eq for Region<T, NDIM> {}

impl<T: ScalarF64, const NDIM: usize> PartialOrd for Region<T, NDIM> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<T: ScalarF64, const NDIM: usize> Ord for Region<T, NDIM> {
    fn cmp(&self, other: &Self) -> Ordering {
        let mut ordering = self.total_cmp_error(other);
        if let Ordering::Equal = ordering {
            ordering = self.total_cmp_interval_length(other);
        }
        if let Ordering::Equal = ordering {
            ordering = self.total_cmp_volume(other);
        }
        ordering
    }
}

impl<T: ScalarF64, const NDIM: usize> Region<T, NDIM> {
    pub(crate) fn unevaluated() -> Self {
        let zero = T::zero();
        Self {
            error: 0.0,
            result: zero,
            limits: [Limits::new(0.0, 0.0); NDIM],
            bisection_axis: 0,
            evaluations: 0,
            volume: 0.0,
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

    pub(crate) fn with_limits(mut self, limits: [Limits; NDIM]) -> Self {
        self.limits = limits;
        self
    }

    pub(crate) fn with_bisection_axis(mut self, bisection_axis: usize) -> Self {
        self.bisection_axis = bisection_axis;
        self
    }

    pub(crate) fn with_evaluations(mut self, evaluations: usize) -> Self {
        self.evaluations = evaluations;
        self
    }

    pub(crate) fn with_volume(mut self, volume: f64) -> Self {
        self.volume = volume;
        self
    }

    pub(crate) fn error(&self) -> f64 {
        self.error
    }

    pub(crate) fn result(&self) -> T {
        self.result
    }

    pub(crate) fn limits(&self) -> &[Limits; NDIM] {
        &self.limits
    }

    pub(crate) fn bisection_axis(&self) -> usize {
        self.bisection_axis
    }

    pub(crate) fn evaluations(&self) -> usize {
        self.evaluations
    }

    pub(crate) fn total_cmp_error(&self, other: &Self) -> Ordering {
        self.error.total_cmp(&other.error)
    }

    pub(crate) fn total_cmp_interval_length(&self, other: &Self) -> Ordering {
        let width = self.limits[self.bisection_axis].width().abs();
        let other_width = other.limits[other.bisection_axis].width().abs();
        let inverse_length = 1.0 / width;
        let other_inverse_length = 1.0 / other_width;
        inverse_length.total_cmp(&other_inverse_length)
    }

    pub(crate) fn total_cmp_volume(&self, other: &Self) -> Ordering {
        let inverse_volume = 1.0 / self.volume.abs();
        let other_inverse_volume = 1.0 / other.volume.abs();
        inverse_volume.total_cmp(&other_inverse_volume)
    }

    #[allow(clippy::needless_borrow)]
    pub(crate) fn bisect<I: MultiDimensionalIntegrand<NDIM>, const J: usize, const K: usize>(
        &self,
        function: &I,
        rule: &Rule<NDIM, J, K>,
    ) -> [Region<I::Scalar, NDIM>; 2] {
        let axis_to_bisect = self.bisection_axis;
        let previous_limits = self.limits();

        let [lower, upper] = previous_limits[axis_to_bisect].bisect();

        let mut lower_limits = *previous_limits;
        lower_limits[axis_to_bisect] = lower;

        let mut upper_limits = *previous_limits;
        upper_limits[axis_to_bisect] = upper;

        let lower_integral = Integrator::new(&function, &rule, lower_limits).integrate();
        let upper_integral = Integrator::new(&function, &rule, upper_limits).integrate();

        [lower_integral, upper_integral]
    }
}
