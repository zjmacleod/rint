use std::cmp::Ordering;

use crate::multi::Integrator;
use crate::multi::Rule;
use crate::Limits;
use crate::MultiDimensionalIntegrand;
use crate::ScalarF64;

#[derive(Debug)]
pub(crate) struct Region<T: ScalarF64, const NDIM: usize> {
    pub(crate) error: f64,
    pub(crate) result: T,
    pub(crate) limits: [Limits; NDIM],
    pub(crate) bisection_axis: usize,
    pub(crate) volume: f64,
}

impl<T: ScalarF64, const NDIM: usize> PartialEq for Region<T, NDIM> {
    fn eq(&self, other: &Self) -> bool {
        (self.result == other.result)
            && (self.error == other.error)
            && (self.bisection_axis == other.bisection_axis)
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
            if let Ordering::Equal = ordering {
                ordering = self.total_cmp_volume(other);
            }
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
    pub(crate) fn bisect<
        I: MultiDimensionalIntegrand<NDIM, Scalar = T>,
        const J: usize,
        const K: usize,
    >(
        &self,
        function: &I,
        rule: &Rule<NDIM, J, K>,
    ) -> [Region<T, NDIM>; 2] {
        let axis_to_bisect = self.bisection_axis;
        let previous_limits = self.limits();
        let previous_result = self.result();

        let [lower, upper] = previous_limits[axis_to_bisect].bisect();

        let mut upper_limits = *previous_limits;
        upper_limits[axis_to_bisect] = upper;

        let mut lower_limits = *previous_limits;
        lower_limits[axis_to_bisect] = lower;

        let mut upper_integral = Integrator::new(&function, &rule, upper_limits).integrate();
        let mut lower_integral = Integrator::new(&function, &rule, lower_limits).integrate();

        // adjust error estimates according to Bernsten & Espelid 1991
        {
            let new_result = upper_integral.result() + lower_integral.result();

            let mut upper_error = upper_integral.error();
            let mut lower_error = lower_integral.error();

            let est1 = (previous_result - new_result).abs();
            let est2 = lower_error + upper_error;

            let error_coeff = rule.adaptive_error_coeff();

            if est2 > 0.0 {
                lower_error *= 1.0 + error_coeff.c5() * est1 / est2;
                upper_error *= 1.0 + error_coeff.c5() * est1 / est2;
            }

            upper_error += error_coeff.c6() * est1;
            lower_error += error_coeff.c6() * est1;

            upper_integral.error = upper_error;
            lower_integral.error = lower_error;
        }

        [lower_integral, upper_integral]
    }
}
