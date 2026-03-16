use std::cmp::Ordering;

use crate::multi::Integrator;
use crate::multi::Rule;
use crate::Limits;
use crate::MultiDimensionalIntegrand;
use crate::ScalarF64;

/// The estimated integral value of a calculated region.
///
/// This is the output of the core [`Integrator`], used internally by the different numerical
/// integration routines in [`multi`]. It contains an estimate of the numerical integration
/// `result`, the estimated `error`, and then further information about the integration which is
/// used in the [`Adaptive`] routines for error estimation and analysis. A key feature of the
/// adaptive routines is the use of a [`BinaryHeap`] for storing the various [`Region`]s which have
/// been integrated. As such, the [`Region`] type implements the [`Ord`] and [`PartialOrd`] traits
/// and order against the `error` (primary) and `limits` (secondary) fields. Thus, when pushed on
/// to the [`BinaryHeap`], the [`Region`] with the largest `error` estimate _or_ if all have equal
/// `error` the region with the smallest [`Limits::length`] will be obtained first upon using
/// [`BinaryHeap::pop`]. This offloads the ordering of the [`Region`]s to the standard libary
/// implementation of the [`BinaryHeap`].
///
/// [`Integrator`]: crate::quadrature::Integrator
/// [`Adaptive`]: crate::quadrature::Adaptive
#[derive(Debug, Clone)]
pub(crate) struct Region<T: ScalarF64, const N: usize> {
    pub(crate) error: f64,
    pub(crate) result: T,
    pub(crate) limits: [Limits; N],
    pub(crate) bisection_axis: usize,
    pub(crate) volume: f64,
}

impl<T: ScalarF64, const N: usize> PartialEq for Region<T, N> {
    fn eq(&self, other: &Self) -> bool {
        (self.result == other.result)
            && (self.error == other.error)
            && (self.bisection_axis == other.bisection_axis)
            && (self.limits == other.limits)
            && (self.volume == other.volume)
    }
}

impl<T: ScalarF64, const N: usize> Eq for Region<T, N> {}

impl<T: ScalarF64, const N: usize> PartialOrd for Region<T, N> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<T: ScalarF64, const N: usize> Ord for Region<T, N> {
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

impl<T: ScalarF64, const N: usize> Region<T, N> {
    pub(crate) fn unevaluated() -> Self {
        let zero = T::zero();
        Self {
            error: 0.0,
            result: zero,
            limits: [Limits::new(0.0, 0.0); N],
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

    pub(crate) fn with_limits(mut self, limits: [Limits; N]) -> Self {
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

    pub(crate) fn limits(&self) -> &[Limits; N] {
        &self.limits
    }

    pub(crate) fn bisection_axis(&self) -> usize {
        self.bisection_axis
    }

    /// Primary comparison function for the [`Region`].
    ///
    /// Uses [`f64::total_cmp`] to determine the ordering of two [`Region`]s based on the `error`
    /// estimates.
    pub(crate) fn total_cmp_error(&self, other: &Self) -> Ordering {
        self.error.total_cmp(&other.error)
    }

    /// Secondary comparison function for the [`Region`].
    ///
    /// In the event that two regions have the same `error` estimate we use the _inverse_ length of
    /// the integration region to determine the ordering. The rational here is that if the `error`
    /// due to each region is the same then the smaller region is likely to be more problematic and
    /// should be tackled first.
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
        I: MultiDimensionalIntegrand<N, Scalar = T>,
        const J: usize,
        const K: usize,
    >(
        &self,
        function: &I,
        rule: &Rule<N, J, K>,
    ) -> [Region<T, N>; 2] {
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
