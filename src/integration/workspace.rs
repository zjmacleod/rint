use std::collections::binary_heap::{BinaryHeap, IntoIter};

use crate::integration::adaptive::{Adaptive, Error, Kind};
use crate::integration::basic::BasicInternal;
use crate::integration::subinterval_too_small;

pub(crate) enum Routine {
    Adaptive(AdaptiveData),
    Singularity(SingularityData),
}

pub(crate) struct AdaptiveData {
    roundoff_type1: usize,
    roundoff_type2: usize,
}

impl AdaptiveData {
    pub(crate) fn new() -> AdaptiveData {
        Self {
            roundoff_type1: 0,
            roundoff_type2: 0,
        }
    }
}

pub(crate) struct SingularityData {
    roundoff_type1: usize,
    roundoff_type2: usize,
    roundoff_type3: usize,
    extrapolate: bool,
}

impl SingularityData {
    pub(crate) fn new() -> SingularityData {
        Self {
            roundoff_type1: 0,
            roundoff_type2: 0,
            roundoff_type3: 0,
            extrapolate: false,
        }
    }
}

pub(crate) struct Workspace {
    heap: BinaryHeap<BasicInternal>,
    iteration: usize,
    routine: Routine,
    result: f64,
    error: f64,
    lower_limit: f64,
    upper_limit: f64,
}

impl Workspace {
    fn new(max_iterations: usize, initial: BasicInternal, routine: Routine) -> Self {
        let mut heap = BinaryHeap::with_capacity(2 * max_iterations + 1);

        let result = initial.result();
        let error = initial.error();
        let lower_limit = initial.lower();
        let upper_limit = initial.upper();

        heap.push(initial);

        Self {
            heap,
            iteration: 1,
            routine,
            result,
            error,
            lower_limit,
            upper_limit,
        }
    }

    pub(crate) fn adaptive(max_iterations: usize, initial: BasicInternal) -> Self {
        let routine = Routine::Adaptive(AdaptiveData::new());
        Self::new(max_iterations, initial, routine)
    }

    pub(crate) fn singularity(max_iterations: usize, initial: BasicInternal) -> Self {
        let routine = Routine::Singularity(SingularityData::new());
        Self::new(max_iterations, initial, routine)
    }

    #[inline]
    pub(crate) fn iteration(&self) -> usize {
        self.iteration
    }

    pub(crate) fn retrieve_largest_error(&mut self) -> Result<BasicInternal, Error> {
        self.iteration += 1;
        if let Some(previous) = self.pop() {
            Ok(previous)
        } else {
            let output = Adaptive::empty();
            let kind = Kind::UninitialisedWorkspace;
            Err(Error::new(kind, output))
        }
    }

    pub(crate) fn update_limits(&mut self, lower: f64, upper: f64) {
        self.lower_limit = lower;
        self.upper_limit = upper;
    }

    pub(crate) fn update_result(&mut self, diff: f64) -> f64 {
        self.result += diff;
        self.result
    }

    pub(crate) fn update_error(&mut self, diff: f64) -> f64 {
        self.error += diff;
        self.error
    }

    pub(crate) fn iteration_result_error(
        &mut self,
        previous: &BasicInternal,
        lower: &BasicInternal,
        upper: &BasicInternal,
    ) -> (f64, f64) {
        let prev_result = previous.result();
        let prev_error = previous.error();
        let new_result = lower.result() + upper.result();
        let new_error = lower.error() + upper.error();

        match self.routine {
            Routine::Adaptive(ref mut state) => {
                if lower.result_asc().to_bits() != lower.error().to_bits()
                    && upper.result_asc().to_bits() != upper.error().to_bits()
                {
                    let delta = (prev_result - new_result).abs();

                    if delta <= 1e-5 * new_result.abs() && new_error >= 0.99 * prev_error {
                        state.roundoff_type1 += 1;
                    }
                    if self.iteration >= 10 && new_error >= prev_error {
                        state.roundoff_type2 += 1;
                    }
                }
            }
            Routine::Singularity(ref mut state) => {
                if lower.result_asc().to_bits() != lower.error().to_bits()
                    && upper.result_asc().to_bits() != upper.error().to_bits()
                {
                    let delta = (prev_result - new_result).abs();

                    if delta <= 1e-5 * new_result.abs() && new_error >= 0.99 * prev_error {
                        if state.extrapolate {
                            state.roundoff_type2 += 1;
                        } else {
                            state.roundoff_type1 += 1;
                        }
                    }
                    if self.iteration >= 10 && new_error >= prev_error {
                        state.roundoff_type3 += 1;
                    }
                }
            }
        }

        let result = self.update_result(new_result - prev_result);
        let error = self.update_error(new_error - prev_error);

        (result, error)
    }

    pub(crate) fn check_roundoff(&self) -> Result<(), Error> {
        match self.routine {
            Routine::Adaptive(ref state) => {
                if state.roundoff_type1 >= 6 || state.roundoff_type2 >= 20 {
                    let output = self.integral_estimate();
                    let kind = Kind::FailedToReachToleranceRoundoff;
                    return Err(Error::new(kind, output));
                }
            }
            Routine::Singularity(ref state) => {
                if state.roundoff_type1 + state.roundoff_type2 >= 10 || state.roundoff_type3 >= 20 {
                    let output = self.integral_estimate();
                    let kind = Kind::FailedToReachToleranceRoundoff;
                    return Err(Error::new(kind, output));
                }

                // TODO
                // if state.roundoff_type2 >= 5 {
                //     error_type2 = 1; // only used to inc
                // }
            }
        }
        Ok(())
    }

    pub(crate) fn check_singularity(&self) -> Result<(), Error> {
        let lower_limit = self.lower_limit;
        let upper_limit = self.upper_limit;
        let midpoint = (lower_limit + upper_limit) * 0.5;
        if subinterval_too_small(lower_limit, midpoint, upper_limit) {
            let output = self.integral_estimate();
            let kind = Kind::PossibleSingularity {
                lower: lower_limit,
                upper: upper_limit,
            };
            Err(Error::new(kind, output))
        } else {
            Ok(())
        }
    }

    pub(crate) fn sum_results(&self) -> f64 {
        self.iter().fold(0.0f64, |a, v| a + v.result())
    }

    pub(crate) fn integral_estimate(&self) -> Adaptive {
        let error = self.error;
        let iterations = self.iteration;
        let result = self.sum_results();
        Adaptive::new(result, error, iterations)
    }
}

impl std::ops::Deref for Workspace {
    type Target = BinaryHeap<BasicInternal>;
    fn deref(&self) -> &Self::Target {
        &self.heap
    }
}

impl std::ops::DerefMut for Workspace {
    fn deref_mut(&mut self) -> &mut BinaryHeap<BasicInternal> {
        &mut self.heap
    }
}

impl IntoIterator for Workspace {
    type Item = BasicInternal;
    type IntoIter = IntoIter<BasicInternal>;

    fn into_iter(self) -> Self::IntoIter {
        self.heap.into_iter()
    }
}
