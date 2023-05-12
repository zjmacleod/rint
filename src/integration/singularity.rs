use std::collections::binary_heap::{BinaryHeap, IntoIter};

use crate::integration::adaptive::{Error, Kind};
use crate::integration::basic::BasicInternal;
use crate::integration::{subinterval_too_small, Adaptive, ErrorBound, GaussKronrodBasic};
use crate::rule::Rule;
use crate::Integrand;

/// An integral with singularities to be evaluated with an adaptive Gauss-Kronrod quadrature.
///
/// The user constructs a `function` implementing [`Integrand`], provides `upper`
/// and `lower` integration limits, and provides an `error_bound`, which can be
/// [`ErrorBound::Absolute`] to work to a specified absolute error,
/// [`ErrorBound::Relative`] to work to a specified relative error,
/// or [`ErrorBound::Either`] to return a result as soon as _either_ the relative
/// or absolute error bound has been satisfied.
pub struct GaussKronrodAdaptiveSingularity<'a, I, R>
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

impl<'a, I, R> GaussKronrodAdaptiveSingularity<'a, I, R>
where
    I: Integrand,
    R: Rule,
{
    /// Create a new [`GaussKronrodAdaptiveSingularity`].
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

    fn roundoff(result_abs: f64) -> f64 {
        100.0 * f64::EPSILON * result_abs
    }

    fn check_initial_integration(
        &self,
        initial: &BasicInternal,
    ) -> Result<Option<Adaptive>, Error> {
        let tolerance = self.error_bound.tolerance(initial.result());
        let roundoff = Self::roundoff(initial.result_abs());

        if initial.error() <= roundoff && initial.error() > tolerance {
            let output = Adaptive::from_basic(initial, 1);
            let kind = Kind::FailedToReachToleranceRoundoff;

            Err(Error::new(kind, output))
        } else if (initial.error() <= tolerance
            && initial.error().to_bits() != initial.result_asc().to_bits())
            || initial.error() == 0.0
        {
            let output = Adaptive::from_basic(initial, 1);

            Ok(Some(output))
        } else if self.max_iterations == 1 {
            let output = Adaptive::from_basic(initial, 1);
            let kind = Kind::MaximumIterationsReached;

            Err(Error::new(kind, output))
        } else {
            Ok(None)
        }
    }

    /// # Errors
    // TODO add extrapolation
    // XXX goto compute_result = use workspace.integral_estimate()
    // XXX goto return_error = use table.integral_estimate()
    pub fn integrate(&self) -> Result<Adaptive, Error> {
        let initial = GaussKronrodBasic::new(self.lower, self.upper, self.rule, self.function)
            .integrate_internal();

        if let Some(output) = self.check_initial_integration(&initial)? {
            return Ok(output);
        }

        let mut workspace = Workspace::new(self.max_iterations, initial);

        while workspace.iteration < self.max_iterations {
            let previous = workspace.retrieve_largest_error()?;

            let half_current_interval = previous.abs_interval_length() / 2.0;
            workspace.update_smallest_interval(half_current_interval);

            let [lower, upper] = previous.bisect(self.function, self.rule);

            let previous_error = previous.error();
            let iteration_error = lower.error() + upper.error();

            let [result, error] = workspace.improved_result_error(&previous, &lower, &upper);

            let iteration_tolerance = self.error_bound.tolerance(result);

            workspace.push(lower);
            workspace.push(upper);

            if error <= iteration_tolerance {
                return workspace.go_to_115();
            }

            workspace.check_roundoff()?;
            workspace.check_singularity()?;

            if workspace.iteration >= self.max_iterations - 1 {
                // error type 1
                break;
            } else if workspace.iteration == 2 {
                // XXX quadpack multiplies this by 0.375
                workspace.smallest_interval = half_current_interval;
                workspace.error_over_large_intervals = workspace.error;
                workspace.ertest = iteration_tolerance;
                workspace.table.append_table(workspace.result);
                continue;
            }

            // probably don't need this. disallow_extrapolation _only_ set to true
            // if table.number == 1, which is only satisfied _before_ bisecting the
            // initial result outside the loop. By the time we reach this point,
            // table.number >= 2.
            //if disallow_extrapolation {
            //    continue;
            //}

            workspace.update_error_over_large_intervals(
                half_current_interval,
                previous_error,
                iteration_error,
            );

            // peek the next interval with the largest error and check if it
            // is the smallest interval.
            if !workspace.extrapolate {
                if let Some(next) = workspace.heap.peek() {
                    if next.abs_interval_length() > workspace.smallest_interval {
                        continue;
                    }
                };

                // the next interval to be bisected has the largest error and
                // smallest interval. Store it to reduce error in large intervals
                if let Some(smallest) = workspace.pop() {
                    workspace.store(smallest);
                };
                workspace.extrapolate = true;
            }

            if !workspace.error_type2 && workspace.error_over_large_intervals > workspace.ertest {
                // XXX not quite right. This should loop through to find the next
                // iteration with a larger interval. Store all the other?
                while let Some(next) = workspace.heap.peek() {
                    if next.abs_interval_length() > workspace.smallest_interval {
                        continue;
                    }
                    if let Some(next) = workspace.pop() {
                        workspace.store(next);
                    };
                }
            }

            let (ext_result, ext_error) = workspace.table.extrapolate(result);
        }

        let final_iteration = workspace.iteration;
        let output = workspace.integral_estimate();

        if final_iteration == self.max_iterations {
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

enum State {
    Initialised,
    Extrapolate,
    Bisect,
    ComputeResult,
    CheckError,
}

struct Workspace {
    heap: BinaryHeap<BasicInternal>,
    iteration: usize,
    roundoff_type1: usize,
    roundoff_type2: usize,
    roundoff_type3: usize,
    error_type2: bool,
    result: f64,
    error: f64,
    lower_limit: f64,
    upper_limit: f64,
    table: ExtrapolationTable,
    smallest_interval: f64,
    error_over_large_intervals: f64,
    ertest: f64,
    state: State,
    positive_integrand: bool,
    initial_abs: f64,
    correction: f64,
    extrapolate: bool,
    store: BinaryHeap<BasicInternal>,
    num_stored: usize,
}

impl Workspace {
    fn new(max_iterations: usize, initial: BasicInternal) -> Self {
        let mut heap = BinaryHeap::with_capacity(2 * max_iterations + 1);

        let result = initial.result();
        let error = initial.error();
        let lower_limit = initial.lower();
        let upper_limit = initial.upper();
        let table = ExtrapolationTable::initialise(&initial);
        let smallest_interval = initial.abs_interval_length();
        let error_over_large_intervals = error;
        let ertest = f64::MAX;
        let state = State::Initialised;
        let positive_integrand = initial.positivity();
        let initial_abs = initial.result_abs();
        let correction = 0.0;
        let store = BinaryHeap::with_capacity(2 * max_iterations + 1);

        heap.push(initial);

        Self {
            heap,
            iteration: 1,
            roundoff_type1: 0,
            roundoff_type2: 0,
            roundoff_type3: 0,
            error_type2: false,
            result,
            error,
            lower_limit,
            upper_limit,
            table,
            smallest_interval,
            error_over_large_intervals,
            ertest,
            state,
            positive_integrand,
            initial_abs,
            correction,
            extrapolate: false,
            store,
            num_stored: 0,
        }
    }

    fn retrieve_largest_error(&mut self) -> Result<BasicInternal, Error> {
        self.state = State::Bisect;
        self.iteration += 1;
        if let Some(previous) = self.pop() {
            self.lower_limit = previous.lower;
            self.upper_limit = previous.upper;

            Ok(previous)
        } else {
            let output = Adaptive::empty();
            let kind = Kind::UninitialisedWorkspace;
            Err(Error::new(kind, output))
        }
    }

    fn next_interval_length(&self) -> f64 {
        if let Some(next) = self.heap.peek() {
            next.abs_interval_length()
        } else {
            f64::MAX
        }
    }

    fn update_smallest_interval(&mut self, half_current_interval: f64) {
        if half_current_interval < self.smallest_interval {
            self.smallest_interval = half_current_interval;
        }
    }

    fn improved_result_error(
        &mut self,
        previous: &BasicInternal,
        lower: &BasicInternal,
        upper: &BasicInternal,
    ) -> [f64; 2] {
        let prev_result = previous.result();
        let prev_error = previous.error();
        let new_result = lower.result() + upper.result();
        let new_error = lower.error() + upper.error();

        if lower.result_asc().to_bits() != lower.error().to_bits()
            && upper.result_asc().to_bits() != upper.error().to_bits()
        {
            let delta = (prev_result - new_result).abs();

            if delta <= 1e-5 * new_result.abs() && new_error >= 0.99 * prev_error {
                if let State::Bisect = self.state {
                    self.roundoff_type1 += 1;
                } else if let State::Extrapolate = self.state {
                    self.roundoff_type2 += 1;
                }
                // in GSL this next conditional is _always_ run, in QUADPACK it is not
            } else if self.iteration > 10 && new_error > prev_error {
                self.roundoff_type3 += 1;
            }
        }

        self.result += new_result - prev_result;
        self.error += new_error - prev_error;

        [self.result(), self.error()]
    }

    // TODO Change these to go_to_100
    fn check_roundoff(&mut self) -> Result<(), Error> {
        if self.roundoff_type1 + self.roundoff_type2 >= 10 || self.roundoff_type3 >= 20 {
            let output = self.integral_estimate();
            let kind = Kind::FailedToReachToleranceRoundoff;
            return Err(Error::new(kind, output));
        } // goto 100

        if self.roundoff_type2 >= 5 {
            self.error_type2 = true;
        }
        Ok(())
    }

    // TODO Change these to go_to_100
    fn check_singularity(&self) -> Result<(), Error> {
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

    fn sum_results(&self) -> f64 {
        let store = self.store.iter().fold(0.0f64, |a, v| a + v.result());
        let heap = self.heap.iter().fold(0.0f64, |a, v| a + v.result());
        store + heap
    }

    fn store(&mut self, integral: BasicInternal) {
        self.store.push(integral);
        self.num_stored += 1;
    }

    fn pop(&mut self) -> Option<BasicInternal> {
        self.heap.pop()
    }

    fn push(&mut self, integral: BasicInternal) {
        self.heap.push(integral);
    }

    fn integral_estimate(&self) -> Adaptive {
        let error = self.error;
        let iterations = self.iteration;
        let result = self.sum_results();
        Adaptive::new(result, error, iterations)
    }

    fn update_error_over_large_intervals(
        &mut self,
        iteration_interval: f64,
        previous_error: f64,
        iteration_error: f64,
    ) {
        self.error_over_large_intervals += -previous_error;

        if iteration_interval > self.smallest_interval {
            self.error_over_large_intervals += iteration_error;
        }
    }

    // XXX This function is entered if:
    //  - There is an error ANYWHERE in the bisection
    //      + If _no_ extrapolation has occurred, this goes to 115
    //  - The error output from extrapolation < extrapolation tolerance
    //  - If there is an error in the extrapolation
    fn go_to_100(&self) -> Result<Adaptive, Error> {
        if self.table.extrapolation_error() == f64::MAX {
            self.go_to_115()
        } else {
            // TODO
            Ok(Adaptive::empty())
        }
    }

    fn go_to_105(&self) -> Result<Adaptive, Error> {
        if self.table.extrapolation_error() / self.table.extrapolation_result().abs()
            > self.error / self.result.abs()
        {
            self.go_to_115()
        } else {
            self.go_to_110()
        }
    }

    // test on divergence
    fn go_to_110(&self) -> Result<Adaptive, Error> {
        let output = self.table.integral_estimate(self.iteration);
        let max_area = f64::max(self.table.extrapolation_result().abs(), self.result.abs());
        if !self.positive_integrand && max_area < 0.01 * self.initial_abs {
            return Ok(output);
        }

        let ratio = self.table.extrapolation_result() / self.result;

        if ratio < 0.01 || ratio > 100.0 || self.error > self.result.abs() {
            let kind = Kind::Divergence;
            return Err(Error::new(kind, output));
        }

        Ok(output)
    }

    fn go_to_115(&self) -> Result<Adaptive, Error> {
        let output = self.integral_estimate();
        Ok(output)
    }

    fn go_to_130(&self) -> Result<Adaptive, Error> {
        let output = self.table.integral_estimate(self.iteration);
        Ok(output)
    }

    fn error(&self) -> f64 {
        self.error
    }

    fn result(&self) -> f64 {
        self.result
    }
}

impl IntoIterator for Workspace {
    type Item = BasicInternal;
    type IntoIter = IntoIter<BasicInternal>;

    fn into_iter(self) -> Self::IntoIter {
        self.heap.into_iter()
    }
}

struct ExtrapolationTable {
    number: usize,
    results: [f64; 52],
    cached_value_0: f64,
    cached_value_1: f64,
    cached_value_2: f64,
    cached: Cached,
    res_ext: f64,
    err_ext: f64,
    ertest: f64,
    ktmin: usize,
}

#[derive(Debug, PartialEq)]
enum Cached {
    Zero,
    One,
    Two,
    All,
}

impl ExtrapolationTable {
    fn new() -> Self {
        let number = 0;
        let results = [0.0; 52];
        let cached_value_0 = 0.0;
        let cached_value_1 = 0.0;
        let cached_value_2 = 0.0;
        let cached = Cached::Zero;
        let res_ext = 0.0;
        let err_ext = f64::MAX;
        let ertest = f64::MAX;
        let ktmin = 0;

        Self {
            number,
            results,
            cached_value_0,
            cached_value_1,
            cached_value_2,
            cached,
            res_ext,
            err_ext,
            ertest,
            ktmin,
        }
    }

    fn initialise(initial: &BasicInternal) -> Self {
        let mut table = Self::new();
        table.res_ext = initial.result();
        table.append_table(initial.result());
        table
    }

    fn cache(&mut self, value: f64) -> &mut Self {
        match self.cached {
            Cached::All => {
                std::mem::swap(&mut self.cached_value_0, &mut self.cached_value_1);
                std::mem::swap(&mut self.cached_value_1, &mut self.cached_value_2);
                self.cached_value_2 = value;
            }
            Cached::Two => {
                self.cached_value_2 = value;
                self.cached = Cached::All;
            }
            Cached::One => {
                self.cached_value_1 = value;
                self.cached = Cached::Two;
            }
            Cached::Zero => {
                self.cached_value_0 = value;
                self.cached = Cached::One;
            }
        }
        self
    }

    fn append_table(&mut self, result: f64) -> &mut Self {
        self.results[self.number] = result;
        self.number += 1;
        self
    }

    // Adapted directly from GSL
    #[must_use]
    fn extrapolate(&mut self, value: f64) -> (f64, f64) {
        self.append_table(value);
        let n_current = self.number - 1;

        let current = self.results[n_current];

        let mut absolute = f64::MAX;
        let mut relative = 5.0 * f64::EPSILON * f64::abs(current);

        let new_element = n_current / 2;
        let n_original = n_current;
        let mut n_final = n_current;

        let mut result = current;
        let mut abserr = f64::MAX;

        if n_current < 2 {
            result = current;
            abserr = f64::max(absolute, relative);
            return (result, abserr);
        }

        self.results[n_current + 2] = self.results[n_current];
        self.results[n_current] = f64::MAX;

        for i in 0..new_element {
            let mut res = self.results[n_current - 2 * i + 2];
            let e0 = self.results[n_current - 2 * i - 2];
            let e1 = self.results[n_current - 2 * i - 1];
            let e2 = res;

            let e1abs = f64::abs(e1);
            let delta2 = e2 - e1;
            let err2 = f64::abs(delta2);
            let tol2 = f64::max(f64::abs(e2), e1abs) * f64::EPSILON;
            let delta3 = e1 - e0;
            let err3 = f64::abs(delta3);
            let tol3 = f64::max(e1abs, f64::abs(e0)) * f64::EPSILON;

            if err2 <= tol2 && err3 <= tol3 {
                // If e0, e1 and e2 are equal to within machine accuracy,
                // convergence is assumed.
                result = res;
                absolute = err2 + err3;
                relative = 5.0 * f64::EPSILON * f64::abs(res);
                abserr = f64::max(absolute, relative);
                return (result, abserr);
            }

            let e3 = self.results[n_current - 2 * i];
            self.results[n_current - 2 * i] = e1;
            let delta1 = e1 - e3;
            let err1 = f64::abs(delta1);
            let tol1 = f64::max(e1abs, f64::abs(e3)) * f64::EPSILON;

            // If two elements are very close to each other, omit a part
            // of the table by adjusting the value of n.
            if err1 <= tol1 || err2 <= tol2 || err3 <= tol3 {
                n_final = 2 * i;
                break;
            }

            let ss = (1.0 / delta1 + 1.0 / delta2) - 1.0 / delta3;

            // Test to detect irregular behaviour in the table, and eventually
            // omit a part of the table by adjusting the value of n.
            if f64::abs(ss * e1) <= 0.0001 {
                n_final = 2 * i;
                break;
            }

            // Compute a new element and eventually adjust the value of result.
            res = e1 + 1.0 / ss;
            self.results[n_current - 2 * i] = res;

            let error = err2 + f64::abs(res - e2) + err3;

            if error <= abserr {
                abserr = error;
                result = res;
            }
        }

        let limexp = 50 - 1;
        if n_final == limexp {
            n_final = 2 * (limexp / 2);
        }

        self.shift_table(n_original, n_final, new_element);
        self.number = n_final + 1;

        if let Cached::All = self.cached {
            abserr = f64::abs(result - self.cached_value_2)
                + f64::abs(result - self.cached_value_1)
                + f64::abs(result - self.cached_value_0);
        } else {
            abserr = f64::MAX;
        }

        self.cache(result);
        abserr = f64::max(abserr, 5.0 * f64::EPSILON * f64::abs(result));

        self.ktmin += 1;

        (result, abserr)
    }

    fn shift_table(&mut self, n_original: usize, n_final: usize, new_element: usize) -> &mut Self {
        if n_original % 2 == 1 {
            for i in 0..=new_element {
                self.results[1 + i * 2] = self.results[i * 2 + 3];
            }
        } else {
            for i in 0..=new_element {
                self.results[i * 2] = self.results[i * 2 + 2];
            }
        }

        if n_original != n_final {
            for i in 0..=n_final {
                self.results[i] = self.results[n_original - n_final + i];
            }
        }
        self
    }

    fn check_extrapolation_roundoff(
        &self,
        result: f64,
        error: f64,
        iterations: usize,
    ) -> Result<(), Error> {
        if self.ktmin > 5 && self.err_ext < 0.001 * error {
            let output = Adaptive::new(result, error, iterations);
            let kind = Kind::FailedToReachToleranceRoundoffExtrapolation;
            Err(Error::new(kind, output))
        } else {
            Ok(())
        }
    }

    fn integral_estimate(&self, iterations: usize) -> Adaptive {
        let result = self.res_ext;
        let error = self.err_ext;
        Adaptive::new(result, error, iterations)
    }

    #[inline]
    fn extrapolation_error(&self) -> f64 {
        self.err_ext
    }

    #[inline]
    fn extrapolation_result(&self) -> f64 {
        self.res_ext
    }
}

struct InfiniteInterval<I: Integrand> {
    function: I,
}

impl<I: Integrand> InfiniteInterval<I> {
    fn transform_evaluate(&self, t: &f64) -> f64 {
        let x = (1.0 - t) / t;
        let y = self.function.evaluate(&x) + self.function.evaluate(&(-x));
        y / (t.powi(2))
    }
}

impl<I: Integrand> Integrand for InfiniteInterval<I> {
    fn evaluate(&self, x: &f64) -> f64 {
        self.transform_evaluate(x)
    }
}

struct SemiInfiniteIntervalPositive<I: Integrand> {
    lower: f64,
    function: I,
}

impl<I: Integrand> SemiInfiniteIntervalPositive<I> {
    fn transform_evaluate(&self, t: &f64) -> f64 {
        let x = self.lower + (1.0 - t) / t;
        let y = self.function.evaluate(&x);
        y / (t.powi(2))
    }
}

impl<I: Integrand> Integrand for SemiInfiniteIntervalPositive<I> {
    fn evaluate(&self, x: &f64) -> f64 {
        self.transform_evaluate(x)
    }
}

struct SemiInfiniteIntervalNegative<I: Integrand> {
    upper: f64,
    function: I,
}

impl<I: Integrand> SemiInfiniteIntervalNegative<I> {
    fn transform_evaluate(&self, t: &f64) -> f64 {
        let x = self.upper - (1.0 - t) / t;
        let y = self.function.evaluate(&x);
        y / (t.powi(2))
    }
}

impl<I: Integrand> Integrand for SemiInfiniteIntervalNegative<I> {
    fn evaluate(&self, x: &f64) -> f64 {
        self.transform_evaluate(x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_swapping_cached() {
        let mut table = ExtrapolationTable::new();
        let value0 = 1.0;
        let value1 = 10.0;
        let value2 = 100.0;
        let value = 999.0;

        assert_eq!(table.cached, Cached::Zero);
        assert_eq!(table.cached_value_0, 0.0);
        assert_eq!(table.cached_value_1, 0.0);
        assert_eq!(table.cached_value_2, 0.0);

        table.cache(value0);

        assert_eq!(table.cached, Cached::One);
        assert_eq!(table.cached_value_0, 1.0);
        assert_eq!(table.cached_value_1, 0.0);
        assert_eq!(table.cached_value_2, 0.0);

        table.cache(value1);

        assert_eq!(table.cached, Cached::Two);
        assert_eq!(table.cached_value_0, 1.0);
        assert_eq!(table.cached_value_1, 10.0);
        assert_eq!(table.cached_value_2, 0.0);

        table.cache(value2);

        assert_eq!(table.cached, Cached::All);
        assert_eq!(table.cached_value_0, 1.0);
        assert_eq!(table.cached_value_1, 10.0);
        assert_eq!(table.cached_value_2, 100.0);

        table.cache(value);

        assert_eq!(table.cached, Cached::All);
        assert_eq!(table.cached_value_0, 10.0);
        assert_eq!(table.cached_value_1, 100.0);
        assert_eq!(table.cached_value_2, 999.0);
    }

    #[test]
    fn test_intitialising() {
        let error = 0.0;
        let result = 0.0;
        let result_abs = 0.0;
        let result_asc = 0.0;
        let lower = 0.0;
        let upper = 0.0;
        let value = BasicInternal {
            error,
            result,
            result_abs,
            result_asc,
            lower,
            upper,
        };

        let table = ExtrapolationTable::initialise(&value);

        assert!(table.number == 1);
        assert!((table.results[0] - value.result()) < 1e-10);
    }
}
