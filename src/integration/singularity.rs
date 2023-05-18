use std::collections::binary_heap::BinaryHeap;

use crate::integration::basic::BasicInternal;
use crate::integration::{subinterval_too_small, Basic, Error, ErrorBound, IntegralEstimate, Kind};
use crate::rule::{GaussKronrod15, GaussKronrod21, Rule};
use crate::Integrand;

/// An integral with singularities to be evaluated with an adaptive Gauss-Kronrod quadrature.
///
/// The user constructs a `function` implementing [`Integrand`], provides `upper`
/// and `lower` integration limits, and provides an `error_bound`, which can be
/// [`ErrorBound::Absolute`] to work to a specified absolute error,
/// [`ErrorBound::Relative`] to work to a specified relative error,
/// or [`ErrorBound::Either`] to return a result as soon as _either_ the relative
/// or absolute error bound has been satisfied.
pub struct AdaptiveSingularity<I, R>
where
    I: Integrand,
    R: Rule,
{
    lower: f64,
    upper: f64,
    error_bound: ErrorBound,
    rule: R,
    function: I,
    max_iterations: usize,
}

impl<I, R> AdaptiveSingularity<I, R>
where
    I: Integrand,
    R: Rule,
{
    /// Create a new [`AdaptiveSingularity`].
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
    fn new(
        lower: f64,
        upper: f64,
        error_bound: ErrorBound,
        rule: R,
        function: I,
        max_iterations: usize,
    ) -> Result<Self, Error> {
        match error_bound {
            ErrorBound::Absolute(v) => {
                if v <= 0.0 {
                    let kind = Kind::RelativeBoundNegativeOrZero;
                    return Err(Error::unevaluated(kind));
                }
            }
            ErrorBound::Relative(v) => {
                if v < 50.0 * f64::EPSILON {
                    let kind = Kind::AbsoluteBoundTooSmall;
                    return Err(Error::unevaluated(kind));
                }
            }
            ErrorBound::Either { absolute, relative } => {
                if absolute <= 0.0 && relative < 50.0 * f64::EPSILON {
                    let kind = Kind::InvalidTolerance;
                    return Err(Error::unevaluated(kind));
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
    ) -> Result<Option<IntegralEstimate>, Error> {
        let tolerance = self.error_bound.tolerance(initial.result());
        let roundoff = Self::roundoff(initial.result_abs());

        if initial.error() <= roundoff && initial.error() > tolerance {
            let output = IntegralEstimate::from_basic(initial, 1, self.rule.evaluations());
            let kind = Kind::RoundoffErrorDetected;

            Err(Error::new(kind, output))
        } else if (initial.error() <= tolerance
            && initial.error().to_bits() != initial.result_asc().to_bits())
            || initial.error() == 0.0
        {
            let output = IntegralEstimate::from_basic(initial, 1, self.rule.evaluations());

            Ok(Some(output))
        } else if self.max_iterations == 1 {
            let output = IntegralEstimate::from_basic(initial, 1, self.rule.evaluations());
            let kind = Kind::MaximumIterationsReached;

            Err(Error::new(kind, output))
        } else {
            Ok(None)
        }
    }

    /// # Errors
    pub fn integrate(&self) -> Result<IntegralEstimate, Error> {
        let initial =
            Basic::new(self.lower, self.upper, self.rule, &self.function).integrate_internal();

        if let Some(output) = self.check_initial_integration(&initial)? {
            return Ok(output);
        }

        let mut workspace = self.initialise_workspace(initial);

        while workspace.iteration < self.max_iterations {
            let previous = workspace.retrieve_largest_error()?;

            let current_interval = previous.abs_interval_length();

            let [lower, upper] = previous.bisect(&self.function, self.rule);

            let previous_error = previous.error();
            let iteration_error = lower.error() + upper.error();

            let [result, error] = workspace.improved_result_error(&previous, &lower, &upper);

            let iteration_tolerance = self.error_bound.tolerance(result);

            workspace.update(lower, upper);

            if error <= iteration_tolerance {
                return workspace.compute_result();
            }

            if workspace.is_err() {
                break;
            }

            if workspace.iteration >= self.max_iterations - 1 {
                workspace.set_error_kind(Kind::MaximumIterationsReached);
                break;
            }

            if workspace.iteration == 2 {
                workspace.smallest_interval *= 0.375;
                workspace.large_interval_error = workspace.error;
                workspace.table.tolerance = iteration_tolerance;
                workspace.table.append_table(workspace.result);
                continue;
            }

            workspace.update_large_interval_error(
                current_interval,
                previous_error,
                iteration_error,
            );

            // peek the next interval with the largest error and check if it
            // is the smallest interval.
            if !workspace.extrapolate {
                if let Some(next) = workspace.peek() {
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

            if !workspace.table.error_detected
                && workspace.large_interval_error > workspace.table.tolerance
            {
                if workspace.remaining_large_intervals() {
                    continue;
                }

                // If we reach here then _all_ the results are now in the store
                // so we mem::swap the store and the heap and continue on to the
                // extrapolation.
                std::mem::swap(&mut workspace.heap, &mut workspace.store);
            }

            let [ext_result, ext_error] = workspace.table.extrapolate(result);

            workspace.check_convergence();

            if ext_error < workspace.table.error {
                let ext_tolerance = self.error_bound.tolerance(ext_result);
                workspace.update_table_values(ext_result, ext_error, ext_tolerance);
                if workspace.table.error <= workspace.table.tolerance {
                    break;
                }
            }

            if let Some(Kind::DoesNotConverge) = workspace.error_kind {
                break;
            }

            workspace.prepare_next_iteration();
        }

        workspace.check_error_and_compute()
    }

    /// Return the value of the `upper` integration limit.
    pub fn upper(&self) -> f64 {
        self.upper
    }

    /// Return the value of the `lower` integration limit.
    pub fn lower(&self) -> f64 {
        self.lower
    }

    fn initialise_workspace(&self, initial: BasicInternal) -> Workspace {
        let mut heap = BinaryHeap::with_capacity(2 * self.max_iterations + 1);

        let iteration = 1;
        let result = initial.result();
        let error = initial.error();
        let large_interval_error = error;
        let table = ExtrapolationTable::initialise(&initial);
        let smallest_interval = initial.abs_interval_length();
        let extrapolate = false;
        let store = BinaryHeap::with_capacity(2 * self.max_iterations + 1);
        let error_kind = None;
        let initial_absolute_result = initial.result_abs();
        let positive_integrand = initial.positivity();
        let roundoff_count = 0;
        let roundoff_on_high_iteration_count = 0;
        let function_evaluations_per_integration = self.rule.evaluations();

        heap.push(initial);

        Workspace {
            heap,
            iteration,
            result,
            error,
            large_interval_error,
            table,
            smallest_interval,
            extrapolate,
            store,
            error_kind,
            initial_absolute_result,
            positive_integrand,
            roundoff_count,
            roundoff_on_high_iteration_count,
            function_evaluations_per_integration,
        }
    }
}

impl<I> AdaptiveSingularity<I, GaussKronrod21>
where
    I: Integrand,
{
    /// # Errors
    pub fn general(
        lower: f64,
        upper: f64,
        error_bound: ErrorBound,
        function: I,
        max_iterations: usize,
    ) -> Result<Self, Error> {
        let rule = GaussKronrod21;
        Self::new(lower, upper, error_bound, rule, function, max_iterations)
    }
}

pub struct InfiniteInterval<I: Integrand> {
    function: I,
}

impl<I: Integrand> InfiniteInterval<I> {
    fn new(function: I) -> Self {
        Self { function }
    }

    // TODO the calculation of the number of function evaluations is wrong for this type, since
    // there are two evaluations per call to .evaluate()
    fn transform_evaluate(&self, t: f64) -> f64 {
        let x = (1.0 - t) / t;
        let y = self.function.evaluate(x) + self.function.evaluate(-x);
        (y / t) / t
    }
}

impl<I: Integrand> Integrand for InfiniteInterval<I> {
    fn evaluate(&self, x: f64) -> f64 {
        self.transform_evaluate(x)
    }
}

impl<I> AdaptiveSingularity<InfiniteInterval<I>, GaussKronrod15>
where
    I: Integrand,
{
    // TODO the calculation of the number of function evaluations is wrong for this type, since
    // there are two evaluations per call to .evaluate()
    /// # Errors
    pub fn infinite(
        error_bound: ErrorBound,
        function: I,
        max_iterations: usize,
    ) -> Result<Self, Error> {
        let rule = GaussKronrod15;
        let transformed = InfiniteInterval::new(function);
        Self::new(0.0, 1.0, error_bound, rule, transformed, max_iterations)
    }
}

pub struct SemiInfiniteIntervalPositive<I: Integrand> {
    lower: f64,
    function: I,
}

impl<I: Integrand> SemiInfiniteIntervalPositive<I> {
    fn new(function: I, lower: f64) -> Self {
        Self { lower, function }
    }

    fn transform_evaluate(&self, t: f64) -> f64 {
        let x = self.lower + (1.0 - t) / t;
        let y = self.function.evaluate(x);
        y / (t.powi(2))
    }
}

impl<I: Integrand> Integrand for SemiInfiniteIntervalPositive<I> {
    fn evaluate(&self, x: f64) -> f64 {
        self.transform_evaluate(x)
    }
}

impl<I> AdaptiveSingularity<SemiInfiniteIntervalPositive<I>, GaussKronrod15>
where
    I: Integrand,
{
    /// # Errors
    pub fn semi_infinite_positive(
        lower: f64,
        error_bound: ErrorBound,
        function: I,
        max_iterations: usize,
    ) -> Result<Self, Error> {
        let rule = GaussKronrod15;
        let transformed = SemiInfiniteIntervalPositive::new(function, lower);
        Self::new(0.0, 1.0, error_bound, rule, transformed, max_iterations)
    }
}

pub struct SemiInfiniteIntervalNegative<I: Integrand> {
    upper: f64,
    function: I,
}

impl<I: Integrand> SemiInfiniteIntervalNegative<I> {
    fn new(function: I, upper: f64) -> Self {
        Self { upper, function }
    }

    fn transform_evaluate(&self, t: f64) -> f64 {
        let x = self.upper - (1.0 - t) / t;
        let y = self.function.evaluate(x);
        y / (t.powi(2))
    }
}

impl<I: Integrand> Integrand for SemiInfiniteIntervalNegative<I> {
    fn evaluate(&self, x: f64) -> f64 {
        self.transform_evaluate(x)
    }
}

impl<I> AdaptiveSingularity<SemiInfiniteIntervalNegative<I>, GaussKronrod15>
where
    I: Integrand,
{
    /// # Errors
    pub fn semi_infinite_negative(
        upper: f64,
        error_bound: ErrorBound,
        function: I,
        max_iterations: usize,
    ) -> Result<Self, Error> {
        let rule = GaussKronrod15;
        let transformed = SemiInfiniteIntervalNegative::new(function, upper);
        Self::new(0.0, 1.0, error_bound, rule, transformed, max_iterations)
    }
}

struct Workspace {
    heap: BinaryHeap<BasicInternal>,
    iteration: usize,
    result: f64,
    error: f64,
    large_interval_error: f64,
    table: ExtrapolationTable,
    smallest_interval: f64,
    extrapolate: bool,
    store: BinaryHeap<BasicInternal>,
    error_kind: Option<Kind>,
    initial_absolute_result: f64,
    positive_integrand: bool,
    roundoff_count: usize,
    roundoff_on_high_iteration_count: usize,
    function_evaluations_per_integration: usize,
}

impl Workspace {
    fn retrieve_largest_error(&mut self) -> Result<BasicInternal, Error> {
        self.iteration += 1;
        if let Some(previous) = self.pop() {
            Ok(previous)
        } else {
            let kind = Kind::UninitialisedWorkspace;
            Err(Error::unevaluated(kind))
        }
    }

    fn pop(&mut self) -> Option<BasicInternal> {
        self.heap.pop()
    }

    fn push(&mut self, integral: BasicInternal) {
        self.heap.push(integral);
    }

    fn peek(&self) -> Option<&BasicInternal> {
        self.heap.peek()
    }

    fn store(&mut self, integral: BasicInternal) {
        self.store.push(integral);
    }

    fn sum_results(&self) -> f64 {
        let store = self.store.iter().fold(0.0f64, |a, v| a + v.result());
        let heap = self.heap.iter().fold(0.0f64, |a, v| a + v.result());
        store + heap
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
                if self.extrapolate {
                    self.table.roundoff_count += 1;
                } else {
                    self.roundoff_count += 1;
                }
            } else if self.iteration > 10 && new_error > prev_error {
                self.roundoff_on_high_iteration_count += 1;
            }
        }

        self.result = (self.result + new_result) - prev_result;
        self.error = (self.error + new_error) - prev_error;

        [self.result, self.error]
    }

    fn integral_estimate(&self) -> IntegralEstimate {
        let result = self.sum_results();
        let error = self.error;
        let iterations = self.iteration;
        let function_evaluations = (2 * iterations - 1) * self.function_evaluations_per_integration;
        IntegralEstimate::new()
            .with_result(result)
            .with_error(error)
            .with_iterations(iterations)
            .with_function_evaluations(function_evaluations)
    }

    fn prepare_next_iteration(&mut self) {
        self.extrapolate = false;
        self.large_interval_error = self.error;
        self.heap.append(&mut self.store);
        self.smallest_interval *= 0.5;
    }

    fn compute_result(self) -> Result<IntegralEstimate, Error> {
        let output = self.integral_estimate();

        if let Some(kind) = self.error_kind {
            Err(Error::new(kind, output))
        } else {
            Ok(output)
        }
    }

    fn compute_extrapolated_result(self) -> Result<IntegralEstimate, Error> {
        let function_evaluations =
            (2 * self.iteration - 1) * self.function_evaluations_per_integration;
        let output = self
            .table
            .integral_estimate(self.iteration, function_evaluations);

        if let Some(kind) = self.error_kind {
            Err(Error::new(kind, output))
        } else {
            Ok(output)
        }
    }

    fn check_error_and_compute(mut self) -> Result<IntegralEstimate, Error> {
        if self.table.error == f64::MAX {
            return self.compute_result();
        }

        if self.is_err() || self.table.is_err() {
            if self.table.is_err() {
                self.table.error += self.table.correction;
            }

            if !self.is_err() {
                self.set_error_kind(Kind::DoesNotConverge);
            }

            if self.table.result != 0.0 && self.result != 0.0 {
                if self.table.error / self.table.result.abs() > self.error / self.result.abs() {
                    return self.compute_result();
                }
            } else if self.table.error > self.error {
                return self.compute_result();
            } else if self.result == 0.0 {
                return self.compute_extrapolated_result();
            }
        }

        let max_area = f64::max(self.table.result.abs(), self.result.abs());

        if !self.positive_integrand && max_area < 0.01 * self.initial_absolute_result {
            return self.compute_extrapolated_result();
        }

        let ratio = self.table.result / self.result;

        if !(0.01..=100.0).contains(&ratio) || self.error > self.result.abs() {
            self.set_error_kind(Kind::DivergentOrSlowlyConverging);
        }

        self.compute_extrapolated_result()
    }

    fn check_roundoff(&mut self) {
        if self.roundoff_count + self.table.roundoff_count >= 10
            || self.roundoff_on_high_iteration_count >= 20
        {
            self.set_error_kind(Kind::RoundoffErrorDetected);
        }

        if self.table.roundoff_count >= 5 {
            self.table.error_detected = true;
        }
    }

    fn check_convergence(&mut self) {
        if self.table.ktmin > 5 && self.table.error < 0.001 * self.error {
            self.set_error_kind(Kind::DoesNotConverge);
        }
    }

    fn remaining_large_intervals(&mut self) -> bool {
        while let Some(next) = self.peek() {
            if next.abs_interval_length() > self.smallest_interval {
                return true;
            }
            if let Some(next) = self.pop() {
                self.store(next);
            };
        }
        false
    }

    fn is_err(&self) -> bool {
        self.error_kind.is_some()
    }

    fn set_error_kind(&mut self, kind: Kind) {
        self.error_kind = Some(kind);
    }

    fn update(&mut self, lower: BasicInternal, upper: BasicInternal) {
        let lower_limit = lower.lower();
        let upper_limit = upper.upper();
        let midpoint = (lower_limit + upper_limit) * 0.5;

        self.check_roundoff();

        if subinterval_too_small(lower_limit, midpoint, upper_limit) {
            self.set_error_kind(Kind::BadIntegrandBehaviour {
                lower: lower_limit,
                upper: upper_limit,
            });
        }

        self.push(lower);
        self.push(upper);
    }

    fn update_large_interval_error(
        &mut self,
        iteration_interval: f64,
        previous_error: f64,
        iteration_error: f64,
    ) {
        self.large_interval_error += -previous_error;

        let half_iteration_interval = iteration_interval * 0.5;
        if half_iteration_interval > self.smallest_interval {
            self.large_interval_error += iteration_error;
        }
    }

    fn update_table_values(&mut self, ext_result: f64, ext_error: f64, tolerance: f64) {
        self.table.ktmin = 0;
        self.table.error = ext_error;
        self.table.result = ext_result;
        self.table.correction = self.large_interval_error;
        self.table.tolerance = tolerance;
    }
}

struct ExtrapolationTable {
    count: usize,
    results: [f64; 52],
    cached: Cached,
    error_detected: bool,
    result: f64,
    error: f64,
    tolerance: f64,
    ktmin: usize,
    roundoff_count: usize,
    correction: f64,
}

#[derive(Debug, PartialEq)]
enum Cached {
    Empty,
    One(f64),
    Two(f64, f64),
    Full(f64, f64, f64),
}

impl ExtrapolationTable {
    fn new() -> Self {
        let count = 0;
        let results = [0.0; 52];
        let cached = Cached::Empty;
        let error_detected = false;
        let result = 0.0;
        let error = f64::MAX;
        let tolerance = f64::MAX;
        let ktmin = 0;
        let roundoff_count = 0;
        let correction = 0.0;

        Self {
            count,
            results,
            cached,
            error_detected,
            result,
            error,
            tolerance,
            ktmin,
            roundoff_count,
            correction,
        }
    }

    fn initialise(initial: &BasicInternal) -> Self {
        let mut table = Self::new();
        table.result = initial.result();
        table.append_table(initial.result());
        table
    }

    fn is_err(&self) -> bool {
        self.error_detected
    }

    fn cache(&mut self, value: f64) -> &mut Self {
        match self.cached {
            Cached::Full(_, a, b) | Cached::Two(a, b) => {
                self.cached = Cached::Full(a, b, value);
            }
            Cached::One(a) => {
                self.cached = Cached::Two(a, value);
            }
            Cached::Empty => {
                self.cached = Cached::One(value);
            }
        }
        self
    }

    fn append_table(&mut self, result: f64) -> &mut Self {
        self.results[self.count] = result;
        self.count += 1;
        self
    }

    // Adapted directly from GSL
    #[must_use]
    fn extrapolate(&mut self, value: f64) -> [f64; 2] {
        self.append_table(value);
        let n_current = self.count - 1;

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
            return [result, abserr];
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
                return [result, abserr];
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

            let ss = ((1.0 / delta1) + (1.0 / delta2)) - (1.0 / delta3);

            // Test to detect irregular behaviour in the table, and eventually
            // omit a part of the table by adjusting the value of n.
            if f64::abs(ss * e1) <= 0.0001 {
                n_final = 2 * i;
                break;
            }

            // Compute a new element and eventually adjust the value of result.
            res = e1 + (1.0 / ss);
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
        self.count = n_final + 1;

        if let Cached::Full(a, b, c) = self.cached {
            abserr = f64::abs(result - c) + f64::abs(result - b) + f64::abs(result - a);
        } else {
            abserr = f64::MAX;
        }

        self.cache(result);
        abserr = f64::max(abserr, 5.0 * f64::EPSILON * f64::abs(result));

        self.ktmin += 1;

        [result, abserr]
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

    fn integral_estimate(
        &self,
        iterations: usize,
        function_evaluations: usize,
    ) -> IntegralEstimate {
        let result = self.result;
        let error = self.error;
        IntegralEstimate::new()
            .with_result(result)
            .with_error(error)
            .with_iterations(iterations)
            .with_function_evaluations(function_evaluations)
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

        assert_eq!(table.cached, Cached::Empty);

        table.cache(value0);

        assert_eq!(table.cached, Cached::One(1.0));

        table.cache(value1);

        assert_eq!(table.cached, Cached::Two(1.0, 10.0));

        table.cache(value2);

        assert_eq!(table.cached, Cached::Full(1.0, 10.0, 100.0));

        table.cache(value);

        assert_eq!(table.cached, Cached::Full(10.0, 100.0, 999.0));
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

        assert!(table.count == 1);
        assert!((table.results[0] - value.result()) < 1e-10);
    }
}
