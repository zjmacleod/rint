use std::collections::BinaryHeap;

use crate::integration::adaptive::{self, Error, Kind, Workspace};
use crate::integration::basic::BasicInternal;
use crate::integration::{Adaptive, ErrorBound, GaussKronrodBasic};
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

    pub fn integrate(&self) -> Result<Adaptive, Error> {
        todo!()
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

struct ExtrapolationTable {
    number: usize,
    results: [f64; 52],
    cached_value_0: f64,
    cached_value_1: f64,
    cached_value_2: f64,
    cached: Cached,
    state: ExtrapolationState,
}

enum ExtrapolationState {
    On,
    Off,
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
        let state = ExtrapolationState::Off;

        Self {
            number,
            results,
            cached_value_0,
            cached_value_1,
            cached_value_2,
            cached,
            state,
        }
    }

    fn initialise(result: f64) -> Self {
        let mut table = Self::new();
        table.append_table(result);
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
    fn extrapolate(&mut self, value: f64) -> Data {
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
            return Data::new(result, abserr);
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
                return Data::new(result, abserr);
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

        Data::new(result, abserr)
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
}

struct Data {
    result: f64,
    error: f64,
}

impl Data {
    fn new(result: f64, error: f64) -> Self {
        Self { result, error }
    }

    fn result(&self) -> f64 {
        self.result
    }

    fn error(&self) -> f64 {
        self.error
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
        let value = 999.0;

        let table = ExtrapolationTable::initialise(value);

        assert!(table.number == 1);
        assert!((table.results[0] - value) < 1e-10);
    }
}

pub(crate) fn test_positivity(result: f64, result_abs: f64) -> bool {
    result.abs() >= (1.0 - 50.0 * f64::EPSILON) * result_abs
}
