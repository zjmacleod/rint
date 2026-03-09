mod fully_symmetric_07;
mod fully_symmetric_09;
mod fully_symmetric_09_2d;
mod fully_symmetric_11_3d;
mod fully_symmetric_13_2d;

use crate::multi::generator::Generator;
use crate::multi::two_pow_n_f64;
use crate::InitialisationError;
use crate::ScalarF64;

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct Rule<const NDIM: usize, const FINAL: usize, const TOTAL: usize> {
    initial_data: [Data<NDIM>; 3],
    final_data: [Data<NDIM>; FINAL],
    scales_norms: [ScalesNorms<TOTAL>; 3],
    basic_error_coeff: BasicErrorCoeff,
    adaptive_error_coeff: AdaptiveErrorCoeff,
    evaluations: usize,
    ratio: f64,
}

impl<const NDIM: usize, const FINAL: usize, const TOTAL: usize> Rule<NDIM, FINAL, TOTAL> {
    pub(crate) const fn initial_data(&self) -> &[Data<NDIM>; 3] {
        &self.initial_data
    }

    pub(crate) const fn final_data(&self) -> &[Data<NDIM>; FINAL] {
        &self.final_data
    }

    pub(crate) const fn scales_norms(&self) -> &[ScalesNorms<TOTAL>] {
        &self.scales_norms
    }

    pub(crate) const fn basic_error_coeff(&self) -> &BasicErrorCoeff {
        &self.basic_error_coeff
    }

    pub(crate) const fn adaptive_error_coeff(&self) -> &AdaptiveErrorCoeff {
        &self.adaptive_error_coeff
    }

    pub(crate) const fn evaluations(&self) -> usize {
        self.evaluations
    }

    pub(crate) const fn ratio(&self) -> f64 {
        self.ratio
    }
}

/// A 7-point fully-symmetric integration rule valid for `2 <= NDIM <= 15`.
///
/// This is the recommended rule to use when a significant amount of adaptivity is required.
pub type Rule07<const NDIM: usize> = Rule<NDIM, 3, 6>;

impl<const NDIM: usize> Rule07<NDIM> {
    /// Generate a fully-symmetric 7-point integration rule.
    ///
    /// # Errors
    /// Will fail if `NDIM < 2` or `NDIM > 15`.
    pub const fn generate() -> Result<Self, InitialisationError> {
        fully_symmetric_07::generate_rule::<NDIM>()
    }
}

/// A 9-point fully-symmetric integration rule valid for `3 <= NDIM <= 15`.
pub type Rule09<const NDIM: usize> = Rule<NDIM, 6, 9>;

impl<const NDIM: usize> Rule09<NDIM> {
    /// Generate a fully-symmetric 9-point integration rule for NDIM > 3 dimensional integration.
    ///
    /// # Errors
    /// Will fail if `NDIM < 3` or `NDIM > 15`.
    pub const fn generate() -> Result<Self, InitialisationError> {
        fully_symmetric_09::generate_rule::<NDIM>()
    }
}

/// A 9-point fully-symmetric integration rule valid for `NDIM == 2`.
pub type Rule09N2 = Rule<2, 5, 8>;

impl Rule09N2 {
    /// Generate a fully-symmetric 9-point integration rule for NDIM == 3 dimensional integration.
    #[must_use]
    pub const fn generate() -> Self {
        fully_symmetric_09_2d::generate_rule()
    }
}

/// A 11-point fully-symmetric integration rule valid for `NDIM == 3`.
pub type Rule11 = Rule<3, 10, 13>;

impl Rule11 {
    /// Generate a fully-symmetric 11-point integration rule for `NDIM = 3` dimensional integration.
    #[must_use]
    pub const fn generate() -> Self {
        fully_symmetric_11_3d::generate_rule()
    }
}

/// A 13-point fully-symmetric integration rule valid for `NDIM == 2`.
pub type Rule13 = Rule<2, 11, 14>;

impl Rule13 {
    /// Generate a fully-symmetric 13-point integration rule for `NDIM = 2` dimensional integration.
    #[must_use]
    pub const fn generate() -> Self {
        fully_symmetric_13_2d::generate_rule()
    }
}

const ADAPTIVE_ERROR_COEFF: AdaptiveErrorCoeff = AdaptiveErrorCoeff::new(0.5, 0.25);

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub(crate) struct Data<const NDIM: usize> {
    generator: Generator<NDIM>,
    weight: f64,
    null1: f64,
    null2: f64,
    null3: f64,
    null4: f64,
}

impl<const NDIM: usize> Data<NDIM> {
    pub(crate) const fn new(generator: Generator<NDIM>, weights: [f64; 5]) -> Self {
        let [weight, null1, null2, null3, null4] = weights;
        Self {
            generator,
            weight,
            null1,
            null2,
            null3,
            null4,
        }
    }

    pub(crate) const fn generator(&self) -> &Generator<NDIM> {
        &self.generator
    }

    pub(crate) fn get_first_value(&self) -> f64 {
        // unwrap is fine as generators are never empty
        *self.generator.generator().first().unwrap()
    }

    pub(crate) const fn weight(&self) -> f64 {
        self.weight
    }

    pub(crate) const fn null1(&self) -> f64 {
        self.null1
    }

    pub(crate) const fn null2(&self) -> f64 {
        self.null2
    }

    pub(crate) const fn null3(&self) -> f64 {
        self.null3
    }

    pub(crate) const fn null4(&self) -> f64 {
        self.null4
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub(crate) struct ScalesNorms<const TOTAL: usize> {
    scales: [f64; TOTAL],
    norms: [f64; TOTAL],
}

impl<const TOTAL: usize> ScalesNorms<TOTAL> {
    const fn new(scales: [f64; TOTAL], norms: [f64; TOTAL]) -> Self {
        Self { scales, norms }
    }

    const fn scales(&self) -> &[f64] {
        &self.scales
    }

    const fn norms(&self) -> &[f64] {
        &self.norms
    }

    pub(crate) fn max_consecutive_null<T: ScalarF64>(&self, null1: T, null2: T) -> T {
        let scales = self.scales();
        let norms = self.norms();

        scales
            .iter()
            .zip(norms)
            .map(|(s, n)| (null1 * s + null2) * n)
            .fold(T::zero(), |a, b| if a.abs() < b.abs() { b } else { a })
    }
}

pub(crate) const fn scales_norms<const NDIM: usize, const TOTAL: usize>(
    weights: &[[f64; 5]; TOTAL],
    rule_points: [f64; TOTAL],
) -> [ScalesNorms<TOTAL>; 3] {
    let two_ndim = two_pow_n_f64(NDIM);

    let mut scales = [[0f64; TOTAL]; 3];
    let mut norms = [[0f64; TOTAL]; 3];
    let mut we = [0f64; 14];

    let mut i = 0;
    while i < TOTAL {
        let mut k = 0;
        while k < 3 {
            scales[k][i] = if weights[i][k + 1] != 0.0 {
                -weights[i][k + 2] / weights[i][k + 1]
            } else {
                100.0
            };
            let mut j = 0;
            while j < TOTAL {
                we[j] = weights[j][k + 2] + scales[k][i] * weights[j][k + 1];
                j += 1;
            }
            norms[k][i] = 0.0;
            let mut j = 0;
            while j < TOTAL {
                let weabs = we[j].abs();

                norms[k][i] += rule_points[j] * weabs;
                j += 1;
            }
            norms[k][i] = two_ndim / norms[k][i];

            k += 1;
        }
        i += 1;
    }

    let scales_norms0 = ScalesNorms::new(scales[0], norms[0]);
    let scales_norms1 = ScalesNorms::new(scales[1], norms[1]);
    let scales_norms2 = ScalesNorms::new(scales[2], norms[2]);

    [scales_norms0, scales_norms1, scales_norms2]
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub(crate) struct BasicErrorCoeff {
    c1: f64,
    c2: f64,
    c3: f64,
    c4: f64,
}

impl BasicErrorCoeff {
    pub(crate) const fn new(c1: f64, c2: f64, c3: f64, c4: f64) -> Self {
        Self { c1, c2, c3, c4 }
    }
    #[inline]
    pub(crate) const fn c1(&self) -> f64 {
        self.c1
    }

    #[inline]
    pub(crate) const fn c2(&self) -> f64 {
        self.c2
    }

    #[inline]
    pub(crate) const fn c3(&self) -> f64 {
        self.c3
    }

    #[inline]
    pub(crate) const fn c4(&self) -> f64 {
        self.c4
    }
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub(crate) struct AdaptiveErrorCoeff {
    c5: f64,
    c6: f64,
}

impl AdaptiveErrorCoeff {
    pub(crate) const fn new(c5: f64, c6: f64) -> Self {
        Self { c5, c6 }
    }

    #[inline]
    pub(crate) const fn c5(&self) -> f64 {
        self.c5
    }

    #[inline]
    pub(crate) const fn c6(&self) -> f64 {
        self.c6
    }
}

#[cfg(test)]
mod util {
    use super::*;

    pub(crate) fn rel_or_abs_diff(a: f64, b: f64) -> f64 {
        if a == 0.0 {
            (a - b).abs()
        } else {
            (a - b).abs() / a.abs()
        }
    }

    pub(crate) fn assert_check_vec_tol<const WL: usize, const TY: usize>(
        calc: &[[f64; TY]; WL],
        should_be: &[[f64; TY]; WL],
        tol: f64,
    ) {
        for (x, y) in calc.iter().zip(should_be.iter()) {
            for (a, b) in x.iter().zip(y.iter()) {
                let val = rel_or_abs_diff(*a, *b);

                assert!(val < tol);
            }
        }
    }

    pub(crate) fn assert_check_slice_tol(calc: &[f64], should_be: &[f64], tol: f64) {
        for (x, y) in calc.iter().zip(should_be.iter()) {
            let val = rel_or_abs_diff(*x, *y);

            assert!(val < tol);
        }
    }

    #[allow(clippy::similar_names)]
    pub(crate) fn assert_check_vec_data_tol<const WL: usize, const TY: usize>(
        calc: &[Data<TY>; WL],
        should_be: &[Data<TY>; WL],
        tol: f64,
    ) {
        for (x, y) in calc.iter().zip(should_be.iter()) {
            let genx = x.generator().generator();
            let geny = y.generator().generator();
            for (gx, gy) in genx.iter().zip(geny.iter()) {
                let val = rel_or_abs_diff(*gx, *gy);
                assert!(val < tol);
            }

            let wx = x.weight();
            let wy = y.weight();
            let val = rel_or_abs_diff(wx, wy);
            assert!(val < tol);

            let n1x = x.null1();
            let n1y = y.null1();
            let val = rel_or_abs_diff(n1x, n1y);
            assert!(val < tol);

            let n2x = x.null2();
            let n2y = y.null2();
            let val = rel_or_abs_diff(n2x, n2y);
            assert!(val < tol);

            let n3x = x.null3();
            let n3y = y.null3();
            let val = rel_or_abs_diff(n3x, n3y);
            assert!(val < tol);

            let n4x = x.null4();
            let n4y = y.null4();
            let val = rel_or_abs_diff(n4x, n4y);
            assert!(val < tol);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn error_when_wrong_dimension() {
        {
            const NDIM: usize = 1;
            let err = Rule07::<NDIM>::generate();
            assert!(err.is_err());
        }
        {
            const NDIM: usize = 2;
            let err = Rule07::<NDIM>::generate();
            assert!(err.is_ok());
        }
        {
            const NDIM: usize = 8;
            let err = Rule07::<NDIM>::generate();
            assert!(err.is_ok());
        }
        {
            const NDIM: usize = 15;
            let err = Rule07::<NDIM>::generate();
            assert!(err.is_ok());
        }
        {
            const NDIM: usize = 16;
            let err = Rule07::<NDIM>::generate();
            assert!(err.is_err());
        }

        {
            const NDIM: usize = 1;
            let err = Rule09::<NDIM>::generate();
            assert!(err.is_err());
        }
        {
            const NDIM: usize = 2;
            let err = Rule09::<NDIM>::generate();
            assert!(err.is_err());
        }
        {
            const NDIM: usize = 8;
            let err = Rule09::<NDIM>::generate();
            assert!(err.is_ok());
        }
        {
            const NDIM: usize = 15;
            let err = Rule09::<NDIM>::generate();
            assert!(err.is_ok());
        }
        {
            const NDIM: usize = 16;
            let err = Rule09::<NDIM>::generate();
            assert!(err.is_err());
        }
    }
}
