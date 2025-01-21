mod fully_symmetric_07;
mod fully_symmetric_11_3d;
mod fully_symmetric_13_2d;

use num_traits::Zero;

use crate::multi::generator::Generator;
use crate::multi::geometry::Geometry;
use crate::multi::two_pow_n_f64;
use crate::Limits;
use crate::MultiDimensionalIntegrand;
use crate::ScalarF64;

const LAMBDA0: f64 = 6.860_757_975_617_562_914_002_852e-1;
const LAMBDA1: f64 = 9.559_073_158_045_390_123_857_208e-1;
const LAMBDA2: f64 = 4.060_571_747_382_397_355_996_069e-1;
const LAMBDAP: f64 = 7.500_000_000_000_000_000_000_000e-1;

const LAMBDA0_SQ: f64 = LAMBDA0 * LAMBDA0;
const LAMBDA1_SQ: f64 = LAMBDA1 * LAMBDA1;
const LAMBDA2_SQ: f64 = LAMBDA2 * LAMBDA2;
const LAMBDAP_SQ: f64 = LAMBDAP * LAMBDAP;

const ADAPTIVE_ERROR_COEFF: AdaptiveErrorCoeff = AdaptiveErrorCoeff::new(0.5, 0.25);

pub struct Rule<const NDIM: usize, const FINAL: usize, const TOTAL: usize> {
    initial_data: [Data<NDIM>; 3],
    final_data: [Data<NDIM>; FINAL],
    scales: Scales<TOTAL>,
    norms: Norms<TOTAL>,
    basic_error_coeff: BasicErrorCoeff,
    adaptive_error_coeff: AdaptiveErrorCoeff,
    evaluations: usize,
    ratio: f64,
}

impl<const NDIM: usize, const FINAL: usize, const TOTAL: usize> Rule<NDIM, FINAL, TOTAL> {
    pub(crate) const fn total() -> usize {
        TOTAL
    }

    pub(crate) const fn initial_data(&self) -> &[Data<NDIM>; 3] {
        &self.initial_data
    }

    pub(crate) const fn final_data(&self) -> &[Data<NDIM>; FINAL] {
        &self.final_data
    }

    pub(crate) const fn scales(&self) -> &Scales<TOTAL> {
        &self.scales
    }

    pub(crate) const fn norms(&self) -> &Norms<TOTAL> {
        &self.norms
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

pub type Rule07<const NDIM: usize> = Rule<NDIM, 3, 6>;

impl<const NDIM: usize> Rule07<NDIM> {
    /// Generate a fully-symmetric 7-point integration rule.
    pub const fn fs07() -> Self {
        fully_symmetric_07::generate_rule::<NDIM>()
    }
}

pub type Rule13 = Rule<2, 11, 14>;

impl Rule13 {
    /// Generate a fully-symmetric 13-point integration rule for NDIM = 2 dimensional integration.
    pub const fn fs13() -> Self {
        fully_symmetric_13_2d::generate_rule()
    }
}

pub type Rule11 = Rule<3, 10, 13>;

impl Rule11 {
    /// Generate a fully-symmetric 11-point integration rule for NDIM = 3 dimensional integration.
    pub const fn fs11() -> Self {
        fully_symmetric_11_3d::generate_rule()
    }
}

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

    pub(crate) const fn generator_inner(&self) -> &[f64; NDIM] {
        self.generator.generator()
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

pub(crate) struct Scales<const TOTAL: usize>(pub(crate) [[f64; TOTAL]; 3]);
pub(crate) struct Norms<const TOTAL: usize>(pub(crate) [[f64; TOTAL]; 3]);

pub(crate) const fn scales_norms<const NDIM: usize, const TOTAL: usize>(
    weights: &[[f64; 5]; TOTAL],
    rule_points: [f64; TOTAL],
) -> (Scales<TOTAL>, Norms<TOTAL>) {
    let two_ndim = two_pow_n_f64(NDIM);

    let mut scales = [[0f64; 3]; TOTAL];
    let mut norms = [[0f64; 3]; TOTAL];
    let mut we = [0f64; 14];

    let mut i = 0;
    while i < TOTAL {
        let mut k = 0;
        while k < 3 {
            scales[i][k] = if weights[i][k + 1] != 0.0 {
                -weights[i][k + 2] / weights[i][k + 1]
            } else {
                100.0
            };
            let mut j = 0;
            while j < TOTAL {
                we[j] = weights[j][k + 2] + scales[i][k] * weights[j][k + 1];
                j += 1;
            }
            norms[i][k] = 0.0;
            let mut j = 0;
            while j < TOTAL {
                // FIXME TODO remove this check when abs() is const in stable
                let weabs = if we[j].is_sign_negative() {
                    -we[j]
                } else {
                    we[j]
                };

                norms[i][k] += rule_points[j] * weabs;
                j += 1;
            }
            norms[i][k] = two_ndim / norms[i][k];

            k += 1;
        }
        i += 1;
    }

    let mut scales_swap = [[0f64; TOTAL]; 3];
    let mut norms_swap = [[0f64; TOTAL]; 3];

    let mut i = 0;
    while i < TOTAL {
        let mut j = 0;
        while j < 3 {
            scales_swap[j][i] = scales[i][j];
            norms_swap[j][i] = norms[i][j];
            j += 1;
        }
        i += 1;
    }

    (Scales(scales_swap), Norms(norms_swap))
}

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
