mod fully_symmetric_07;

use crate::multi::generator::Generator;
use crate::multi::two_pow_n_f64;

const LAMBDA0: f64 = 6.860_757_975_617_562_914_002_852e-1;
const LAMBDA1: f64 = 9.559_073_158_045_390_123_857_208e-1;
const LAMBDA2: f64 = 4.060_571_747_382_397_355_996_069e-1;
const LAMBDAP: f64 = 7.500_000_000_000_000_000_000_000e-1;

const LAMBDA0_SQ: f64 = LAMBDA0 * LAMBDA0;
const LAMBDA1_SQ: f64 = LAMBDA1 * LAMBDA1;
const LAMBDA2_SQ: f64 = LAMBDA2 * LAMBDA2;
const LAMBDAP_SQ: f64 = LAMBDAP * LAMBDAP;

pub struct Rule<const NDIM: usize, const WEIGHTS_LENGTH: usize> {
    weights: Weights<WEIGHTS_LENGTH>,
    scales: Scales<WEIGHTS_LENGTH>,
    norms: Norms<WEIGHTS_LENGTH>,
    generators: [Generator<NDIM>; WEIGHTS_LENGTH],
    rule_points: [f64; WEIGHTS_LENGTH],
    evaluations: usize,
    error_coefficients: [f64; 6],
    ratio: f64,
}

impl<const NDIM: usize, const WEIGHTS_LENGTH: usize> Rule<NDIM, WEIGHTS_LENGTH> {
    pub(crate) const fn weights(&self) -> &Weights<WEIGHTS_LENGTH> {
        &self.weights
    }

    pub(crate) const fn scales(&self) -> &Scales<WEIGHTS_LENGTH> {
        &self.scales
    }

    pub(crate) const fn norms(&self) -> &Norms<WEIGHTS_LENGTH> {
        &self.norms
    }

    pub(crate) const fn generators(&self) -> &[Generator<NDIM>; WEIGHTS_LENGTH] {
        &self.generators
    }

    pub(crate) const fn rule_points(&self) -> &[f64; WEIGHTS_LENGTH] {
        &self.rule_points
    }

    pub(crate) const fn evaluations(&self) -> usize {
        self.evaluations
    }

    pub(crate) const fn error_coefficients(&self) -> &[f64; 6] {
        &self.error_coefficients
    }

    pub(crate) const fn ratio(&self) -> f64 {
        self.ratio
    }
}

impl<const NDIM: usize> Rule<NDIM, { fully_symmetric_07::WEIGHTS_LENGTH }> {
    /// Generate a fully-symmetric 7-point integration rule.
    pub const fn fs07() -> Self {
        let weights = fully_symmetric_07::weights::<NDIM>();
        let rule_points = fully_symmetric_07::rule_points::<NDIM>();
        let (scales, norms) = weights.scales_norms::<NDIM>(&rule_points);
        let generators = fully_symmetric_07::generators::<NDIM>();
        let evaluations = fully_symmetric_07::evaluations::<NDIM>();
        let error_coefficients = fully_symmetric_07::error_coefficients();
        let ratio = (LAMBDA2 / LAMBDA1) * (LAMBDA2 / LAMBDA1);

        Self {
            weights,
            scales,
            norms,
            generators,
            rule_points,
            evaluations,
            error_coefficients,
            ratio,
        }
    }
}

pub(crate) struct Weights<const WEIGHTS_LENGTH: usize>([[f64; 5]; WEIGHTS_LENGTH]);
pub(crate) struct Scales<const WEIGHTS_LENGTH: usize>([[f64; 3]; WEIGHTS_LENGTH]);
pub(crate) struct Norms<const WEIGHTS_LENGTH: usize>([[f64; 3]; WEIGHTS_LENGTH]);

impl<const WEIGHTS_LENGTH: usize> Weights<WEIGHTS_LENGTH> {
    const fn scales_norms<const NDIM: usize>(
        &self,
        rule_points: &[f64; WEIGHTS_LENGTH],
    ) -> (Scales<WEIGHTS_LENGTH>, Norms<WEIGHTS_LENGTH>) {
        let two_ndim = two_pow_n_f64(NDIM);
        let weights = self.0;

        let mut scales = [[0f64; 3]; WEIGHTS_LENGTH];
        let mut norms = [[0f64; 3]; WEIGHTS_LENGTH];
        let mut we = [0f64; 14];

        let mut i = 0;
        while i < WEIGHTS_LENGTH {
            let mut k = 0;
            while k < 3 {
                scales[i][k] = if weights[i][k + 1] != 0.0 {
                    -weights[i][k + 2] / weights[i][k + 1]
                } else {
                    100.0
                };
                let mut j = 0;
                while j < WEIGHTS_LENGTH {
                    we[j] = weights[j][k + 2] + scales[i][k] * weights[j][k + 1];
                    j += 1;
                }
                norms[i][k] = 0.0;
                let mut j = 0;
                while j < WEIGHTS_LENGTH {
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

        (Scales(scales), Norms(norms))
    }
}
