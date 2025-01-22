use crate::multi::generator::Generator;
use crate::multi::rule::{
    scales_norms, BasicErrorCoeff, Data, Norms, Rule, Scales, ADAPTIVE_ERROR_COEFF,
};
use crate::multi::{two_pow_n, two_pow_n_f64};
use crate::{InitialisationError, InitialisationErrorKind};

const LAMBDA0: f64 = 6.860_757_975_617_562_914_002_852e-1;
const LAMBDA1: f64 = 9.559_073_158_045_390_123_857_208e-1;
const LAMBDA2: f64 = 4.060_571_747_382_397_355_996_069e-1;
const LAMBDA3: f64 = 8.952_547_092_523_562_576_415_410e-1;
const LAMBDAP: f64 = 2.500_000_000_000_000_000_000_000e-1;

const LAMBDA0_SQ: f64 = LAMBDA0 * LAMBDA0;
const LAMBDA1_SQ: f64 = LAMBDA1 * LAMBDA1;
const LAMBDA2_SQ: f64 = LAMBDA2 * LAMBDA2;
const LAMBDA3_SQ: f64 = LAMBDA3 * LAMBDA3;
const LAMBDAP_SQ: f64 = LAMBDAP * LAMBDAP;

const TOTAL: usize = 9;
const FINAL: usize = TOTAL - 3;
const RATIO: f64 = (LAMBDA2 / LAMBDA1) * (LAMBDA2 / LAMBDA1);

pub(crate) const fn generate_rule<const NDIM: usize>(
) -> Result<Rule<NDIM, FINAL, TOTAL>, InitialisationError> {
    if NDIM < 3 || NDIM > 15 {
        return Err(InitialisationError::new(
            InitialisationErrorKind::InvalidDimensionForRule09(NDIM),
        ));
    };

    let evaluations = evaluations::<NDIM>();
    let basic_error_coeff = BASIC_ERROR_COEFF;
    let adaptive_error_coeff = ADAPTIVE_ERROR_COEFF;
    let ratio = RATIO;

    let generators = generators::<NDIM>();
    let (weights, scales, norms) = weights_scales_norms::<NDIM>();

    let [g0, g1, g2, g3, g4, g5, g6, g7, g8] = generators;
    let [w0, w1, w2, w3, w4, w5, w6, w7, w8] = weights;

    let data0 = Data::new(g0, w0);
    let data1 = Data::new(g1, w1);
    let data2 = Data::new(g2, w2);
    let data3 = Data::new(g3, w3);
    let data4 = Data::new(g4, w4);
    let data5 = Data::new(g5, w5);
    let data6 = Data::new(g6, w6);
    let data7 = Data::new(g7, w7);
    let data8 = Data::new(g8, w8);

    let initial_data = [data0, data1, data2];
    let final_data = [data3, data4, data5, data6, data7, data8];

    Ok(Rule {
        initial_data,
        final_data,
        scales,
        norms,
        basic_error_coeff,
        adaptive_error_coeff,
        evaluations,
        ratio,
    })
}

const fn evaluations<const NDIM: usize>() -> usize {
    1 + 8 * NDIM
        + 2 * NDIM * (NDIM - 1)
        + 4 * NDIM * (NDIM - 1)
        + 4 * NDIM * (NDIM - 1) * (NDIM - 2) / 3
        + two_pow_n(NDIM)
}

const fn rule_points<const NDIM: usize>() -> [f64; TOTAL] {
    let ndim = NDIM as f64;
    [
        1.0,
        2.0 * ndim,
        2.0 * ndim,
        2.0 * ndim,
        2.0 * ndim,
        2.0 * ndim * (ndim - 1.0),
        4.0 * ndim * (ndim - 1.0),
        4.0 * ndim * (ndim - 1.0) * (ndim - 2.0) / 3.0,
        two_pow_n_f64(NDIM),
    ]
}

const fn initial_weights_1<const NDIM: usize>() -> [f64; TOTAL] {
    let two_ndim = two_pow_n_f64(NDIM);
    let ndim = NDIM as f64;

    let w8 = 1.0
        / ((3.0 * LAMBDA0_SQ) * (3.0 * LAMBDA0_SQ) * (3.0 * LAMBDA0_SQ) * (3.0 * LAMBDA0_SQ))
        / two_ndim;

    let w7 = (1.0 - 1.0 / (3.0 * LAMBDA0_SQ))
        / ((6.0 * LAMBDA1_SQ) * (6.0 * LAMBDA1_SQ) * (6.0 * LAMBDA1_SQ));

    let w6 = (1.0 - 7.0 * (LAMBDA0_SQ + LAMBDA1_SQ) / 5.0 + 7.0 * LAMBDA0_SQ * LAMBDA1_SQ / 3.0)
        / (84.0 * LAMBDA1_SQ * LAMBDA2_SQ * (LAMBDA2_SQ - LAMBDA0_SQ) * (LAMBDA2_SQ - LAMBDA1_SQ));

    let w5 = (1.0 - 7.0 * (LAMBDA0_SQ + LAMBDA2_SQ) / 5.0 + 7.0 * LAMBDA0_SQ * LAMBDA2_SQ / 3.0)
        / (84.0 * LAMBDA1_SQ * LAMBDA1_SQ * (LAMBDA1_SQ - LAMBDA0_SQ) * (LAMBDA1_SQ - LAMBDA2_SQ))
        - w6 * LAMBDA2_SQ / LAMBDA1_SQ
        - 2.0 * (ndim - 2.0) * w7;

    let w4 = 0.0;

    let w3 = (1.0
        - 9.0
            * ((LAMBDA0_SQ + LAMBDA1_SQ + LAMBDA2_SQ) / 7.0
                - (LAMBDA0_SQ * LAMBDA1_SQ + LAMBDA0_SQ * LAMBDA2_SQ + LAMBDA1_SQ * LAMBDA2_SQ)
                    / 5.0)
        - 3.0 * LAMBDA0_SQ * LAMBDA1_SQ * LAMBDA2_SQ)
        / (18.0
            * LAMBDA3_SQ
            * (LAMBDA3_SQ - LAMBDA0_SQ)
            * (LAMBDA3_SQ - LAMBDA1_SQ)
            * (LAMBDA3_SQ - LAMBDA2_SQ));

    let w2 = (1.0
        - 9.0
            * ((LAMBDA0_SQ + LAMBDA1_SQ + LAMBDA3_SQ) / 7.0
                - (LAMBDA0_SQ * LAMBDA1_SQ + LAMBDA0_SQ * LAMBDA3_SQ + LAMBDA1_SQ * LAMBDA3_SQ)
                    / 5.0)
        - 3.0 * LAMBDA0_SQ * LAMBDA1_SQ * LAMBDA3_SQ)
        / (18.0
            * LAMBDA2_SQ
            * (LAMBDA2_SQ - LAMBDA0_SQ)
            * (LAMBDA2_SQ - LAMBDA1_SQ)
            * (LAMBDA2_SQ - LAMBDA3_SQ))
        - 2.0 * (ndim - 1.0) * w6;

    let w1 = (1.0
        - 9.0
            * ((LAMBDA0_SQ + LAMBDA2_SQ + LAMBDA3_SQ) / 7.0
                - (LAMBDA0_SQ * LAMBDA2_SQ + LAMBDA0_SQ * LAMBDA3_SQ + LAMBDA2_SQ * LAMBDA3_SQ)
                    / 5.0)
        - 3.0 * LAMBDA0_SQ * LAMBDA2_SQ * LAMBDA3_SQ)
        / (18.0
            * LAMBDA1_SQ
            * (LAMBDA1_SQ - LAMBDA0_SQ)
            * (LAMBDA1_SQ - LAMBDA2_SQ)
            * (LAMBDA1_SQ - LAMBDA3_SQ))
        - 2.0 * (ndim - 1.0) * (w6 + w5 + (ndim - 2.0) * w7);

    let w0 = 0.0;

    [w0, w1, w2, w3, w4, w5, w6, w7, w8]
}

const fn initial_weights_2<const NDIM: usize>() -> [f64; TOTAL] {
    let two_ndim = two_pow_n_f64(NDIM);
    let ndim = NDIM as f64;

    let w8 = 1.0 / (108.0 * LAMBDA0_SQ * LAMBDA0_SQ * LAMBDA0_SQ * LAMBDA0_SQ) / two_ndim;

    let w7 = (1.0 - 27.0 * two_ndim * w8 * LAMBDA0_SQ * LAMBDA0_SQ * LAMBDA0_SQ)
        / ((6.0 * LAMBDA1_SQ) * (6.0 * LAMBDA1_SQ) * (6.0 * LAMBDA1_SQ));

    let w6 = (1.0
        - 5.0 * LAMBDA1_SQ / 3.0
        - 15.0 * two_ndim * w8 * LAMBDA0_SQ * LAMBDA0_SQ * (LAMBDA0_SQ - LAMBDA1_SQ))
        / (60.0 * LAMBDA1_SQ * LAMBDA2_SQ * (LAMBDA2_SQ - LAMBDA1_SQ));

    let w5 = (1.0
        - 9.0 * (8.0 * LAMBDA1_SQ * LAMBDA2_SQ * w6 + two_ndim * w8 * LAMBDA0_SQ * LAMBDA0_SQ))
        / (36.0 * LAMBDA1_SQ * LAMBDA1_SQ)
        - 2.0 * w7 * (ndim - 2.0);

    let w4 = 0.0;

    let w3 = (1.0
        - 7.0
            * ((LAMBDA1_SQ + LAMBDA2_SQ) / 5.0 - LAMBDA1_SQ * LAMBDA2_SQ / 3.0
                + two_ndim
                    * w8
                    * LAMBDA0_SQ
                    * (LAMBDA0_SQ - LAMBDA1_SQ)
                    * (LAMBDA0_SQ - LAMBDA2_SQ)))
        / (14.0 * LAMBDA3_SQ * (LAMBDA3_SQ - LAMBDA1_SQ) * (LAMBDA3_SQ - LAMBDA2_SQ));

    let w2 = (1.0
        - 7.0
            * ((LAMBDA1_SQ + LAMBDA3_SQ) / 5.0 - LAMBDA1_SQ * LAMBDA3_SQ / 3.0
                + two_ndim
                    * w8
                    * LAMBDA0_SQ
                    * (LAMBDA0_SQ - LAMBDA1_SQ)
                    * (LAMBDA0_SQ - LAMBDA3_SQ)))
        / (14.0 * LAMBDA2_SQ * (LAMBDA2_SQ - LAMBDA1_SQ) * (LAMBDA2_SQ - LAMBDA3_SQ))
        - 2.0 * (ndim - 1.0) * w6;

    let w1 = (1.0
        - 7.0
            * ((LAMBDA2_SQ + LAMBDA3_SQ) / 5.0 - LAMBDA2_SQ * LAMBDA3_SQ / 3.0
                + two_ndim
                    * w8
                    * LAMBDA0_SQ
                    * (LAMBDA0_SQ - LAMBDA2_SQ)
                    * (LAMBDA0_SQ - LAMBDA3_SQ)))
        / (14.0 * LAMBDA1_SQ * (LAMBDA1_SQ - LAMBDA2_SQ) * (LAMBDA1_SQ - LAMBDA3_SQ))
        - 2.0 * (ndim - 1.0) * (w6 + w5 + (ndim - 2.0) * w7);

    let w0 = 0.0;

    [w0, w1, w2, w3, w4, w5, w6, w7, w8]
}

const fn initial_weights_3<const NDIM: usize>() -> [f64; TOTAL] {
    let two_ndim = two_pow_n_f64(NDIM);
    let ndim = NDIM as f64;

    let w8 = 5.0 / (324.0 * LAMBDA0_SQ * LAMBDA0_SQ * LAMBDA0_SQ * LAMBDA0_SQ) / two_ndim;

    let w7 = (1.0 - 27.0 * two_ndim * w8 * LAMBDA0_SQ * LAMBDA0_SQ * LAMBDA0_SQ)
        / ((6.0 * LAMBDA1_SQ) * (6.0 * LAMBDA1_SQ) * (6.0 * LAMBDA1_SQ));

    let w6 = (1.0
        - 5.0 * LAMBDA1_SQ / 3.0
        - 15.0 * two_ndim * w8 * LAMBDA0_SQ * LAMBDA0_SQ * (LAMBDA0_SQ - LAMBDA1_SQ))
        / (60.0 * LAMBDA1_SQ * LAMBDA2_SQ * (LAMBDA2_SQ - LAMBDA1_SQ));

    let w5 = (1.0
        - 9.0 * (8.0 * LAMBDA1_SQ * LAMBDA2_SQ * w6 + two_ndim * w8 * LAMBDA0_SQ * LAMBDA0_SQ))
        / (36.0 * LAMBDA1_SQ * LAMBDA1_SQ)
        - 2.0 * w7 * (ndim - 2.0);

    let w4 = (1.0
        - 7.0
            * ((LAMBDA1_SQ + LAMBDA2_SQ) / 5.0 - LAMBDA1_SQ * LAMBDA2_SQ / 3.0
                + two_ndim
                    * w8
                    * LAMBDA0_SQ
                    * (LAMBDA0_SQ - LAMBDA1_SQ)
                    * (LAMBDA0_SQ - LAMBDA2_SQ)))
        / (14.0 * LAMBDAP_SQ * (LAMBDAP_SQ - LAMBDA1_SQ) * (LAMBDAP_SQ - LAMBDA2_SQ));

    let w3 = 0.0;

    let w2 = (1.0
        - 7.0
            * ((LAMBDA1_SQ + LAMBDAP_SQ) / 5.0 - LAMBDA1_SQ * LAMBDAP_SQ / 3.0
                + two_ndim
                    * w8
                    * LAMBDA0_SQ
                    * (LAMBDA0_SQ - LAMBDA1_SQ)
                    * (LAMBDA0_SQ - LAMBDAP_SQ)))
        / (14.0 * LAMBDA2_SQ * (LAMBDA2_SQ - LAMBDA1_SQ) * (LAMBDA2_SQ - LAMBDAP_SQ))
        - 2.0 * (ndim - 1.0) * w6;

    let w1 = (1.0
        - 7.0
            * ((LAMBDA2_SQ + LAMBDAP_SQ) / 5.0 - LAMBDA2_SQ * LAMBDAP_SQ / 3.0
                + two_ndim
                    * w8
                    * LAMBDA0_SQ
                    * (LAMBDA0_SQ - LAMBDA2_SQ)
                    * (LAMBDA0_SQ - LAMBDAP_SQ)))
        / (14.0 * LAMBDA1_SQ * (LAMBDA1_SQ - LAMBDA2_SQ) * (LAMBDA1_SQ - LAMBDAP_SQ))
        - 2.0 * (ndim - 1.0) * (w6 + w5 + (ndim - 2.0) * w7);

    let w0 = 0.0;

    [w0, w1, w2, w3, w4, w5, w6, w7, w8]
}

const fn initial_weights_4<const NDIM: usize>() -> [f64; TOTAL] {
    let two_ndim = two_pow_n_f64(NDIM);
    let ndim = NDIM as f64;

    let w8 = 2.0 / (81.0 * LAMBDA0_SQ * LAMBDA0_SQ * LAMBDA0_SQ * LAMBDA0_SQ) / two_ndim;

    let w7 = (2.0 - 27.0 * two_ndim * w8 * LAMBDA0_SQ * LAMBDA0_SQ * LAMBDA0_SQ)
        / ((6.0 * LAMBDA1_SQ) * (6.0 * LAMBDA1_SQ) * (6.0 * LAMBDA1_SQ));

    let w6 = (2.0
        - 15.0 * LAMBDA1_SQ / 9.0
        - 15.0 * two_ndim * w8 * LAMBDA0_SQ * (LAMBDA0_SQ - LAMBDA1_SQ))
        / (60.0 * LAMBDA1_SQ * LAMBDA2_SQ * (LAMBDA2_SQ - LAMBDA1_SQ));

    let w5 = (1.0
        - 9.0 * (8.0 * LAMBDA1_SQ * LAMBDA2_SQ * w6 + two_ndim * w8 * LAMBDA0_SQ * LAMBDA0_SQ))
        / (36.0 * LAMBDA1_SQ * LAMBDA1_SQ)
        - 2.0 * w7 * (ndim - 2.0);

    let w4 = 0.0;

    let w3 = (2.0
        - 7.0
            * ((LAMBDA1_SQ + LAMBDA2_SQ) / 5.0 - LAMBDA1_SQ * LAMBDA2_SQ / 3.0
                + two_ndim
                    * w8
                    * LAMBDA0_SQ
                    * (LAMBDA0_SQ - LAMBDA1_SQ)
                    * (LAMBDA0_SQ - LAMBDA2_SQ)))
        / (14.0 * LAMBDA3_SQ * (LAMBDA3_SQ - LAMBDA1_SQ) * (LAMBDA3_SQ - LAMBDA2_SQ));

    let w2 = (2.0
        - 7.0
            * ((LAMBDA1_SQ + LAMBDA3_SQ) / 5.0 - LAMBDA1_SQ * LAMBDA3_SQ / 3.0
                + two_ndim
                    * w8
                    * LAMBDA0_SQ
                    * (LAMBDA0_SQ - LAMBDA1_SQ)
                    * (LAMBDA0_SQ - LAMBDA3_SQ)))
        / (14.0 * LAMBDA2_SQ * (LAMBDA2_SQ - LAMBDA1_SQ) * (LAMBDA2_SQ - LAMBDA3_SQ))
        - 2.0 * (ndim - 1.0) * w6;

    let w1 = (2.0
        - 7.0
            * ((LAMBDA2_SQ + LAMBDA3_SQ) / 5.0 - LAMBDA2_SQ * LAMBDA3_SQ / 3.0
                + two_ndim
                    * w8
                    * LAMBDA0_SQ
                    * (LAMBDA0_SQ - LAMBDA2_SQ)
                    * (LAMBDA0_SQ - LAMBDA3_SQ)))
        / (14.0 * LAMBDA1_SQ * (LAMBDA1_SQ - LAMBDA2_SQ) * (LAMBDA1_SQ - LAMBDA3_SQ))
        - 2.0 * (ndim - 1.0) * (w6 + w5 + (ndim - 2.0) * w7);

    let w0 = 0.0;

    [w0, w1, w2, w3, w4, w5, w6, w7, w8]
}

const fn initial_weights_5<const NDIM: usize>() -> [f64; TOTAL] {
    let two_ndim = two_pow_n_f64(NDIM);
    let ndim = NDIM as f64;

    let mut weights = [0.0; TOTAL];
    weights[1] = 1.0 / (6.0 * LAMBDA1_SQ);

    weights
}

const fn initial_weights<const NDIM: usize>() -> [[f64; 5]; TOTAL] {
    let initial = [
        initial_weights_1::<NDIM>(),
        initial_weights_2::<NDIM>(),
        initial_weights_3::<NDIM>(),
        initial_weights_4::<NDIM>(),
        initial_weights_5::<NDIM>(),
    ];

    let mut output = [[0.0; 5]; TOTAL];

    let mut j = 0;
    while j < 5 {
        let mut i = 0;
        while i < TOTAL {
            output[i][j] = initial[j][i];
            i += 1;
        }
        j += 1;
    }

    output
}

const fn weights_scales_norms<const NDIM: usize>(
) -> ([[f64; 5]; TOTAL], Scales<TOTAL>, Norms<TOTAL>) {
    let mut weights = initial_weights::<NDIM>();
    let rule_points = rule_points::<NDIM>();

    let ndim = NDIM as f64;
    let two_ndim = two_pow_n_f64(NDIM);

    weights[0][0] = two_ndim;

    let mut i = 1;
    while i < TOTAL {
        let mut j = 1;
        while j < 5 {
            weights[i][j] -= weights[i][0];
            weights[0][j] -= rule_points[i] * weights[i][j];
            j += 1;
        }
        i += 1;
    }

    let mut k = 1;
    while k < TOTAL {
        weights[k][0] *= two_ndim;
        weights[0][0] -= rule_points[k] * weights[k][0];
        k += 1;
    }

    let (scales, norms) = scales_norms::<NDIM, TOTAL>(&weights, rule_points);
    (weights, scales, norms)
}

const fn generators<const NDIM: usize>() -> [Generator<NDIM>; TOTAL] {
    let gen_0 = Generator::new([0.0; NDIM]);
    let gen_1 = {
        let mut init = [0.0; NDIM];
        init[0] = LAMBDA1;
        Generator::new(init)
    };
    let gen_2 = {
        let mut init = [0.0; NDIM];
        init[0] = LAMBDA2;
        Generator::new(init)
    };
    let gen_3 = {
        let mut init = [0.0; NDIM];
        init[0] = LAMBDA3;
        Generator::new(init)
    };
    let gen_4 = {
        let mut init = [0.0; NDIM];
        init[0] = LAMBDAP;
        Generator::new(init)
    };
    let gen_5 = {
        let mut init = [0.0; NDIM];
        init[0] = LAMBDA1;
        init[1] = init[0];
        Generator::new(init)
    };
    let gen_6 = {
        let mut init = [0.0; NDIM];
        init[0] = LAMBDA1;
        init[1] = LAMBDA2;
        Generator::new(init)
    };
    let gen_7 = {
        let mut init = [0.0; NDIM];
        init[0] = LAMBDA1;
        init[1] = init[0];
        init[2] = init[0];
        Generator::new(init)
    };
    let gen_8 = Generator::new([LAMBDA0; NDIM]);

    [
        gen_0, gen_1, gen_2, gen_3, gen_4, gen_5, gen_6, gen_7, gen_8,
    ]
}

const BASIC_ERROR_COEFF: BasicErrorCoeff = BasicErrorCoeff::new(5.0, 5.0, 1.0, 5.0);

#[cfg(test)]
mod wtest {
    use super::*;
}

#[cfg(test)]
fn rel_or_abs_diff(a: f64, b: f64) -> f64 {
    if a == 0.0 {
        (a - b).abs()
    } else {
        (a - b).abs() / a.abs()
    }
}

#[cfg(test)]
fn assert_check_vec_tol<const WL: usize, const TY: usize>(
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

#[cfg(test)]
fn assert_check_vec_data_tol<const WL: usize, const TY: usize>(
    calc: &[Data<TY>; WL],
    should_be: &[Data<TY>; WL],
    tol: f64,
) {
    for (x, y) in calc.iter().zip(should_be.iter()) {
        let genx = x.generator_inner();
        let geny = y.generator_inner();
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn check_evaluations_correct() {
        assert_eq!(33, evaluations::<2>());
        assert_eq!(77, evaluations::<3>());
        assert_eq!(153, evaluations::<4>());
        assert_eq!(273, evaluations::<5>());
        assert_eq!(453, evaluations::<6>());
        assert_eq!(717, evaluations::<7>());
        assert_eq!(1105, evaluations::<8>());
        assert_eq!(1689, evaluations::<9>());
        assert_eq!(2605, evaluations::<10>());
    }

    #[test]
    fn check_rule_points_correct() {
        {
            const NDIM: usize = 3;
            let rule_points = rule_points::<NDIM>();
            let should_be = [1.0, 6.0, 6.0, 6.0, 6.0, 12.0, 24.0, 8.0, 8.0];
            let mut sum: usize = 0;
            for (x, y) in rule_points.iter().zip(should_be.iter()) {
                assert_eq!(x, y);
                sum += *x as usize;
            }
            assert_eq!(sum, evaluations::<NDIM>());
        }

        {
            const NDIM: usize = 4;
            let rule_points = rule_points::<NDIM>();
            let should_be = [1.0, 8.0, 8.0, 8.0, 8.0, 24.0, 48.0, 32.0, 16.0];
            let mut sum: usize = 0;
            for (x, y) in rule_points.iter().zip(should_be.iter()) {
                assert_eq!(x, y);
                sum += *x as usize;
            }
            assert_eq!(sum, evaluations::<NDIM>());
        }

        {
            const NDIM: usize = 5;
            let rule_points = rule_points::<NDIM>();
            let should_be = [1.0, 10.0, 10.0, 10.0, 10.0, 40.0, 80.0, 80.0, 32.0];
            let mut sum: usize = 0;
            for (x, y) in rule_points.iter().zip(should_be.iter()) {
                assert_eq!(x, y);
                sum += *x as usize;
            }
            assert_eq!(sum, evaluations::<NDIM>());
        }

        {
            const NDIM: usize = 6;
            let rule_points = rule_points::<NDIM>();
            let should_be = [1.0, 12.0, 12.0, 12.0, 12.0, 60.0, 120.0, 160.0, 64.0];
            let mut sum: usize = 0;
            for (x, y) in rule_points.iter().zip(should_be.iter()) {
                assert_eq!(x, y);
                sum += *x as usize;
            }
            assert_eq!(sum, evaluations::<NDIM>());
        }

        {
            const NDIM: usize = 15;
            let rule_points = rule_points::<NDIM>();
            let should_be = [1.0, 30.0, 30.0, 30.0, 30.0, 420.0, 840.0, 3640.0, 32768.0];
            let mut sum: usize = 0;
            for (x, y) in rule_points.iter().zip(should_be.iter()) {
                assert_eq!(x, y);
                sum += *x as usize;
            }
            assert_eq!(sum, evaluations::<NDIM>());
        }
    }

    #[test]
    fn test_initial_weights_1() {
        {
            const NDIM: usize = 3;
            let w = initial_weights_1::<NDIM>();
            let should_be = [
                0.0000000000000000,
                -8.0377377870386235E-002,
                7.8899967360456499E-002,
                5.7693384490972686E-002,
                0.0000000000000000,
                4.9071479215722618E-003,
                2.2543144647178933E-002,
                1.7708782258391350E-003,
                3.1437514369143479E-002,
            ];
            let tol = 1e-10;
            assert_check_vec_tol(&[w; 1], &[should_be; 1], tol);
        }

        {
            const NDIM: usize = 9;
            let w = initial_weights_1::<NDIM>();
            let should_be = [
                0.0000000000000000,
                -0.26102711772491327,
                -0.19161776840569070,
                5.7693384490972686E-002,
                0.0000000000000000,
                -1.6343390788497357E-002,
                2.2543144647178933E-002,
                1.7708782258391350E-003,
                4.9121116201786687E-004,
            ];
            let tol = 1e-10;
            assert_check_vec_tol(&[w; 1], &[should_be; 1], tol);
        }

        {
            const NDIM: usize = 14;
            let w = initial_weights_1::<NDIM>();
            let should_be = [
                0.0000000000000000,
                -0.21677196276138094,
                -0.41704921487748003,
                5.7693384490972686E-002,
                0.0000000000000000,
                -3.4052173046888706E-002,
                2.2543144647178933E-002,
                1.7708782258391350E-003,
                1.5350348813058340E-005,
            ];
            let tol = 1e-10;
            assert_check_vec_tol(&[w; 1], &[should_be; 1], tol);
        }
    }

    #[test]
    fn test_initial_weights_2() {
        {
            const NDIM: usize = 3;
            let w = initial_weights_2::<NDIM>();
            let should_be = [
                0.0000000000000000,
                -0.14720886256392740,
                5.1786394829528515E-002,
                9.2693011093116021E-002,
                0.0000000000000000,
                1.9941502172517367E-003,
                3.6218917910451759E-002,
                2.8451795033440693E-003,
                2.3578135776857620E-002,
            ];
            let tol = 1e-10;
            assert_check_vec_tol(&[w; 1], &[should_be; 1], tol);
        }

        {
            const NDIM: usize = 9;
            let w = initial_weights_2::<NDIM>();
            let should_be = [
                0.0000000000000000,
                -0.36677060181546750,
                -0.38284062009589259,
                9.2693011093116021E-002,
                0.0000000000000000,
                -3.2148003822877094E-002,
                3.6218917910451759E-002,
                2.8451795033440693E-003,
                3.6840837151340031E-004,
            ];
            let tol = 1e-10;
            assert_check_vec_tol(&[w; 1], &[should_be; 1], tol);
        }

        {
            const NDIM: usize = 14;
            let w = initial_weights_2::<NDIM>();
            let should_be = [
                0.0000000000000000,
                -0.23676897249056994,
                -0.74502979920041013,
                9.2693011093116021E-002,
                0.0000000000000000,
                -6.0599798856317791E-002,
                3.6218917910451759E-002,
                2.8451795033440693E-003,
                1.1512761609793760E-005,
            ];
            let tol = 1e-10;
            assert_check_vec_tol(&[w; 1], &[should_be; 1], tol);
        }
    }

    #[test]
    fn test_initial_weights_3() {
        {
            const NDIM: usize = 3;
            let w = initial_weights_3::<NDIM>();
            let should_be = [
                0.0000000000000000,
                1.1432093192986353E-003,
                0.22538950191708801,
                0.0000000000000000,
                -0.23866687325750005,
                7.8201456258927757E-003,
                8.8673713839060881E-003,
                6.9657694833419933E-004,
                3.9296892961429367E-002,
            ];
            let tol = 1e-10;
            assert_check_vec_tol(&[w; 1], &[should_be; 1], tol);
        }

        {
            const NDIM: usize = 9;
            let w = initial_weights_3::<NDIM>();
            let should_be = [
                0.0000000000000000,
                -0.14059453113821502,
                0.11898104531021494,
                0.0000000000000000,
                -0.23866687325750005,
                -5.3877775411761628E-004,
                8.8673713839060881E-003,
                6.9657694833419933E-004,
                6.1401395252233385E-004,
            ];
            let tol = 1e-10;
            assert_check_vec_tol(&[w; 1], &[should_be; 1], tol);
        }

        {
            const NDIM: usize = 14;
            let w = initial_weights_3::<NDIM>();
            let should_be = [
                0.0000000000000000,
                -0.18208585053604776,
                3.0307331471154053E-002,
                0.0000000000000000,
                -0.23866687325750005,
                -7.5045472374596096E-003,
                8.8673713839060881E-003,
                6.9657694833419933E-004,
                1.9187936016322933E-005,
            ];
            let tol = 1e-10;
            assert_check_vec_tol(&[w; 1], &[should_be; 1], tol);
        }
    }

    #[test]
    fn test_initial_weights_4() {
        {
            const NDIM: usize = 3;
            let w = initial_weights_4::<NDIM>();
            let should_be = [
                0.0000000000000000,
                1.8330509758613034,
                2.1790349474124113,
                -1.3291693874368409,
                0.0000000000000000,
                0.10213149571758670,
                -0.30290445231242258,
                3.5417564516782686E-003,
                6.2875028738286987E-002,
            ];
            let tol = 1e-10;
            assert_check_vec_tol(&[w; 1], &[should_be; 1], tol);
        }

        {
            const NDIM: usize = 9;
            let w = initial_weights_4::<NDIM>();
            let should_be = [
                0.0000000000000000,
                4.5398339969403088,
                5.8138883751614818,
                -1.3291693874368409,
                0.0000000000000000,
                5.9630418297447481E-002,
                -0.30290445231242258,
                3.5417564516782686E-003,
                9.8242232403573417E-004,
            ];
            let tol = 1e-10;
            assert_check_vec_tol(&[w; 1], &[should_be; 1], tol);
        }

        {
            const NDIM: usize = 14;
            let w = initial_weights_4::<NDIM>();
            let should_be = [
                0.0000000000000000,
                7.1850797241907554,
                8.8429328982857083,
                -1.3291693874368409,
                0.0000000000000000,
                2.4212853780664789E-002,
                -0.30290445231242258,
                3.5417564516782686E-003,
                3.0700697626116693E-005,
            ];
            let tol = 1e-10;
            assert_check_vec_tol(&[w; 1], &[should_be; 1], tol);
        }
    }

    #[test]
    fn test_initial_weights_5() {
        {
            const NDIM: usize = 3;
            let w = initial_weights_5::<NDIM>();
            let should_be = [
                0.0000000000000000,
                0.18239678493024578,
                0.0000000000000000,
                0.0000000000000000,
                0.0000000000000000,
                0.0000000000000000,
                0.0000000000000000,
                0.0000000000000000,
                0.0000000000000000,
            ];
            let tol = 1e-10;
            assert_check_vec_tol(&[w; 1], &[should_be; 1], tol);
        }

        {
            const NDIM: usize = 9;
            let w = initial_weights_5::<NDIM>();
            let should_be = [
                0.0000000000000000,
                0.18239678493024578,
                0.0000000000000000,
                0.0000000000000000,
                0.0000000000000000,
                0.0000000000000000,
                0.0000000000000000,
                0.0000000000000000,
                0.0000000000000000,
            ];
            let tol = 1e-10;
            assert_check_vec_tol(&[w; 1], &[should_be; 1], tol);
        }

        {
            const NDIM: usize = 14;
            let w = initial_weights_5::<NDIM>();
            let should_be = [
                0.0000000000000000,
                0.18239678493024578,
                0.0000000000000000,
                0.0000000000000000,
                0.0000000000000000,
                0.0000000000000000,
                0.0000000000000000,
                0.0000000000000000,
                0.0000000000000000,
            ];
            let tol = 1e-10;
            assert_check_vec_tol(&[w; 1], &[should_be; 1], tol);
        }
    }

    #[test]
    fn check_generators_correct() {
        let tol = 1e-14;

        {
            const NDIM: usize = 3;
            let generators = generators::<NDIM>();
            let should_be = [
                Generator::new([0.0000000000000000, 0.0000000000000000, 0.0000000000000000]),
                Generator::new([0.95590731580453892, 0.0000000000000000, 0.0000000000000000]),
                Generator::new([0.40605717473823955, 0.0000000000000000, 0.0000000000000000]),
                Generator::new([0.89525470925235517, 0.0000000000000000, 0.0000000000000000]),
                Generator::new([0.25000000000000000, 0.0000000000000000, 0.0000000000000000]),
                Generator::new([0.95590731580453892, 0.95590731580453892, 0.0000000000000000]),
                Generator::new([0.95590731580453892, 0.40605717473823955, 0.0000000000000000]),
                Generator::new([
                    0.95590731580453892,
                    0.95590731580453892,
                    0.95590731580453892,
                ]),
                Generator::new([
                    0.68607579756175630,
                    0.68607579756175630,
                    0.68607579756175630,
                ]),
            ];
            for (x, y) in generators.iter().zip(should_be.iter()) {
                let gx = x.generator();
                let gy = y.generator();
                for (a, b) in gx.iter().zip(gy.iter()) {
                    let val = rel_or_abs_diff(*a, *b);
                    assert!(val < tol);
                }
            }
        }

        {
            const NDIM: usize = 6;
            let generators = generators::<NDIM>();
            let should_be = [
                Generator::new([
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                ]),
                Generator::new([
                    0.95590731580453892,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                ]),
                Generator::new([
                    0.40605717473823955,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                ]),
                Generator::new([
                    0.89525470925235517,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                ]),
                Generator::new([
                    0.25000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                ]),
                Generator::new([
                    0.95590731580453892,
                    0.95590731580453892,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                ]),
                Generator::new([
                    0.95590731580453892,
                    0.40605717473823955,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                ]),
                Generator::new([
                    0.95590731580453892,
                    0.95590731580453892,
                    0.95590731580453892,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                ]),
                Generator::new([
                    0.68607579756175630,
                    0.68607579756175630,
                    0.68607579756175630,
                    0.68607579756175630,
                    0.68607579756175630,
                    0.68607579756175630,
                ]),
            ];
            for (x, y) in generators.iter().zip(should_be.iter()) {
                let gx = x.generator();
                let gy = y.generator();
                for (a, b) in gx.iter().zip(gy.iter()) {
                    let val = rel_or_abs_diff(*a, *b);
                    assert!(val < tol);
                }
            }
        }

        {
            const NDIM: usize = 15;
            let generators = generators::<NDIM>();
            let should_be = [
                Generator::new([
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                ]),
                Generator::new([
                    0.95590731580453892,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                ]),
                Generator::new([
                    0.40605717473823955,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                ]),
                Generator::new([
                    0.89525470925235517,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                ]),
                Generator::new([
                    0.25000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                ]),
                Generator::new([
                    0.95590731580453892,
                    0.95590731580453892,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                ]),
                Generator::new([
                    0.95590731580453892,
                    0.40605717473823955,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                ]),
                Generator::new([
                    0.95590731580453892,
                    0.95590731580453892,
                    0.95590731580453892,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                    0.0000000000000000,
                ]),
                Generator::new([
                    0.68607579756175630,
                    0.68607579756175630,
                    0.68607579756175630,
                    0.68607579756175630,
                    0.68607579756175630,
                    0.68607579756175630,
                    0.68607579756175630,
                    0.68607579756175630,
                    0.68607579756175630,
                    0.68607579756175630,
                    0.68607579756175630,
                    0.68607579756175630,
                    0.68607579756175630,
                    0.68607579756175630,
                    0.68607579756175630,
                ]),
            ];
            for (x, y) in generators.iter().zip(should_be.iter()) {
                let gx = x.generator();
                let gy = y.generator();
                for (a, b) in gx.iter().zip(gy.iter()) {
                    let val = rel_or_abs_diff(*a, *b);
                    assert!(val < tol);
                }
            }
        }
    }

    #[test]
    fn check_weights_scales_norms_correct() {
        {
            let tol = 1e-12;
            const NDIM: usize = 3;
            let (weights, scales, norms) = weights_scales_norms::<NDIM>();
            let weights_should_be = [
                [
                    -1.6230738498982409,
                    0.11469061638550071,
                    0.64908278336139280,
                    -9.3818203584165847,
                    0.10850352165580529,
                ],
                [
                    -0.64301902296308988,
                    -6.6831484693541165E-002,
                    8.1520587189684870E-002,
                    1.9134283537316896,
                    0.26277416280063204,
                ],
                [
                    0.63119973888365200,
                    -2.7113572530927985E-002,
                    0.14648953455663150,
                    2.1001349800519549,
                    -7.8899967360456499E-002,
                ],
                [
                    0.46154707592778149,
                    3.4999626602143334E-002,
                    -5.7693384490972686E-002,
                    -1.3868627719278135,
                    -5.7693384490972686E-002,
                ],
                [
                    0.0000000000000000,
                    0.0000000000000000,
                    -0.23866687325750005,
                    0.0000000000000000,
                    0.0000000000000000,
                ],
                [
                    3.9257183372578094E-002,
                    -2.9129977043205251E-003,
                    2.9129977043205139E-003,
                    9.7224347796014438E-002,
                    -4.9071479215722618E-003,
                ],
                [
                    0.18034515717743146,
                    1.3675773263272826E-002,
                    -1.3675773263272845E-002,
                    -0.32544759695960152,
                    -2.2543144647178933E-002,
                ],
                [
                    1.4167025806713080E-002,
                    1.0743012775049343E-003,
                    -1.0743012775049356E-003,
                    1.7708782258391337E-003,
                    -1.7708782258391350E-003,
                ],
                [
                    0.25150011495314784,
                    -7.8593785922858594E-003,
                    7.8593785922858872E-003,
                    3.1437514369143507E-002,
                    -3.1437514369143479E-002,
                ],
            ];
            let scales_should_be = Scales([
                [
                    -5.6594236199732419,
                    1.2197931493442238,
                    5.4028119824318024,
                    1.6483999999999890,
                    100.00000000000000,
                    0.99999999999999611,
                    1.0000000000000013,
                    1.0000000000000011,
                    1.0000000000000036,
                ],
                [
                    14.453965809771026,
                    -23.471719472278327,
                    -14.336416498340789,
                    -24.038506046474300,
                    0.0000000000000000,
                    -33.376046830319424,
                    -23.797381741741187,
                    1.6483999999999979,
                    -3.9999999999999947,
                ],
                [
                    1.1565295167740556E-002,
                    -0.13733159242056445,
                    3.7568998235772734E-002,
                    -4.1599922976356049E-002,
                    100.00000000000000,
                    5.0472417998297159E-002,
                    -6.9268124447012774E-002,
                    1.0000000000000007,
                    0.99999999999999911,
                ],
            ]);
            let norms_should_be = Norms([
                [
                    0.76781207907455518,
                    2.5914954991584516,
                    1.1299189578775166,
                    2.3990862545604470,
                    6.0499504513972208E-002,
                    2.5507539668184198,
                    2.5507539668184243,
                    2.5507539668184243,
                    2.5507539668184247,
                ],
                [
                    8.5675998466278228E-002,
                    0.11704928336754453,
                    0.15503283097331050,
                    0.11412437973700976,
                    0.15677827347450360,
                    7.3432577243211505E-002,
                    0.11550210593131525,
                    0.14328726963651109,
                    0.15740251232779079,
                ],
                [
                    2.3399733757982366,
                    1.4678090079104389,
                    1.9920730478289697,
                    2.5033785128094785,
                    1.5673043855019438E-003,
                    1.7251115407760926,
                    2.2598818600894051,
                    0.15213503508669102,
                    0.15213503508669127,
                ],
            ]);
            assert_check_vec_tol(&weights, &weights_should_be, tol);
            assert_check_vec_tol(&scales.0, &scales_should_be.0, tol);
            assert_check_vec_tol(&norms.0, &norms_should_be.0, tol);
        }

        {
            let tol = 1e-12;
            const NDIM: usize = 6;
            let (weights, scales, norms) = weights_scales_norms::<NDIM>();
            let weights_should_be = [
                [
                    33.151498239465923,
                    0.96894426965246350,
                    0.55860252984132153,
                    -38.571794294575213,
                    -1.7067535791546047,
                ],
                [
                    -12.964995575216268,
                    -0.10562490738713654,
                    0.12031400988328041,
                    3.3252689261333517,
                    0.38497484079299993,
                ],
                [
                    -3.6069696334474948,
                    -0.10916821211056493,
                    0.22854417413626860,
                    4.0528205618095638,
                    5.6358900522617106E-002,
                ],
                [
                    3.6923766074222519,
                    3.4999626602143334E-002,
                    -5.7693384490972686E-002,
                    -1.3868627719278135,
                    -5.7693384490972686E-002,
                ],
                [
                    0.0000000000000000,
                    0.0000000000000000,
                    -0.23866687325750005,
                    0.0000000000000000,
                    0.0000000000000000,
                ],
                [
                    -0.36595977174160310,
                    -9.3588053693501299E-003,
                    9.3588053693501282E-003,
                    8.6599078440979652E-002,
                    5.7181214334625485E-003,
                ],
                [
                    1.4427612574194517,
                    1.3675773263272826E-002,
                    -1.3675773263272845E-002,
                    -0.32544759695960152,
                    -2.2543144647178933E-002,
                ],
                [
                    0.11333620645370464,
                    1.0743012775049343E-003,
                    -1.0743012775049356E-003,
                    1.7708782258391337E-003,
                    -1.7708782258391350E-003,
                ],
                [
                    0.25150011495314784,
                    -9.8242232403573243E-004,
                    9.8242232403573590E-004,
                    3.9296892961429384E-003,
                    -3.9296892961429349E-003,
                ],
            ];
            let scales_should_be = Scales([
                [
                    -0.57650635576974774,
                    1.1390685479354350,
                    2.0935047823701680,
                    1.6483999999999890,
                    100.00000000000000,
                    0.99999999999999978,
                    1.0000000000000013,
                    1.0000000000000011,
                    1.0000000000000036,
                ],
                [
                    69.050518452775435,
                    -27.638252015366099,
                    -17.733204432474768,
                    -24.038506046474300,
                    0.0000000000000000,
                    -9.2532192970472096,
                    -23.797381741741187,
                    1.6483999999999979,
                    -3.9999999999999947,
                ],
                [
                    -4.4248747313127841E-002,
                    -0.11577254331746685,
                    -1.3906093216585227E-002,
                    -4.1599922976356049E-002,
                    100.00000000000000,
                    -6.6029818520062553E-002,
                    -6.9268124447012774E-002,
                    1.0000000000000007,
                    0.99999999999999911,
                ],
            ]);
            let norms_should_be = Norms([
                [
                    4.8072378429509444,
                    10.111369549905003,
                    6.7275829740942514,
                    8.1748547032549119,
                    0.10055134204812069,
                    10.203015867273686,
                    10.203015867273693,
                    10.203015867273693,
                    10.203015867273693,
                ],
                [
                    7.5088119053814448E-002,
                    0.34413984358525346,
                    0.46769948560737951,
                    0.40581058418476856,
                    0.33945822334345882,
                    0.40933571314201506,
                    0.40986275239428399,
                    0.31398875090344353,
                    0.36710460234560116,
                ],
                [
                    10.626688399072346,
                    5.5083347238790363,
                    7.3807000568771750,
                    10.367181928889480,
                    3.3927451889699933E-003,
                    9.9022119597539149,
                    9.7525811632034358,
                    0.32202194117196942,
                    0.32202194117196986,
                ],
            ]);
            assert_check_vec_tol(&weights, &weights_should_be, tol);
            assert_check_vec_tol(&scales.0, &scales_should_be.0, tol);
            assert_check_vec_tol(&norms.0, &norms_should_be.0, tol);
        }

        {
            let tol = 1e-12;
            const NDIM: usize = 15;
            let (weights, scales, norms) = weights_scales_norms::<NDIM>();

            let weights_should_be = [
                [
                    291277.61628182390,
                    6.0258428637611914,
                    -2.2069758650267364,
                    -250.32473522703492,
                    -13.360990958851707,
                ],
                [
                    -6116.8154397443632,
                    1.0043900473143297E-002,
                    4.6452020230009494E-003,
                    7.9433003401195883,
                    0.36906717798885058,
                ],
                [
                    -15143.256200702785,
                    -0.35533213084947590,
                    0.47470809287517979,
                    9.9108773070823926,
                    0.46213550417183791,
                ],
                [
                    1890.4968230001930,
                    3.4999626602143334E-002,
                    -5.7693384490972686E-002,
                    -1.3868627719278135,
                    -5.7693384490972686E-002,
                ],
                [
                    0.0000000000000000,
                    0.0000000000000000,
                    -0.23866687325750005,
                    0.0000000000000000,
                    0.0000000000000000,
                ],
                [
                    -1231.8778818090427,
                    -2.8696228364438955E-002,
                    2.8696228364438969E-002,
                    5.4723270375875238E-002,
                    3.7593929498566978E-002,
                ],
                [
                    738.69376379875928,
                    1.3675773263272826E-002,
                    -1.3675773263272845E-002,
                    -0.32544759695960152,
                    -2.2543144647178933E-002,
                ],
                [
                    58.028137704296775,
                    1.0743012775049343E-003,
                    -1.0743012775049356E-003,
                    1.7708782258391337E-003,
                    -1.7708782258391350E-003,
                ],
                [
                    0.25150011495314784,
                    -1.9187936016322899E-006,
                    1.9187936016322967E-006,
                    7.6751744065291766E-006,
                    -7.6751744065291698E-006,
                ],
            ];
            let scales_should_be = Scales([
                [
                    0.36625181155971820,
                    -0.46248984997630171,
                    1.3359560018969221,
                    1.6483999999999890,
                    100.00000000000000,
                    1.0000000000000004,
                    1.0000000000000013,
                    1.0000000000000011,
                    1.0000000000000036,
                ],
                [
                    -113.42431931125915,
                    -1710.0010507159734,
                    -20.877835149291140,
                    -24.038506046474300,
                    0.0000000000000000,
                    -1.9069847675065752,
                    -23.797381741741187,
                    1.6483999999999979,
                    -3.9999999999999947,
                ],
                [
                    -5.3374633340703630E-002,
                    -4.6462699657066495E-002,
                    -4.6629121706671939E-002,
                    -4.1599922976356049E-002,
                    100.00000000000000,
                    -0.68698250744787848,
                    -6.9268124447012774E-002,
                    1.0000000000000007,
                    0.99999999999999911,
                ],
            ]);
            let norms_should_be = Norms([
                [
                    897.02747938104835,
                    444.13291124505685,
                    1417.5127327513665,
                    893.04584515699275,
                    7.2549218010363896,
                    2089.5776496176491,
                    2089.5776496176504,
                    2089.5776496176500,
                    2089.5776496176441,
                ],
                [
                    6.1388350037895334,
                    0.36612995128823062,
                    34.641198384247041,
                    32.107709132362274,
                    28.982502601066287,
                    30.573407450159156,
                    32.451485859286883,
                    27.191716954556536,
                    30.981555069519679,
                ],
                [
                    1088.4500436838150,
                    1094.2096916225107,
                    1097.6097232600287,
                    936.90125600431941,
                    0.28965067384530918,
                    46.239680023357145,
                    847.72153102289974,
                    27.336981838259597,
                    27.336981838259643,
                ],
            ]);
            assert_check_vec_tol(&weights, &weights_should_be, tol);
            assert_check_vec_tol(&scales.0, &scales_should_be.0, tol);
            assert_check_vec_tol(&norms.0, &norms_should_be.0, tol);
        }
    }

    #[test]
    fn check_data09_correct() {
        {
            let tol = 1e-12;
            const NDIM: usize = 3;
            let rule = generate_rule::<NDIM>().unwrap();

            let initial_data = rule.initial_data();
            let final_data = rule.final_data();
            let initial_should_be = [
                Data::new(
                    Generator::new([0.0000000000000000, 0.0000000000000000, 0.0000000000000000]),
                    [
                        -1.6230738498982409,
                        0.11469061638550071,
                        0.64908278336139280,
                        -9.3818203584165847,
                        0.10850352165580529,
                    ],
                ),
                Data::new(
                    Generator::new([0.95590731580453892, 0.0000000000000000, 0.0000000000000000]),
                    [
                        -0.64301902296308988,
                        -6.6831484693541165E-002,
                        8.1520587189684870E-002,
                        1.9134283537316896,
                        0.26277416280063204,
                    ],
                ),
                Data::new(
                    Generator::new([0.40605717473823955, 0.0000000000000000, 0.0000000000000000]),
                    [
                        0.63119973888365200,
                        -2.7113572530927985E-002,
                        0.14648953455663150,
                        2.1001349800519549,
                        -7.8899967360456499E-002,
                    ],
                ),
            ];
            let final_should_be = [
                Data::new(
                    Generator::new([0.89525470925235517, 0.0000000000000000, 0.0000000000000000]),
                    [
                        0.46154707592778149,
                        3.4999626602143334E-002,
                        -5.7693384490972686E-002,
                        -1.3868627719278135,
                        -5.7693384490972686E-002,
                    ],
                ),
                Data::new(
                    Generator::new([0.25000000000000000, 0.0000000000000000, 0.0000000000000000]),
                    [
                        0.0000000000000000,
                        0.0000000000000000,
                        -0.23866687325750005,
                        0.0000000000000000,
                        0.0000000000000000,
                    ],
                ),
                Data::new(
                    Generator::new([0.95590731580453892, 0.95590731580453892, 0.0000000000000000]),
                    [
                        3.9257183372578094E-002,
                        -2.9129977043205251E-003,
                        2.9129977043205139E-003,
                        9.7224347796014438E-002,
                        -4.9071479215722618E-003,
                    ],
                ),
                Data::new(
                    Generator::new([0.95590731580453892, 0.40605717473823955, 0.0000000000000000]),
                    [
                        0.18034515717743146,
                        1.3675773263272826E-002,
                        -1.3675773263272845E-002,
                        -0.32544759695960152,
                        -2.2543144647178933E-002,
                    ],
                ),
                Data::new(
                    Generator::new([
                        0.95590731580453892,
                        0.95590731580453892,
                        0.95590731580453892,
                    ]),
                    [
                        1.4167025806713080E-002,
                        1.0743012775049343E-003,
                        -1.0743012775049356E-003,
                        1.7708782258391337E-003,
                        -1.7708782258391350E-003,
                    ],
                ),
                Data::new(
                    Generator::new([
                        0.68607579756175630,
                        0.68607579756175630,
                        0.68607579756175630,
                    ]),
                    [
                        0.25150011495314784,
                        -7.8593785922858594E-003,
                        7.8593785922858872E-003,
                        3.1437514369143507E-002,
                        -3.1437514369143479E-002,
                    ],
                ),
            ];
            assert_check_vec_data_tol(&initial_data, &initial_should_be, tol);
            assert_check_vec_data_tol(&final_data, &final_should_be, tol);
        }

        {
            let tol = 1e-12;
            const NDIM: usize = 6;
            let rule = generate_rule::<NDIM>().unwrap();

            let initial_data = rule.initial_data();
            let final_data = rule.final_data();

            let initial_should_be = [
                Data::new(
                    Generator::new([
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                    ]),
                    [
                        33.151498239465923,
                        0.96894426965246350,
                        0.55860252984132153,
                        -38.571794294575213,
                        -1.7067535791546047,
                    ],
                ),
                Data::new(
                    Generator::new([
                        0.95590731580453892,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                    ]),
                    [
                        -12.964995575216268,
                        -0.10562490738713654,
                        0.12031400988328041,
                        3.3252689261333517,
                        0.38497484079299993,
                    ],
                ),
                Data::new(
                    Generator::new([
                        0.40605717473823955,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                    ]),
                    [
                        -3.6069696334474948,
                        -0.10916821211056493,
                        0.22854417413626860,
                        4.0528205618095638,
                        5.6358900522617106E-002,
                    ],
                ),
            ];
            let final_should_be = [
                Data::new(
                    Generator::new([
                        0.89525470925235517,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                    ]),
                    [
                        3.6923766074222519,
                        3.4999626602143334E-002,
                        -5.7693384490972686E-002,
                        -1.3868627719278135,
                        -5.7693384490972686E-002,
                    ],
                ),
                Data::new(
                    Generator::new([
                        0.25000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                    ]),
                    [
                        0.0000000000000000,
                        0.0000000000000000,
                        -0.23866687325750005,
                        0.0000000000000000,
                        0.0000000000000000,
                    ],
                ),
                Data::new(
                    Generator::new([
                        0.95590731580453892,
                        0.95590731580453892,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                    ]),
                    [
                        -0.36595977174160310,
                        -9.3588053693501299E-003,
                        9.3588053693501282E-003,
                        8.6599078440979652E-002,
                        5.7181214334625485E-003,
                    ],
                ),
                Data::new(
                    Generator::new([
                        0.95590731580453892,
                        0.40605717473823955,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                    ]),
                    [
                        1.4427612574194517,
                        1.3675773263272826E-002,
                        -1.3675773263272845E-002,
                        -0.32544759695960152,
                        -2.2543144647178933E-002,
                    ],
                ),
                Data::new(
                    Generator::new([
                        0.95590731580453892,
                        0.95590731580453892,
                        0.95590731580453892,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                    ]),
                    [
                        0.11333620645370464,
                        1.0743012775049343E-003,
                        -1.0743012775049356E-003,
                        1.7708782258391337E-003,
                        -1.7708782258391350E-003,
                    ],
                ),
                Data::new(
                    Generator::new([
                        0.68607579756175630,
                        0.68607579756175630,
                        0.68607579756175630,
                        0.68607579756175630,
                        0.68607579756175630,
                        0.68607579756175630,
                    ]),
                    [
                        0.25150011495314784,
                        -9.8242232403573243E-004,
                        9.8242232403573590E-004,
                        3.9296892961429384E-003,
                        -3.9296892961429349E-003,
                    ],
                ),
            ];
            assert_check_vec_data_tol(&initial_data, &initial_should_be, tol);
            assert_check_vec_data_tol(&final_data, &final_should_be, tol);
        }

        {
            let tol = 1e-12;
            const NDIM: usize = 15;
            let rule = generate_rule::<NDIM>().unwrap();

            let initial_data = rule.initial_data();
            let final_data = rule.final_data();

            let initial_should_be = [
                Data::new(
                    Generator::new([
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                    ]),
                    [
                        291277.61628182390,
                        6.0258428637611914,
                        -2.2069758650267364,
                        -250.32473522703492,
                        -13.360990958851707,
                    ],
                ),
                Data::new(
                    Generator::new([
                        0.95590731580453892,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                    ]),
                    [
                        -6116.8154397443632,
                        1.0043900473143297E-002,
                        4.6452020230009494E-003,
                        7.9433003401195883,
                        0.36906717798885058,
                    ],
                ),
                Data::new(
                    Generator::new([
                        0.40605717473823955,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                    ]),
                    [
                        -15143.256200702785,
                        -0.35533213084947590,
                        0.47470809287517979,
                        9.9108773070823926,
                        0.46213550417183791,
                    ],
                ),
            ];
            let final_should_be = [
                Data::new(
                    Generator::new([
                        0.89525470925235517,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                    ]),
                    [
                        1890.4968230001930,
                        3.4999626602143334E-002,
                        -5.7693384490972686E-002,
                        -1.3868627719278135,
                        -5.7693384490972686E-002,
                    ],
                ),
                Data::new(
                    Generator::new([
                        0.25000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                    ]),
                    [
                        0.0000000000000000,
                        0.0000000000000000,
                        -0.23866687325750005,
                        0.0000000000000000,
                        0.0000000000000000,
                    ],
                ),
                Data::new(
                    Generator::new([
                        0.95590731580453892,
                        0.95590731580453892,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                    ]),
                    [
                        -1231.8778818090427,
                        -2.8696228364438955E-002,
                        2.8696228364438969E-002,
                        5.4723270375875238E-002,
                        3.7593929498566978E-002,
                    ],
                ),
                Data::new(
                    Generator::new([
                        0.95590731580453892,
                        0.40605717473823955,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                    ]),
                    [
                        738.69376379875928,
                        1.3675773263272826E-002,
                        -1.3675773263272845E-002,
                        -0.32544759695960152,
                        -2.2543144647178933E-002,
                    ],
                ),
                Data::new(
                    Generator::new([
                        0.95590731580453892,
                        0.95590731580453892,
                        0.95590731580453892,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                        0.0000000000000000,
                    ]),
                    [
                        58.028137704296775,
                        1.0743012775049343E-003,
                        -1.0743012775049356E-003,
                        1.7708782258391337E-003,
                        -1.7708782258391350E-003,
                    ],
                ),
                Data::new(
                    Generator::new([
                        0.68607579756175630,
                        0.68607579756175630,
                        0.68607579756175630,
                        0.68607579756175630,
                        0.68607579756175630,
                        0.68607579756175630,
                        0.68607579756175630,
                        0.68607579756175630,
                        0.68607579756175630,
                        0.68607579756175630,
                        0.68607579756175630,
                        0.68607579756175630,
                        0.68607579756175630,
                        0.68607579756175630,
                        0.68607579756175630,
                    ]),
                    [
                        0.25150011495314784,
                        -1.9187936016322899E-006,
                        1.9187936016322967E-006,
                        7.6751744065291766E-006,
                        -7.6751744065291698E-006,
                    ],
                ),
            ];
            assert_check_vec_data_tol(&initial_data, &initial_should_be, tol);
            assert_check_vec_data_tol(&final_data, &final_should_be, tol);
        }
    }
}
