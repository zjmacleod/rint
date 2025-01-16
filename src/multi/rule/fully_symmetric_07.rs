use crate::multi::generator::Generator;
use crate::multi::rule::{
    Norms, Scales, Weights, LAMBDA0, LAMBDA0_SQ, LAMBDA1, LAMBDA1_SQ, LAMBDA2, LAMBDA2_SQ, LAMBDAP,
    LAMBDAP_SQ,
};
use crate::multi::{two_pow_n, two_pow_n_f64};

pub(crate) const fn evaluations<const NDIM: usize>() -> usize {
    1 + 6 * NDIM + 2 * NDIM * (NDIM - 1) + two_pow_n(NDIM)
}

pub(crate) const fn error_coefficients() -> [f64; 6] {
    ERROR_COEFFICIENTS
}

pub(crate) const fn rule_points<const NDIM: usize>() -> [f64; WEIGHTS_LENGTH] {
    let ndim = NDIM as f64;
    [
        1.0,
        2.0 * ndim,
        2.0 * ndim,
        2.0 * ndim,
        2.0 * ndim * (ndim - 1.0),
        two_pow_n_f64(NDIM),
    ]
}

const fn initial_weights_1<const NDIM: usize>() -> [f64; WEIGHTS_LENGTH] {
    let two_ndim = two_pow_n_f64(NDIM);

    let f1 = 1f64 / ((3f64 * LAMBDA0_SQ) * (3f64 * LAMBDA0_SQ) * (3f64 * LAMBDA0_SQ)) / two_ndim;

    let e1 = (1f64 - 5f64 * LAMBDA0_SQ / 3f64)
        / (60f64 * (LAMBDA1_SQ - LAMBDA0_SQ) * LAMBDA1_SQ * LAMBDA1_SQ);

    let d1 = 0f64;

    let c1 = (1f64
        - (5f64 * LAMBDA2_SQ / 3f64)
        - 5f64 * two_ndim * f1 * LAMBDA0_SQ * (LAMBDA0_SQ - LAMBDA2_SQ))
        / (10f64 * LAMBDA1_SQ * (LAMBDA1_SQ - LAMBDA2_SQ))
        - 2f64 * ((NDIM - 1) as f64) * e1;

    let b1 = (1f64
        - (5f64 * LAMBDA1_SQ / 3f64)
        - 5f64 * two_ndim * f1 * LAMBDA0_SQ * (LAMBDA0_SQ - LAMBDA1_SQ))
        / (10f64 * LAMBDA2_SQ * (LAMBDA2_SQ - LAMBDA1_SQ));

    let a1 = 0.0;

    [a1, b1, c1, d1, e1, f1]
}

const fn initial_weights_2<const NDIM: usize>() -> [f64; WEIGHTS_LENGTH] {
    let two_ndim = two_pow_n_f64(NDIM);

    let f2 = 1f64 / (36f64 * LAMBDA0_SQ * LAMBDA0_SQ * LAMBDA0_SQ) / two_ndim;

    let e2 =
        (1f64 - 9f64 * two_ndim * f2 * LAMBDA0_SQ * LAMBDA0_SQ) / (36f64 * LAMBDA1_SQ * LAMBDA1_SQ);

    let d2 = 0f64;

    let c2 = (1f64
        - (5f64 * LAMBDA2_SQ / 3f64)
        - 5f64 * two_ndim * f2 * LAMBDA0_SQ * (LAMBDA0_SQ - LAMBDA2_SQ))
        / (10f64 * LAMBDA1_SQ * (LAMBDA1_SQ - LAMBDA2_SQ))
        - 2f64 * ((NDIM - 1) as f64) * e2;

    let b2 = (1f64
        - (5f64 * LAMBDA1_SQ / 3f64)
        - 5f64 * two_ndim * f2 * LAMBDA0_SQ * (LAMBDA0_SQ - LAMBDA1_SQ))
        / (10f64 * LAMBDA2_SQ * (LAMBDA2_SQ - LAMBDA1_SQ));

    let a2 = 0f64;

    [a2, b2, c2, d2, e2, f2]
}

const fn initial_weights_3<const NDIM: usize>() -> [f64; WEIGHTS_LENGTH] {
    let two_ndim = two_pow_n_f64(NDIM);

    let f3 = 5f64 / (108f64 * LAMBDA0_SQ * LAMBDA0_SQ * LAMBDA0_SQ) / two_ndim;

    let e3 =
        (1f64 - 9f64 * two_ndim * f3 * LAMBDA0_SQ * LAMBDA0_SQ) / (36f64 * LAMBDA1_SQ * LAMBDA1_SQ);

    let d3 = (1f64
        - (5f64 * LAMBDA1_SQ / 3f64)
        - 5f64 * two_ndim * f3 * LAMBDA0_SQ * (LAMBDA0_SQ - LAMBDA1_SQ))
        / (10f64 * LAMBDAP_SQ * (LAMBDAP_SQ - LAMBDA1_SQ));

    let c3 = (1f64
        - (5f64 * LAMBDAP_SQ / 3f64)
        - 5f64 * two_ndim * f3 * LAMBDA0_SQ * (LAMBDA0_SQ - LAMBDAP_SQ))
        / (10f64 * LAMBDA1_SQ * (LAMBDA1_SQ - LAMBDAP_SQ))
        - 2f64 * ((NDIM - 1) as f64) * e3;

    let b3 = 0f64;

    let a3 = 0f64;

    [a3, b3, c3, d3, e3, f3]
}

const fn initial_weights_4<const NDIM: usize>() -> [f64; WEIGHTS_LENGTH] {
    let two_ndim = two_pow_n_f64(NDIM);

    let f4 = 1f64 / (54f64 * LAMBDA0_SQ * LAMBDA0_SQ * LAMBDA0_SQ) / two_ndim;

    let e4 = (1f64 - 18f64 * two_ndim * f4 * LAMBDA0_SQ * LAMBDA0_SQ)
        / (72f64 * LAMBDA1_SQ * LAMBDA1_SQ);

    let d4 = 0f64;

    let c4 = (1f64
        - (10f64 * LAMBDA2_SQ / 3f64)
        - 10f64 * two_ndim * f4 * LAMBDA0_SQ * (LAMBDA0_SQ - LAMBDA2_SQ))
        / (20f64 * LAMBDA1_SQ * (LAMBDA1_SQ - LAMBDA2_SQ))
        - 2f64 * ((NDIM - 1) as f64) * e4;

    let b4 = (1f64
        - (10f64 * LAMBDA1_SQ / 3f64)
        - 10f64 * two_ndim * f4 * LAMBDA0_SQ * (LAMBDA0_SQ - LAMBDA1_SQ))
        / (20f64 * LAMBDA2_SQ * (LAMBDA2_SQ - LAMBDA1_SQ));

    let a4 = 0f64;

    [a4, b4, c4, d4, e4, f4]
}

const fn initial_weights_5<const NDIM: usize>() -> [f64; WEIGHTS_LENGTH] {
    let a5 = 0.0;
    let b5 = 0.0;
    let c5 = 0.0;
    let d5 = 0.0;
    let e5 = 0.0;
    let f5 = 0.0;

    [a5, b5, c5, d5, e5, f5]
}

const fn initial_weights<const NDIM: usize>() -> [[f64; 5]; WEIGHTS_LENGTH] {
    let initial = [
        initial_weights_1::<NDIM>(),
        initial_weights_2::<NDIM>(),
        initial_weights_3::<NDIM>(),
        initial_weights_4::<NDIM>(),
        initial_weights_5::<NDIM>(),
    ];

    let mut output = [[0.0; 5]; WEIGHTS_LENGTH];

    let mut j = 0;
    while j < 5 {
        let mut i = 0;
        while i < WEIGHTS_LENGTH {
            output[i][j] = initial[j][i];
            i += 1;
        }
        j += 1;
    }

    output
}

pub(crate) const fn weights<const NDIM: usize>() -> Weights<WEIGHTS_LENGTH> {
    let mut weights = initial_weights::<NDIM>();
    let rule_points = rule_points::<NDIM>();

    let ndim = NDIM as f64;
    let two_ndim = two_pow_n_f64(NDIM);

    weights[0][0] = two_ndim;

    let mut i = 1;
    while i < WEIGHTS_LENGTH {
        let mut j = 1;
        while j < 5 {
            weights[i][j] -= weights[i][0];
            weights[0][j] -= rule_points[i] * weights[i][j];
            j += 1;
        }
        i += 1;
    }

    let mut k = 1;
    while k < WEIGHTS_LENGTH {
        weights[k][0] *= two_ndim;
        weights[0][0] -= rule_points[k] * weights[k][0];
        k += 1;
    }

    Weights(weights)
}

pub(crate) const fn generators<const NDIM: usize>() -> [Generator<NDIM>; WEIGHTS_LENGTH] {
    let gen_1 = Generator::new([0.0; NDIM]);
    let gen_2 = {
        let mut init = [0.0; NDIM];
        init[0] = LAMBDA2;
        Generator::new(init)
    };
    let gen_3 = {
        let mut init = [0.0; NDIM];
        init[0] = LAMBDA1;
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
    let gen_6 = Generator::new([LAMBDA0; NDIM]);

    [gen_1, gen_2, gen_3, gen_4, gen_5, gen_6]
}

pub(crate) const WEIGHTS_LENGTH: usize = 6;

pub(crate) const ERROR_COEFFICIENTS: [f64; 6] = [5.0, 5.0, 1.0, 5.0, 0.5, 0.25];

#[cfg(test)]
mod tests {
    use super::*;

    fn rel_or_abs_diff(a: f64, b: f64) -> f64 {
        if a == 0.0 {
            (a - b).abs()
        } else {
            (a - b).abs() / a.abs()
        }
    }

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

    #[test]
    fn check_evaluations_correct() {
        assert_eq!(21, evaluations::<2>());
        assert_eq!(39, evaluations::<3>());
        assert_eq!(65, evaluations::<4>());
        assert_eq!(103, evaluations::<5>());
        assert_eq!(161, evaluations::<6>());
        assert_eq!(255, evaluations::<7>());
        assert_eq!(417, evaluations::<8>());
        assert_eq!(711, evaluations::<9>());
        assert_eq!(1265, evaluations::<10>());
    }

    #[test]
    fn check_rule_points_correct() {
        {
            const NDIM: usize = 2;
            let rule_points = rule_points::<NDIM>();
            let should_be = [1.0, 4.0, 4.0, 4.0, 4.0, 4.0];
            for (x, y) in rule_points.iter().zip(should_be.iter()) {
                assert_eq!(x, y);
            }
        }

        {
            const NDIM: usize = 3;
            let rule_points = rule_points::<NDIM>();
            let should_be = [1.0, 6.0, 6.0, 6.0, 12.0, 8.0];
            let tol = 1e-15;
            for (x, y) in rule_points.iter().zip(should_be.iter()) {
                assert_eq!(x, y);
            }
        }

        {
            const NDIM: usize = 4;
            let rule_points = rule_points::<NDIM>();
            let should_be = [1.0, 8.0, 8.0, 8.0, 24.0, 16.0];
            for (x, y) in rule_points.iter().zip(should_be.iter()) {
                assert_eq!(x, y);
            }
        }

        {
            const NDIM: usize = 5;
            let rule_points = rule_points::<NDIM>();
            let should_be = [1.0, 10.0, 10.0, 10.0, 40.0, 32.0];
            for (x, y) in rule_points.iter().zip(should_be.iter()) {
                assert_eq!(x, y);
            }
        }

        {
            const NDIM: usize = 6;
            let rule_points = rule_points::<NDIM>();
            let should_be = [1.0, 12.0, 12.0, 12.0, 60.0, 64.0];
            for (x, y) in rule_points.iter().zip(should_be.iter()) {
                assert_eq!(x, y);
            }
        }

        {
            const NDIM: usize = 15;
            let rule_points = rule_points::<NDIM>();
            let should_be = [1.0, 30.0, 30.0, 30.0, 420.0, 32768.0];
            for (x, y) in rule_points.iter().zip(should_be.iter()) {
                assert_eq!(x, y);
            }
        }
    }

    #[test]
    fn check_generators_correct() {
        let tol = 1e-15;

        {
            const NDIM: usize = 2;
            let generators = generators::<NDIM>();
            let should_be = [
                Generator::new([0.0000000000000000, 0.0000000000000000]),
                Generator::new([0.40605717473823955, 0.0000000000000000]),
                Generator::new([0.95590731580453892, 0.0000000000000000]),
                Generator::new([0.75000000000000000, 0.0000000000000000]),
                Generator::new([0.95590731580453892, 0.95590731580453892]),
                Generator::new([0.68607579756175630, 0.68607579756175630]),
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
            const NDIM: usize = 3;
            let generators = generators::<NDIM>();
            let should_be = [
                Generator::new([0.0000000000000000, 0.0000000000000000, 0.0000000000000000]),
                Generator::new([0.40605717473823955, 0.0000000000000000, 0.0000000000000000]),
                Generator::new([0.95590731580453892, 0.0000000000000000, 0.0000000000000000]),
                Generator::new([0.75000000000000000, 0.0000000000000000, 0.0000000000000000]),
                Generator::new([0.95590731580453892, 0.95590731580453892, 0.0000000000000000]),
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
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                ]),
                Generator::new([
                    0.40605717473823955,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                ]),
                Generator::new([
                    0.95590731580453892,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                ]),
                Generator::new([
                    0.75000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                ]),
                Generator::new([
                    0.95590731580453892,
                    0.95590731580453892,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
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
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                ]),
                Generator::new([
                    0.40605717473823955,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                ]),
                Generator::new([
                    0.95590731580453892,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                ]),
                Generator::new([
                    0.75000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                ]),
                Generator::new([
                    0.95590731580453892,
                    0.95590731580453892,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
                    0.00000000000000000,
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
    fn check_weights_correct() {
        {
            let tol = 1e-14;
            const NDIM: usize = 2;
            let weights = weights::<NDIM>();
            let should_be = Weights([
                [
                    -0.34070799859740308,
                    -0.22492123910230477,
                    0.33757853071887722,
                    -1.8438476127784753,
                    1.0851769996493508,
                ],
                [
                    0.49437592128172986,
                    7.4978148702033634E-002,
                    -0.12359398032043246,
                    0.55489147051423582,
                    -0.12359398032043246,
                ],
                [
                    0.19682203269278409,
                    -2.4412953600027321E-003,
                    -7.4889249639626052E-003,
                    -4.4682186610262520E-002,
                    -4.9205508173196022E-002,
                ],
                [
                    0.0000000000000000,
                    0.0000000000000000,
                    3.0381729038221013E-002,
                    0.0000000000000000,
                    0.0000000000000000,
                ],
                [
                    3.8835733349496797E-002,
                    5.8899134538790317E-003,
                    -5.8899134538790369E-003,
                    -4.8544666686871005E-003,
                    -9.7089333373741992E-003,
                ],
                [
                    0.35514331232534008,
                    -2.2196457020333737E-002,
                    2.2196457020333793E-002,
                    -4.4392914040667496E-002,
                    -8.8785828081335019E-002,
                ],
            ]);
            assert_check_vec_tol(&weights.0, &should_be.0, tol);
        }

        {
            let tol = 1e-14;
            const NDIM: usize = 3;
            let weights = weights::<NDIM>();
            let should_be = Weights([
                [
                    -3.1355219462968886,
                    -0.34643529197085043,
                    0.51542122939570922,
                    -2.8836840472611702,
                    1.3919402432871111,
                ],
                [
                    0.98875184256345972,
                    7.4978148702033634E-002,
                    -0.12359398032043246,
                    0.55489147051423582,
                    -0.12359398032043246,
                ],
                [
                    0.23830113198758102,
                    -1.4221122267760802E-002,
                    4.2909019437954651E-003,
                    -3.4973253272888326E-002,
                    -2.9787641498447627E-002,
                ],
                [
                    0.0000000000000000,
                    0.0000000000000000,
                    3.0381729038221013E-002,
                    0.0000000000000000,
                    0.0000000000000000,
                ],
                [
                    7.7671466698993594E-002,
                    5.8899134538790317E-003,
                    -5.8899134538790369E-003,
                    -4.8544666686871005E-003,
                    -9.7089333373741992E-003,
                ],
                [
                    0.35514331232534008,
                    -1.1098228510166869E-002,
                    1.1098228510166896E-002,
                    -2.2196457020333748E-002,
                    -4.4392914040667510E-002,
                ],
            ]);
            assert_check_vec_tol(&weights.0, &should_be.0, tol);
        }

        {
            let tol = 1e-14;
            const NDIM: usize = 6;
            let weights = weights::<NDIM>();
            let should_be = Weights([
                [
                    -69.069796742618294,
                    -0.56961952768339064,
                    0.90759140253310810,
                    -6.1197005507577442,
                    2.0792155741034111,
                ],
                [
                    7.9100147405076777,
                    7.4978148702033634E-002,
                    -0.12359398032043246,
                    0.55489147051423582,
                    -0.12359398032043246,
                ],
                [
                    -1.8218213456510446,
                    -4.9560602991034999E-002,
                    3.9630382667069697E-002,
                    -5.8464532607657230E-003,
                    2.8465958525797572E-002,
                ],
                [
                    0.0000000000000000,
                    0.0000000000000000,
                    3.0381729038221013E-002,
                    0.0000000000000000,
                    0.0000000000000000,
                ],
                [
                    0.62137173359194875,
                    5.8899134538790317E-003,
                    -5.8899134538790369E-003,
                    -4.8544666686871005E-003,
                    -9.7089333373741992E-003,
                ],
                [
                    0.35514331232534008,
                    -1.3872785637708586E-003,
                    1.3872785637708621E-003,
                    -2.7745571275417185E-003,
                    -5.5491142550834387E-003,
                ],
            ]);
            assert_check_vec_tol(&weights.0, &should_be.0, tol);
        }

        {
            let tol = 1e-13;
            const NDIM: usize = 15;
            let weights = weights::<NDIM>();
            let should_be = Weights([
                [
                    -34206.907291385309,
                    3.3049071216860737E-002,
                    0.81188061590743388,
                    -16.876314861683888,
                    2.0439119656794835,
                ],
                [
                    4049.9275471399310,
                    7.4978148702033634E-002,
                    -0.12359398032043246,
                    0.55489147051423582,
                    -0.12359398032043246,
                ],
                [
                    -6659.3344257567342,
                    -0.15557904516085760,
                    0.14564882483689234,
                    8.1533946775602092E-002,
                    0.20322675859853315,
                ],
                [
                    0.0000000000000000,
                    0.0000000000000000,
                    3.0381729038221013E-002,
                    0.0000000000000000,
                    0.0000000000000000,
                ],
                [
                    318.14232759907776,
                    5.8899134538790317E-003,
                    -5.8899134538790369E-003,
                    -4.8544666686871005E-003,
                    -9.7089333373741992E-003,
                ],
                [
                    0.35514331232534008,
                    -2.7095284448649582E-006,
                    2.7095284448649650E-006,
                    -5.4190568897299189E-006,
                    -1.0838113779459841E-005,
                ],
            ]);
            assert_check_vec_tol(&weights.0, &should_be.0, tol);
        }
    }

    #[test]
    fn check_scales_norms_correct() {
        {
            let tol = 1e-14;
            const NDIM: usize = 2;
            let weights = weights::<NDIM>();
            let rule_points = rule_points::<NDIM>();
            let (scales, norms) = weights.scales_norms::<NDIM>(&rule_points);
            let scales_should_be = Scales([
                [1.5008744041523381, 5.4619812724819363, 0.58853941731882520],
                [1.6484000000000030, 4.4896318499947325, 0.22273541203632835],
                [
                    -3.0676029974325698,
                    -5.9664353462315765,
                    -1.1012332185617477,
                ],
                [100.00000000000000, -0.0000000000000000, 100.00000000000000],
                [
                    1.0000000000000009,
                    -0.82419999999999971,
                    -1.9999999999999996,
                ],
                [1.0000000000000024, 1.9999999999999960, -2.0000000000000004],
            ]);
            let norms_should_be = Norms([
                [15.000673861797612, 2.0594085530307291, 2.4632784012088300],
                [14.619562897422956, 2.6087348779812332, 2.9652172226381066],
                [1.3243103526083444, 0.37804488185591240, 0.64191377794942761],
                [
                    6.2596209976816711E-002,
                    0.90107710528805529,
                    9.0308860340991368E-003,
                ],
                [8.5402855258664925, 0.76131586120970718, 0.40539107822219211],
                [8.5402855258665067, 1.3569737411797032, 0.40539107822219189],
            ]);
            assert_check_vec_tol(&scales.0, &scales_should_be.0, tol);
            assert_check_vec_tol(&norms.0, &norms_should_be.0, tol);
        }

        {
            let tol = 1e-14;
            const NDIM: usize = 3;
            let weights = weights::<NDIM>();
            let rule_points = rule_points::<NDIM>();
            let (scales, norms) = weights.scales_norms::<NDIM>(&rule_points);
            let scales_should_be = Scales([
                [1.4877849957592586, 5.5948103857539246, 0.48269512903437911],
                [1.6484000000000030, 4.4896318499947325, 0.22273541203632835],
                [
                    0.30172737868395333,
                    8.1505598895026630,
                    -0.85172635402321450,
                ],
                [100.00000000000000, -0.0000000000000000, 100.00000000000000],
                [
                    1.0000000000000009,
                    -0.82419999999999971,
                    -1.9999999999999996,
                ],
                [1.0000000000000024, 1.9999999999999960, -2.0000000000000004],
            ]);
            let norms_should_be = Norms([
                [18.453034922164381, 2.9871962141915933, 4.6216242818616493],
                [17.534737522530300, 3.8481493502257038, 5.3358825416997862],
                [6.1051989517188385, 1.1943369141825844, 1.0394875480397585],
                [
                    7.7783432141282596E-002,
                    1.2014361403840739,
                    1.2041181378798850E-002,
                ],
                [11.387047367821989, 1.0150878149462763, 0.54052143762958937],
                [11.387047367822007, 1.8092983215729377, 0.54052143762958926],
            ]);
            assert_check_vec_tol(&scales.0, &scales_should_be.0, tol);
            assert_check_vec_tol(&norms.0, &norms_should_be.0, tol);
        }

        {
            let tol = 1e-14;
            const NDIM: usize = 6;
            let weights = weights::<NDIM>();
            let rule_points = rule_points::<NDIM>();
            let (scales, norms) = weights.scales_norms::<NDIM>(&rule_points);
            let scales_should_be = Scales([
                [1.5933291581913094, 6.7427925536508191, 0.33975773109453233],
                [1.6484000000000030, 4.4896318499947325, 0.22273541203632835],
                [0.79963479609476151, 0.14752452202849281, 4.8689277509197639],
                [100.00000000000000, -0.0000000000000000, 100.00000000000000],
                [
                    1.0000000000000009,
                    -0.82419999999999971,
                    -1.9999999999999996,
                ],
                [1.0000000000000024, 1.9999999999999960, -2.0000000000000004],
            ]);
            let norms_should_be = Norms([
                [55.723869508898815, 5.3192717756037737, 29.171230471345883],
                [53.897283837715740, 8.1573764597229559, 30.707034067695144],
                [38.347503188392174, 4.9278648878824640, 1.0343404256680435],
                [
                    0.25839518414799995,
                    4.8057445615362964,
                    4.8164725515195399E-002,
                ],
                [45.548189471287955, 4.0603512597851052, 2.1620857505183579],
                [45.548189471287984, 6.0348034352727939, 2.1620857505183575],
            ]);
            assert_check_vec_tol(&scales.0, &scales_should_be.0, tol);
            assert_check_vec_tol(&norms.0, &norms_should_be.0, tol);
        }

        {
            let tol = 1e-13;
            const NDIM: usize = 15;
            let weights = weights::<NDIM>();
            let rule_points = rule_points::<NDIM>();
            let (scales, norms) = weights.scales_norms::<NDIM>(&rule_points);
            let scales_should_be = Scales([
                [-24.565913231874259, 20.786695150765901, 0.12111127236195365],
                [1.6484000000000030, 4.4896318499947325, 0.22273541203632835],
                [
                    0.93617250759124970,
                    -0.55979817802793441,
                    -2.4925416545559149,
                ],
                [100.00000000000000, -0.0000000000000000, 100.00000000000000],
                [
                    1.0000000000000009,
                    -0.82419999999999971,
                    -1.9999999999999996,
                ],
                [1.0000000000000024, 1.9999999999999960, -2.0000000000000004],
            ]);
            let norms_should_be = Norms([
                [134.06536084307373, 143.86213251118596, 2562.7857645452909],
                [4844.7587668888145, 621.16054365500952, 2466.8688867813166],
                [9309.4130587031723, 875.10266687973365, 362.47368576565020],
                [34.773936788842335, 858.12621321598749, 8.5705382397315226],
                [9328.2692037197648, 831.55993800398949, 442.79516170615960],
                [9328.2692037197430, 736.72520029221721, 442.79516170615943],
            ]);
            assert_check_vec_tol(&scales.0, &scales_should_be.0, tol);
            assert_check_vec_tol(&norms.0, &norms_should_be.0, tol);
        }
    }
}
