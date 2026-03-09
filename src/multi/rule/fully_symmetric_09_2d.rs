use crate::multi::generator::Generator;
use crate::multi::rule::{
    scales_norms, BasicErrorCoeff, Data, Rule, ScalesNorms, ADAPTIVE_ERROR_COEFF,
};

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

const TOTAL: usize = 8;
const FINAL: usize = TOTAL - 3;
const RATIO: f64 = (LAMBDA2 / LAMBDA1) * (LAMBDA2 / LAMBDA1);
const EVALUATIONS: usize = 33;

pub(crate) const fn generate_rule() -> Rule<2, FINAL, TOTAL> {
    let evaluations = EVALUATIONS;
    let basic_error_coeff = BASIC_ERROR_COEFF;
    let adaptive_error_coeff = ADAPTIVE_ERROR_COEFF;
    let ratio = RATIO;

    let (weights, scales_norms) = weights_scales_norms();

    let [g0, g1, g2, g3, g4, g5, g6, g7] = GENERATORS;
    let [w0, w1, w2, w3, w4, w5, w6, w7] = weights;

    let data0 = Data::new(g0, w0);
    let data1 = Data::new(g1, w1);
    let data2 = Data::new(g2, w2);
    let data3 = Data::new(g3, w3);
    let data4 = Data::new(g4, w4);
    let data5 = Data::new(g5, w5);
    let data6 = Data::new(g6, w6);
    let data7 = Data::new(g7, w7);

    let initial_data = [data0, data1, data2];
    let final_data = [data3, data4, data5, data6, data7];

    Rule {
        initial_data,
        final_data,
        scales_norms,
        basic_error_coeff,
        adaptive_error_coeff,
        evaluations,
        ratio,
    }
}

const fn rule_points() -> [f64; TOTAL] {
    [1.0, 4.0, 4.0, 4.0, 4.0, 4.0, 8.0, 4.0]
}

const fn initial_weights_1() -> [f64; TOTAL] {
    let two_ndim = 4.0;
    let ndim = 2.0;

    let w7 = 1.0
        / ((3.0 * LAMBDA0_SQ) * (3.0 * LAMBDA0_SQ) * (3.0 * LAMBDA0_SQ) * (3.0 * LAMBDA0_SQ))
        / two_ndim;

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

    [w0, w1, w2, w3, w4, w5, w6, w7]
}

const fn initial_weights_2() -> [f64; TOTAL] {
    let two_ndim = 4.0;
    let ndim = 2.0;

    let w7 = 1.0 / (108.0 * LAMBDA0_SQ * LAMBDA0_SQ * LAMBDA0_SQ * LAMBDA0_SQ) / two_ndim;

    let w6 = (1.0
        - 5.0 * LAMBDA1_SQ / 3.0
        - 15.0 * two_ndim * w7 * LAMBDA0_SQ * LAMBDA0_SQ * (LAMBDA0_SQ - LAMBDA1_SQ))
        / (60.0 * LAMBDA1_SQ * LAMBDA2_SQ * (LAMBDA2_SQ - LAMBDA1_SQ));

    let w5 = (1.0
        - 9.0 * (8.0 * LAMBDA1_SQ * LAMBDA2_SQ * w6 + two_ndim * w7 * LAMBDA0_SQ * LAMBDA0_SQ))
        / (36.0 * LAMBDA1_SQ * LAMBDA1_SQ)
        - 2.0 * w7 * (ndim - 2.0);

    let w4 = 0.0;

    let w3 = (1.0
        - 7.0
            * ((LAMBDA1_SQ + LAMBDA2_SQ) / 5.0 - LAMBDA1_SQ * LAMBDA2_SQ / 3.0
                + two_ndim
                    * w7
                    * LAMBDA0_SQ
                    * (LAMBDA0_SQ - LAMBDA1_SQ)
                    * (LAMBDA0_SQ - LAMBDA2_SQ)))
        / (14.0 * LAMBDA3_SQ * (LAMBDA3_SQ - LAMBDA1_SQ) * (LAMBDA3_SQ - LAMBDA2_SQ));

    let w2 = (1.0
        - 7.0
            * ((LAMBDA1_SQ + LAMBDA3_SQ) / 5.0 - LAMBDA1_SQ * LAMBDA3_SQ / 3.0
                + two_ndim
                    * w7
                    * LAMBDA0_SQ
                    * (LAMBDA0_SQ - LAMBDA1_SQ)
                    * (LAMBDA0_SQ - LAMBDA3_SQ)))
        / (14.0 * LAMBDA2_SQ * (LAMBDA2_SQ - LAMBDA1_SQ) * (LAMBDA2_SQ - LAMBDA3_SQ))
        - 2.0 * (ndim - 1.0) * w6;

    let w1 = (1.0
        - 7.0
            * ((LAMBDA2_SQ + LAMBDA3_SQ) / 5.0 - LAMBDA2_SQ * LAMBDA3_SQ / 3.0
                + two_ndim
                    * w7
                    * LAMBDA0_SQ
                    * (LAMBDA0_SQ - LAMBDA2_SQ)
                    * (LAMBDA0_SQ - LAMBDA3_SQ)))
        / (14.0 * LAMBDA1_SQ * (LAMBDA1_SQ - LAMBDA2_SQ) * (LAMBDA1_SQ - LAMBDA3_SQ))
        - 2.0 * (ndim - 1.0) * (w6 + w5 + (ndim - 2.0) * w7);

    let w0 = 0.0;

    [w0, w1, w2, w3, w4, w5, w6, w7]
}

const fn initial_weights_3() -> [f64; TOTAL] {
    let two_ndim = 4.0;
    let ndim = 2.0;

    let w7 = 5.0 / (324.0 * LAMBDA0_SQ * LAMBDA0_SQ * LAMBDA0_SQ * LAMBDA0_SQ) / two_ndim;

    let w6 = (1.0
        - 5.0 * LAMBDA1_SQ / 3.0
        - 15.0 * two_ndim * w7 * LAMBDA0_SQ * LAMBDA0_SQ * (LAMBDA0_SQ - LAMBDA1_SQ))
        / (60.0 * LAMBDA1_SQ * LAMBDA2_SQ * (LAMBDA2_SQ - LAMBDA1_SQ));

    let w5 = (1.0
        - 9.0 * (8.0 * LAMBDA1_SQ * LAMBDA2_SQ * w6 + two_ndim * w7 * LAMBDA0_SQ * LAMBDA0_SQ))
        / (36.0 * LAMBDA1_SQ * LAMBDA1_SQ)
        - 2.0 * w7 * (ndim - 2.0);

    let w4 = (1.0
        - 7.0
            * ((LAMBDA1_SQ + LAMBDA2_SQ) / 5.0 - LAMBDA1_SQ * LAMBDA2_SQ / 3.0
                + two_ndim
                    * w7
                    * LAMBDA0_SQ
                    * (LAMBDA0_SQ - LAMBDA1_SQ)
                    * (LAMBDA0_SQ - LAMBDA2_SQ)))
        / (14.0 * LAMBDAP_SQ * (LAMBDAP_SQ - LAMBDA1_SQ) * (LAMBDAP_SQ - LAMBDA2_SQ));

    let w3 = 0.0;

    let w2 = (1.0
        - 7.0
            * ((LAMBDA1_SQ + LAMBDAP_SQ) / 5.0 - LAMBDA1_SQ * LAMBDAP_SQ / 3.0
                + two_ndim
                    * w7
                    * LAMBDA0_SQ
                    * (LAMBDA0_SQ - LAMBDA1_SQ)
                    * (LAMBDA0_SQ - LAMBDAP_SQ)))
        / (14.0 * LAMBDA2_SQ * (LAMBDA2_SQ - LAMBDA1_SQ) * (LAMBDA2_SQ - LAMBDAP_SQ))
        - 2.0 * (ndim - 1.0) * w6;

    let w1 = (1.0
        - 7.0
            * ((LAMBDA2_SQ + LAMBDAP_SQ) / 5.0 - LAMBDA2_SQ * LAMBDAP_SQ / 3.0
                + two_ndim
                    * w7
                    * LAMBDA0_SQ
                    * (LAMBDA0_SQ - LAMBDA2_SQ)
                    * (LAMBDA0_SQ - LAMBDAP_SQ)))
        / (14.0 * LAMBDA1_SQ * (LAMBDA1_SQ - LAMBDA2_SQ) * (LAMBDA1_SQ - LAMBDAP_SQ))
        - 2.0 * (ndim - 1.0) * (w6 + w5 + (ndim - 2.0) * w7);

    let w0 = 0.0;

    [w0, w1, w2, w3, w4, w5, w6, w7]
}

const fn initial_weights_4() -> [f64; TOTAL] {
    let two_ndim = 4.0;
    let ndim = 2.0;

    let w7 = 2.0 / (81.0 * LAMBDA0_SQ * LAMBDA0_SQ * LAMBDA0_SQ * LAMBDA0_SQ) / two_ndim;

    let w6 = (2.0
        - 15.0 * LAMBDA1_SQ / 9.0
        - 15.0 * two_ndim * w7 * LAMBDA0_SQ * (LAMBDA0_SQ - LAMBDA1_SQ))
        / (60.0 * LAMBDA1_SQ * LAMBDA2_SQ * (LAMBDA2_SQ - LAMBDA1_SQ));

    let w5 = (1.0
        - 9.0 * (8.0 * LAMBDA1_SQ * LAMBDA2_SQ * w6 + two_ndim * w7 * LAMBDA0_SQ * LAMBDA0_SQ))
        / (36.0 * LAMBDA1_SQ * LAMBDA1_SQ)
        - 2.0 * w7 * (ndim - 2.0);

    let w4 = 0.0;

    let w3 = (2.0
        - 7.0
            * ((LAMBDA1_SQ + LAMBDA2_SQ) / 5.0 - LAMBDA1_SQ * LAMBDA2_SQ / 3.0
                + two_ndim
                    * w7
                    * LAMBDA0_SQ
                    * (LAMBDA0_SQ - LAMBDA1_SQ)
                    * (LAMBDA0_SQ - LAMBDA2_SQ)))
        / (14.0 * LAMBDA3_SQ * (LAMBDA3_SQ - LAMBDA1_SQ) * (LAMBDA3_SQ - LAMBDA2_SQ));

    let w2 = (2.0
        - 7.0
            * ((LAMBDA1_SQ + LAMBDA3_SQ) / 5.0 - LAMBDA1_SQ * LAMBDA3_SQ / 3.0
                + two_ndim
                    * w7
                    * LAMBDA0_SQ
                    * (LAMBDA0_SQ - LAMBDA1_SQ)
                    * (LAMBDA0_SQ - LAMBDA3_SQ)))
        / (14.0 * LAMBDA2_SQ * (LAMBDA2_SQ - LAMBDA1_SQ) * (LAMBDA2_SQ - LAMBDA3_SQ))
        - 2.0 * (ndim - 1.0) * w6;

    let w1 = (2.0
        - 7.0
            * ((LAMBDA2_SQ + LAMBDA3_SQ) / 5.0 - LAMBDA2_SQ * LAMBDA3_SQ / 3.0
                + two_ndim
                    * w7
                    * LAMBDA0_SQ
                    * (LAMBDA0_SQ - LAMBDA2_SQ)
                    * (LAMBDA0_SQ - LAMBDA3_SQ)))
        / (14.0 * LAMBDA1_SQ * (LAMBDA1_SQ - LAMBDA2_SQ) * (LAMBDA1_SQ - LAMBDA3_SQ))
        - 2.0 * (ndim - 1.0) * (w6 + w5 + (ndim - 2.0) * w7);

    let w0 = 0.0;

    [w0, w1, w2, w3, w4, w5, w6, w7]
}

const fn initial_weights_5() -> [f64; TOTAL] {
    let mut weights = [0.0; TOTAL];
    weights[1] = 1.0 / (6.0 * LAMBDA1_SQ);

    weights
}

const fn initial_weights() -> [[f64; 5]; TOTAL] {
    let initial = [
        initial_weights_1(),
        initial_weights_2(),
        initial_weights_3(),
        initial_weights_4(),
        initial_weights_5(),
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

const fn weights_scales_norms() -> ([[f64; 5]; TOTAL], [ScalesNorms<TOTAL>; 3]) {
    let mut weights = initial_weights();
    let rule_points = rule_points();

    let two_ndim = 4.0;

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

    let scales_norms = scales_norms::<2, TOTAL>(&weights, rule_points);
    (weights, scales_norms)
}

pub const GENERATORS: [Generator<2>; TOTAL] = [
    Generator::new([0.0; 2]),
    Generator::new([LAMBDA1, 0.0]),
    Generator::new([LAMBDA2, 0.0]),
    Generator::new([LAMBDA3, 0.0]),
    Generator::new([LAMBDAP, 0.0]),
    Generator::new([LAMBDA1; 2]),
    Generator::new([LAMBDA1, LAMBDA2]),
    Generator::new([LAMBDA0; 2]),
];

const BASIC_ERROR_COEFF: BasicErrorCoeff = BasicErrorCoeff::new(5.0, 5.0, 1.0, 5.0);

#[cfg(test)]
#[allow(clippy::unreadable_literal)]
mod tests {
    use super::*;
    use crate::multi::rule::util;

    #[test]
    fn test_initial_weights_1() {
        let w = initial_weights_1();
        let should_be = [
            0.0000000000000000,
            -2.5476792732883849E-002,
            0.12398625665481436,
            5.7693384490972686E-002,
            0.0000000000000000,
            8.4489043732505313E-003,
            2.2543144647178933E-002,
            6.2875028738286959E-002,
        ];
        let tol = 1e-10;
        util::assert_check_vec_tol(&[w; 1], &[should_be; 1], tol);
        println!("{w:?}");
        println!("{should_be:?}");
    }

    #[test]
    fn test_initial_weights_2() {
        let w = initial_weights_2();
        let should_be = [
            0.0000000000000000,
            -7.0782726308520394E-002,
            0.12422423065043203,
            9.2693011093116021E-002,
            0.0000000000000000,
            7.6845092239398753E-003,
            3.6218917910451759E-002,
            4.7156271553715240E-002,
        ];
        let tol = 1e-10;
        util::assert_check_vec_tol(&[w; 1], &[should_be; 1], tol);
    }

    #[test]
    fn test_initial_weights_3() {
        let w = initial_weights_3();
        let should_be = [
            0.0000000000000000,
            3.4518243338896359E-002,
            0.24312424468490018,
            0.0000000000000000,
            -0.23866687325750005,
            9.2132995225611743E-003,
            8.8673713839060881E-003,
            7.8593785922858733E-002,
        ];
        let tol = 1e-10;
        util::assert_check_vec_tol(&[w; 1], &[should_be; 1], tol);
    }

    #[test]
    fn test_initial_weights_4() {
        let w = initial_weights_4();
        let should_be = [
            0.0000000000000000,
            1.4315050626716317,
            1.5732260427875659,
            -1.3291693874368409,
            0.0000000000000000,
            0.10921500862094324,
            -0.30290445231242258,
            0.12575005747657397,
        ];
        let tol = 1e-10;
        util::assert_check_vec_tol(&[w; 1], &[should_be; 1], tol);
    }

    #[test]
    fn test_initial_weights_5() {
        let w = initial_weights_5();
        let should_be = [
            0.0000000000000000,
            0.18239678493024578,
            0.0000000000000000,
            0.0000000000000000,
            0.0000000000000000,
            0.0000000000000000,
            0.0000000000000000,
            0.0000000000000000,
        ];
        let tol = 1e-10;
        util::assert_check_vec_tol(&[w; 1], &[should_be; 1], tol);
    }

    #[test]
    fn check_weights_scales_norms_correct() {
        let tol = 1e-11;
        let (weights, scales_norms) = weights_scales_norms();
        let weights_should_be = [
            [
                -0.36180913310077722,
                -3.2002448591509636E-003,
                0.51238251135707991,
                -4.1284192347049213,
                0.36086514355421107,
            ],
            [
                -0.10190717093153540,
                -4.5305933575636545E-002,
                5.9995036071780208E-002,
                1.4569818554045155,
                0.20787357766312964,
            ],
            [
                0.49594502661925743,
                2.3797399561767407E-004,
                0.11913798803008582,
                1.4492397861327515,
                -0.12398625665481436,
            ],
            [
                0.23077353796389075,
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
                3.3795617493002125E-002,
                -7.6439514931065603E-004,
                7.6439514931064302E-004,
                0.10076610424769271,
                -8.4489043732505313E-003,
            ],
            [
                9.0172578588715732E-002,
                1.3675773263272826E-002,
                -1.3675773263272845E-002,
                -0.32544759695960152,
                -2.2543144647178933E-002,
            ],
            [
                0.25150011495314784,
                -1.5718757184571719E-002,
                1.5718757184571774E-002,
                6.2875028738287014E-002,
                -6.2875028738286959E-002,
            ],
        ];
        let scales_norms_should_be = [
            ScalesNorms::new(
                [
                    160.10728363235830,
                    1.3242202805868852,
                    -500.63448201916719,
                    1.6483999999999890,
                    100.00000000000000,
                    0.99999999999998301,
                    1.0000000000000013,
                    1.0000000000000036,
                ],
                [
                    4.9726151183683655E-002,
                    1.9580414074405414,
                    1.5853458710738299E-002,
                    1.8936394668033427,
                    7.8850706471243232E-002,
                    1.9130654751138099,
                    1.9130654751138183,
                    1.9130654751138187,
                ],
            ),
            ScalesNorms::new(
                [
                    8.0572992699741501,
                    -24.285040076671180,
                    -12.164380229140484,
                    -24.038506046474300,
                    0.0000000000000000,
                    -131.82462544217731,
                    -23.797381741741187,
                    -3.9999999999999947,
                ],
                [
                    0.10761943790598244,
                    8.4662550736617992E-002,
                    0.13434064708344193,
                    8.5602187887356643E-002,
                    0.16287372978020184,
                    1.2304533447304313E-002,
                    8.6333728651175307E-002,
                    0.15422266776148044,
                ],
            ),
            ScalesNorms::new(
                [
                    8.7410004420252124E-002,
                    -0.14267410187165011,
                    8.5552617200544565E-002,
                    -4.1599922976356049E-002,
                    100.00000000000000,
                    8.3846690673704302E-002,
                    -6.9268124447012774E-002,
                    0.99999999999999911,
                ],
                [
                    1.4780688775799582,
                    1.1752910762669440,
                    1.4942719782757945,
                    1.7830714482050403,
                    1.6286706445912172E-003,
                    1.4982871677150911,
                    1.6284977942816639,
                    0.16220988716480675,
                ],
            ),
        ];
        util::assert_check_vec_tol(&weights, &weights_should_be, tol);
        for (calc, should) in scales_norms.iter().zip(scales_norms_should_be.iter()) {
            let scales = calc.scales();
            let scales_should_be = should.scales();

            let norms = calc.norms();
            let norms_should_be = should.norms();
            util::assert_check_slice_tol(scales, scales_should_be, tol);
            util::assert_check_slice_tol(norms, norms_should_be, tol);
        }
    }

    #[test]
    fn check_data_correct() {
        let tol = 1e-11;
        let rule = generate_rule();

        let initial_data = rule.initial_data();
        let final_data = rule.final_data();
        let initial_should_be = [
            Data::new(
                Generator::new([0.0000000000000000, 0.0000000000000000]),
                [
                    -0.36180913310077722,
                    -3.2002448591509636E-003,
                    0.51238251135707991,
                    -4.1284192347049213,
                    0.36086514355421107,
                ],
            ),
            Data::new(
                Generator::new([0.95590731580453892, 0.0000000000000000]),
                [
                    -0.10190717093153540,
                    -4.5305933575636545E-002,
                    5.9995036071780208E-002,
                    1.4569818554045155,
                    0.20787357766312964,
                ],
            ),
            Data::new(
                Generator::new([0.40605717473823955, 0.0000000000000000]),
                [
                    0.49594502661925743,
                    2.3797399561767407E-004,
                    0.11913798803008582,
                    1.4492397861327515,
                    -0.12398625665481436,
                ],
            ),
        ];
        let final_should_be = [
            Data::new(
                Generator::new([0.89525470925235517, 0.0000000000000000]),
                [
                    0.23077353796389075,
                    3.4999626602143334E-002,
                    -5.7693384490972686E-002,
                    -1.3868627719278135,
                    -5.7693384490972686E-002,
                ],
            ),
            Data::new(
                Generator::new([0.25000000000000000, 0.0000000000000000]),
                [
                    0.0000000000000000,
                    0.0000000000000000,
                    -0.23866687325750005,
                    0.0000000000000000,
                    0.0000000000000000,
                ],
            ),
            Data::new(
                Generator::new([0.95590731580453892, 0.95590731580453892]),
                [
                    3.3795617493002125E-002,
                    -7.6439514931065603E-004,
                    7.6439514931064302E-004,
                    0.10076610424769271,
                    -8.4489043732505313E-003,
                ],
            ),
            Data::new(
                Generator::new([0.95590731580453892, 0.40605717473823955]),
                [
                    9.0172578588715732E-002,
                    1.3675773263272826E-002,
                    -1.3675773263272845E-002,
                    -0.32544759695960152,
                    -2.2543144647178933E-002,
                ],
            ),
            Data::new(
                Generator::new([0.68607579756175630, 0.68607579756175630]),
                [
                    0.25150011495314784,
                    -1.5718757184571719E-002,
                    1.5718757184571774E-002,
                    6.2875028738287014E-002,
                    -6.2875028738286959E-002,
                ],
            ),
        ];
        util::assert_check_vec_data_tol(&initial_data, &initial_should_be, tol);
        util::assert_check_vec_data_tol(&final_data, &final_should_be, tol);
    }
}
