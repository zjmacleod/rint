// Test the basic Gauss-Kronrod integration rules with a smooth positive function.
// Ported from gsl-2.6/integration/test.c

use rint::quadrature::Basic;
use rint::quadrature::Rule;
use rint::Limits;

mod util;

basic_test! {
        name: basic_gauss_kronrod_smooth_function_15_point,
        function: util::Function1,
        alpha: f64 => 2.6,
        rule: gk15,
        lower: 0.0,
        upper: 1.0,
        exp_result: 7.716049357767090777E-02,
        exp_error: 2.990224871000550874E-06,
        abs_tolerance: 1.0e-15,
        rel_tolerance: 1.0e-7,
        "Function1 Basic 15-point"
}

basic_test! {
        name: basic_gauss_kronrod_smooth_function_21_point,
        function: util::Function1,
        alpha: f64 => 2.6,
        rule: gk21,
        lower: 0.0,
        upper: 1.0,
        exp_result: 7.716049379303084599E-02,
        exp_error: 9.424302194248481445E-08,
        abs_tolerance: 1.0e-15,
        rel_tolerance: 1.0e-7,
        "Function1 Basic 21-point"
}

basic_test! {
        name: basic_gauss_kronrod_smooth_function_31_point,
        function: util::Function1,
        alpha: f64 => 2.6,
        rule: gk31,
        lower: 0.0,
        upper: 1.0,
        exp_result: 7.716049382494900855E-02,
        exp_error: 1.713503193600029893E-09,
        abs_tolerance: 1.0e-15,
        rel_tolerance: 1.0e-7,
        "Function1 Basic 31-point"
}

basic_test! {
        name: basic_gauss_kronrod_smooth_function_41_point,
        function: util::Function1,
        alpha: f64 => 2.6,
        rule: gk41,
        lower: 0.0,
        upper: 1.0,
        exp_result: 7.716049382681375302E-02,
        exp_error: 9.576386660975511224E-11,
        abs_tolerance: 1.0e-15,
        rel_tolerance: 1.0e-7,
        "Function1 Basic 41-point"
}

basic_test! {
        name: basic_gauss_kronrod_smooth_function_51_point,
        function: util::Function1,
        alpha: f64 => 2.6,
        rule: gk51,
        lower: 0.0,
        upper: 1.0,
        exp_result: 7.716049382708510540E-02,
        exp_error: 1.002079980317363772E-11,
        abs_tolerance: 1.0e-15,
        rel_tolerance: 1.0e-5,
        "Function1 Basic 51-point"
}

basic_test! {
        name: basic_gauss_kronrod_smooth_function_61_point,
        function: util::Function1,
        alpha: f64 => 2.6,
        rule: gk61,
        lower: 0.0,
        upper: 1.0,
        exp_result: 7.716049382713800753E-02,
        exp_error: 1.566060362296155616E-12,
        abs_tolerance: 1.0e-15,
        rel_tolerance: 1.0e-5,
        "Function1 Basic 61-point"
}

// Test the basic Gauss-Kronrod integration rules with a smooth oscillating
// function over an unsymmetric range. This should find any discrepancies in
// the abscissae.
//
// Ported from gsl-2.6/integration/test.c

basic_test! {
        name: basic_gauss_kronrod_smooth_oscillating_function_unsymmetric_range_15_point,
        function: util::Function3,
        alpha: f64 => 1.3,
        rule: gk15,
        lower: 0.3,
        upper: 2.71,
        exp_result: -7.238969575483799046E-01,
        exp_error:  8.760080200939757174E-06,
        abs_tolerance: 1.0e-15,
        rel_tolerance: 1.0e-7,
        "Function3 Basic 15-point"
}

basic_test! {
        name: basic_gauss_kronrod_smooth_oscillating_function_unsymmetric_range_21_point,
        function: util::Function3,
        alpha: f64 => 1.3,
        rule: gk21,
        lower: 0.3,
        upper: 2.71,
        exp_result: -7.238969575482959717E-01,
        exp_error:  7.999213141433641888E-11,
        abs_tolerance: 1.0e-15,
        rel_tolerance: 1.0e-7,
        "Function3 Basic 21-point"
}

basic_test! {
        name: basic_gauss_kronrod_smooth_oscillating_function_unsymmetric_range_31_point,
        function: util::Function3,
        alpha: f64 => 1.3,
        rule: gk31,
        lower: 0.3,
        upper: 2.71,
        exp_result: -7.238969575482959717E-01,
        exp_error:  1.285805464427459261E-14,
        abs_tolerance: 1.0e-15,
        rel_tolerance: 1.0e-7,
        "Function3 Basic 31-point"
}

basic_test! {
        name: basic_gauss_kronrod_smooth_oscillating_function_unsymmetric_range_41_point,
        function: util::Function3,
        alpha: f64 => 1.3,
        rule: gk41,
        lower: 0.3,
        upper: 2.71,
        exp_result: -7.238969575482959717E-01,
        exp_error:  1.286535726271015626E-14,
        abs_tolerance: 1.0e-15,
        rel_tolerance: 1.0e-7,
        "Function3 Basic 41-point"
}

basic_test! {
        name: basic_gauss_kronrod_smooth_oscillating_function_unsymmetric_range_51_point,
        function: util::Function3,
        alpha: f64 => 1.3,
        rule: gk51,
        lower: 0.3,
        upper: 2.71,
        exp_result: -7.238969575482961938E-01,
        exp_error:  1.285290995039385778E-14,
        abs_tolerance: 1.0e-15,
        rel_tolerance: 1.0e-7,
        "Function3 Basic 51-point"
}

basic_test! {
        name: basic_gauss_kronrod_smooth_oscillating_function_unsymmetric_range_61_point,
        function: util::Function3,
        alpha: f64 => 1.3,
        rule: gk61,
        lower: 0.3,
        upper: 2.71,
        exp_result: -7.238969575482959717E-01,
        exp_error:  1.286438572027470736E-14,
        abs_tolerance: 1.0e-15,
        rel_tolerance: 1.0e-7,
        "Function3 Basic 61-point"
}

// Test the basic Gauss-Kronrod integration rules with a positive function that
// has a singularity. This should give large values of the absolute error which
// would find discrepancies in the absolute error calculation.
//
// Ported from gsl-2.6/integration/test.c

basic_test! {
        name: basic_gauss_kronrod_positive_function_with_singularity_15_point,
        function: util::Function1,
        alpha: f64 => -0.9,
        rule: gk15,
        lower: 0.0,
        upper: 1.0,
        exp_result:     1.555688196612745777E+01,
        exp_error:      2.350164577239293706E+01,
        abs_tolerance: 1.0e-15,
        rel_tolerance: 1.0e-7,
        "Function1 with singularity Basic 15-point"
}

basic_test! {
        name: basic_gauss_kronrod_positive_function_with_singularity_21_point,
        function: util::Function1,
        alpha: f64 => -0.9,
        rule: gk21,
        lower: 0.0,
        upper: 1.0,
        exp_result:     1.799045317938126232E+01,
        exp_error:      2.782360287710622515E+01,
        abs_tolerance: 1.0e-15,
        rel_tolerance: 1.0e-7,
        "Function1 with singularity Basic 21-point"
}

basic_test! {
        name: basic_gauss_kronrod_positive_function_with_singularity_31_point,
        function: util::Function1,
        alpha: f64 => -0.9,
        rule: gk31,
        lower: 0.0,
        upper: 1.0,
        exp_result:     2.081873305159121657E+01,
        exp_error:      3.296500137482590276E+01,
        abs_tolerance: 1.0e-15,
        rel_tolerance: 1.0e-7,
        "Function1 with singularity Basic 31-point"
}

basic_test! {
        name: basic_gauss_kronrod_positive_function_with_singularity_41_point,
        function: util::Function1,
        alpha: f64 => -0.9,
        rule: gk41,
        lower: 0.0,
        upper: 1.0,
        exp_result:     2.288677623903126701E+01,
        exp_error:      3.671538820274916048E+01,
        abs_tolerance: 1.0e-15,
        rel_tolerance: 1.0e-7,
        "Function1 with singularity Basic 41-point"
}

basic_test! {
        name: basic_gauss_kronrod_positive_function_with_singularity_51_point,
        function: util::Function1,
        alpha: f64 => -0.9,
        rule: gk51,
        lower: 0.0,
        upper: 1.0,
        exp_result:     2.449953612016972215E+01,
        exp_error:      3.967771249391228849E+01,
        abs_tolerance: 1.0e-15,
        rel_tolerance: 1.0e-7,
        "Function1 with singularity Basic 51-point"
}

basic_test! {
        name: basic_gauss_kronrod_positive_function_with_singularity_61_point,
        function: util::Function1,
        alpha: f64 => -0.9,
        rule: gk61,
        lower: 0.0,
        upper: 1.0,
        exp_result:     2.583030240976628988E+01,
        exp_error:      4.213750493076978643E+01,
        abs_tolerance: 1.0e-15,
        rel_tolerance: 1.0e-7,
        "Function1 with singularity Basic 61-point"
}
