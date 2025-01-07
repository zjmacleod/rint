// Test the basic Gauss-Kronrod integration rules with a smooth positive function.
// Ported from gsl-2.6/integration/test.c

use num_complex::Complex;
use num_complex::ComplexFloat;
use rint::quadrature::Basic;
use rint::quadrature::Rule;
use rint::Limits;

mod util;

basic_test_complex_equal! {
        name: basic_gauss_kronrod_smooth_function_equal_re_im_15_point,
        function: util::ComplexF1F1,
        alpha: f64 => 2.6,
        rule: gk15,
        lower: 0.0,
        upper: 1.0,
        exp_result: 7.716049357767090777E-02,
        exp_error: 2.990224871000550874E-06,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-7,
        "ComplexF1F1 alpha1 = alpha2 Basic 15-point"
}

basic_test_complex_equal! {
        name: basic_gauss_kronrod_smooth_function_equal_re_im_21_point,
        function: util::ComplexF1F1,
        alpha: f64 => 2.6,
        rule: gk21,
        lower: 0.0,
        upper: 1.0,
        exp_result: 7.716049379303084599E-02,
        exp_error: 9.424302194248481445E-08,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-7,
        "ComplexF1F1 alpha1 = alpha2 Basic 21-point"
}

basic_test_complex_equal! {
        name: basic_gauss_kronrod_smooth_function_equal_re_im_31_point,
        function: util::ComplexF1F1,
        alpha: f64 => 2.6,
        rule: gk31,
        lower: 0.0,
        upper: 1.0,
        exp_result: 7.716049382494900855E-02,
        exp_error: 1.713503193600029893E-09,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-7,
        "ComplexF1F1 alpha1 = alpha2 Basic 31-point"
}

basic_test_complex_equal! {
        name: basic_gauss_kronrod_smooth_function_equal_re_im_41_point,
        function: util::ComplexF1F1,
        alpha: f64 => 2.6,
        rule: gk41,
        lower: 0.0,
        upper: 1.0,
        exp_result: 7.716049382681375302E-02,
        exp_error: 9.576386660975511224E-11,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-7,
        "ComplexF1F1 alpha1 = alpha2 Basic 41-point"
}

basic_test_complex_equal! {
        name: basic_gauss_kronrod_smooth_function_equal_re_im_51_point,
        function: util::ComplexF1F1,
        alpha: f64 => 2.6,
        rule: gk51,
        lower: 0.0,
        upper: 1.0,
        exp_result: 7.716049382708510540E-02,
        exp_error: 1.002079980317363772E-11,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-7,
        "ComplexF1F1 alpha1 = alpha2 Basic 51-point"
}

basic_test_complex_equal! {
        name: basic_gauss_kronrod_smooth_function_equal_re_im_61_point,
        function: util::ComplexF1F1,
        alpha: f64 => 2.6,
        rule: gk61,
        lower: 0.0,
        upper: 1.0,
        exp_result: 7.716049382713800753E-02,
        exp_error: 1.566060362296155616E-12,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-5,
        "ComplexF1F1 alpha1 = alpha2 Basic 61-point"
}
