// Test the basic Gauss-Kronrod integration rules with a smooth positive function.
// Ported from gsl-2.6/integration/test.c

use crate::quadrature::rule::{
    GaussKronrod15, GaussKronrod21, GaussKronrod31, GaussKronrod41, GaussKronrod51, GaussKronrod61,
};
use crate::quadrature::Basic;
use crate::Limits;

basic_test! {
        name: basic_gauss_kronrod_smooth_function_15_point,
        function: util::Function1,
        alpha: f64 => 2.6,
        rule: GaussKronrod15,
        lower: 0.0,
        upper: 1.0,
        exp_result: 7.716049357767090777E-02,
        exp_error: 2.990224871000550874E-06,
        exp_result_abs: 7.716049357767090777E-02,
        exp_result_asc: 4.434273814139995384E-02,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-7,
        "Function1 Basic 15-point"
}

basic_test! {
        name: basic_gauss_kronrod_smooth_function_21_point,
        function: util::Function1,
        alpha: f64 => 2.6,
        rule: GaussKronrod21,
        lower: 0.0,
        upper: 1.0,
        exp_result: 7.716049379303084599E-02,
        exp_error: 9.424302194248481445E-08,
        exp_result_abs: 7.716049379303084599E-02,
        exp_result_asc: 4.434311425038358484E-02,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-7,
        "Function1 Basic 21-point"
}

basic_test! {
        name: basic_gauss_kronrod_smooth_function_31_point,
        function: util::Function1,
        alpha: f64 => 2.6,
        rule: GaussKronrod31,
        lower: 0.0,
        upper: 1.0,
        exp_result: 7.716049382494900855E-02,
        exp_error: 1.713503193600029893E-09,
        exp_result_abs: 7.716049382494900855E-02,
        exp_result_asc: 4.427995051868838933E-02,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-7,
        "Function1 Basic 31-point"
}

basic_test! {
        name: basic_gauss_kronrod_smooth_function_41_point,
        function: util::Function1,
        alpha: f64 => 2.6,
        rule: GaussKronrod41,
        lower: 0.0,
        upper: 1.0,
        exp_result: 7.716049382681375302E-02,
        exp_error: 9.576386660975511224E-11,
        exp_result_abs: 7.716049382681375302E-02,
        exp_result_asc: 4.421521169637691873E-02,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-7,
        "Function1 Basic 41-point"
}

basic_test! {
        name: basic_gauss_kronrod_smooth_function_51_point,
        function: util::Function1,
        alpha: f64 => 2.6,
        rule: GaussKronrod51,
        lower: 0.0,
        upper: 1.0,
        exp_result: 7.716049382708510540E-02,
        exp_error: 1.002079980317363772E-11,
        exp_result_abs: 7.716049382708510540E-02,
        exp_result_asc: 4.416474291216854892E-02,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-5,
        "Function1 Basic 51-point"
}

basic_test! {
        name: basic_gauss_kronrod_smooth_function_61_point,
        function: util::Function1,
        alpha: f64 => 2.6,
        rule: GaussKronrod61,
        lower: 0.0,
        upper: 1.0,
        exp_result: 7.716049382713800753E-02,
        exp_error: 1.566060362296155616E-12,
        exp_result_abs: 7.716049382713800753E-02,
        exp_result_asc: 4.419287685934316506E-02,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-5,
        "Function1 Basic 61-point"
}
