// Test the basic Gauss-Kronrod integration rules with a smooth oscillating
// function over an unsymmetric range. This should find any discrepancies in
// the abscissae.
//
// Ported from gsl-2.6/integration/test.c

use rint::quadrature::rule::{
    GaussKronrod15, GaussKronrod21, GaussKronrod31, GaussKronrod41, GaussKronrod51, GaussKronrod61,
};
use rint::quadrature::Basic;

mod util;

basic_test! {
        name: basic_gauss_kronrod_smooth_oscillating_function_unsymmetric_range_15_point,
        function: util::Function3,
        alpha: f64 => 1.3,
        rule: GaussKronrod15,
        lower: 0.3,
        upper: 2.71,
        exp_result: -7.238969575483799046E-01,
        exp_error:  8.760080200939757174E-06,
        exp_result_abs:  1.165564172429140788E+00,
        exp_result_asc:  9.334560307787327371E-01,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-7,
        "Function3 Basic 15-point"
}

basic_test! {
        name: basic_gauss_kronrod_smooth_oscillating_function_unsymmetric_range_21_point,
        function: util::Function3,
        alpha: f64 => 1.3,
        rule: GaussKronrod21,
        lower: 0.3,
        upper: 2.71,
        exp_result: -7.238969575482959717E-01,
        exp_error:  7.999213141433641888E-11,
        exp_result_abs:  1.150829032708484023E+00,
        exp_result_asc:  9.297591249133687619E-01,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-7,
        "Function3 Basic 21-point"
}

basic_test! {
        name: basic_gauss_kronrod_smooth_oscillating_function_unsymmetric_range_31_point,
        function: util::Function3,
        alpha: f64 => 1.3,
        rule: GaussKronrod31,
        lower: 0.3,
        upper: 2.71,
        exp_result: -7.238969575482959717E-01,
        exp_error:  1.285805464427459261E-14,
        exp_result_abs:  1.158150602093290571E+00,
        exp_result_asc:  9.277828092501518853E-01,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-7,
        "Function3 Basic 31-point"
}

basic_test! {
        name: basic_gauss_kronrod_smooth_oscillating_function_unsymmetric_range_41_point,
        function: util::Function3,
        alpha: f64 => 1.3,
        rule: GaussKronrod41,
        lower: 0.3,
        upper: 2.71,
        exp_result: -7.238969575482959717E-01,
        exp_error:  1.286535726271015626E-14,
        exp_result_abs:  1.158808363486595328E+00,
        exp_result_asc:  9.264382258645686985E-01,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-7,
        "Function3 Basic 41-point"
}

basic_test! {
        name: basic_gauss_kronrod_smooth_oscillating_function_unsymmetric_range_51_point,
        function: util::Function3,
        alpha: f64 => 1.3,
        rule: GaussKronrod51,
        lower: 0.3,
        upper: 2.71,
        exp_result: -7.238969575482961938E-01,
        exp_error:  1.285290995039385778E-14,
        exp_result_abs:  1.157687209264406381E+00,
        exp_result_asc:  9.264666884071264263E-01,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-7,
        "Function3 Basic 51-point"
}

basic_test! {
        name: basic_gauss_kronrod_smooth_oscillating_function_unsymmetric_range_61_point,
        function: util::Function3,
        alpha: f64 => 1.3,
        rule: GaussKronrod61,
        lower: 0.3,
        upper: 2.71,
        exp_result: -7.238969575482959717E-01,
        exp_error:  1.286438572027470736E-14,
        exp_result_abs:  1.158720854723590099E+00,
        exp_result_asc:  9.270469641771273972E-01,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-7,
        "Function3 Basic 61-point"
}
