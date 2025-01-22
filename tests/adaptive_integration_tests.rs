use rint::quadrature::Adaptive;
use rint::quadrature::Rule;
use rint::Limits;
use rint::Tolerance;

mod util;

// Test the smooth Function1 with 15 point adaptive integration and relative error
// bound.
adaptive_test_passing! {
        name: adaptive_gauss_kronrod_smooth_positive_function_relative_error_bound_15_point,
        function: util::Function1,
        alpha: f64 => 2.6,
        rule: gk15,
        lower: 0.0,
        upper: 1.0,
        exp_result:      7.716049382715854665E-02,
        exp_error:       6.679384885865053037E-12,
        exp_iterations:  6,
        exp_evaluations: 165,
        tolerance_rel: 1e-10,
        test_abs_error_bound: 1e-15,
        test_rel_error_bound: 1e-6,
        "Function1 Adaptive relative bound 15-point"
}

// Test the smooth Function1 with 21 point adaptive integration and absolute error
// bound.
adaptive_test_passing! {
        name: adaptive_gauss_kronrod_smooth_positive_function_absolute_error_bound_21_point,
        function: util::Function1,
        alpha: f64 => 2.6,
        rule: gk21,
        lower: 0.0,
        upper: 1.0,
        exp_result:      7.716049382716050342E-02,
        exp_error:       2.227969521869139532E-15,
        exp_iterations:  8,
        exp_evaluations: 315,
        tolerance_abs: 1e-14,
        test_abs_error_bound: 1e-15,
        test_rel_error_bound: 1e-6,
        "Function1 Adaptive absolute bound 21-point"
}
