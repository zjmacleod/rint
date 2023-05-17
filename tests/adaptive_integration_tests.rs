use rint::integration::{ErrorBound, GaussKronrodAdaptive, Kind};
use rint::rule::{GaussKronrod15, GaussKronrod21, GaussKronrod31, GaussKronrod51, GaussKronrod61};

mod util;

// Test the smooth Function1 with 15 point adaptive integration and relative error
// bound.
adaptive_test_passing! {
        name: adaptive_gauss_kronrod_smooth_positive_function_relative_error_bound_15_point,
        function: util::Function1,
        alpha: f64 => 2.6,
        rule: GaussKronrod15,
        lower: 0.0,
        upper: 1.0,
        exp_result:      7.716049382715854665E-02,
        exp_error:       6.679384885865053037E-12,
        exp_iterations:  6,
        exp_evaluations: 165,
        tolerance_rel: 1e-10,
        test_abs_error_bound: 1e-15,
        test_rel_error_bound: 1e-6,
        "Function1 GaussKronrodAdaptive relative bound 15-point"
}

// Test the smooth Function1 with 21 point adaptive integration and absolute error
// bound.
adaptive_test_passing! {
        name: adaptive_gauss_kronrod_smooth_positive_function_absolute_error_bound_21_point,
        function: util::Function1,
        alpha: f64 => 2.6,
        rule: GaussKronrod21,
        lower: 0.0,
        upper: 1.0,
        exp_result:      7.716049382716050342E-02,
        exp_error:       2.227969521869139532E-15,
        exp_iterations:  8,
        exp_evaluations: 315,
        tolerance_abs: 1e-14,
        test_abs_error_bound: 1e-15,
        test_rel_error_bound: 1e-6,
        "Function1 GaussKronrodAdaptive absolute bound 21-point"
}

// Test oscillating Function3 with 31 point adaptive integration and absolute error
// bound. Should terminate due to roundoff error.
adaptive_test_error! {
        name: adaptive_gauss_kronrod_terminates_due_to_roundoff_error_31_point,
        function: util::Function3,
        alpha: f64 => 1.3,
        rule: GaussKronrod31,
        lower: 0.3,
        upper: 2.71,
        exp_result:      -7.238969575482959717E-01,
        exp_error:        1.285805464427459261E-14,
        exp_iterations:  1,
        exp_evaluations: 31,
        tolerance_abs: 1e-14,
        test_abs_error_bound: 1e-15,
        test_rel_error_bound: 1e-6,
        kind: RoundoffErrorDetected,
        "Function3 GaussKronrodAdaptive absolute bound 31-point (terminates due to roundoff error)"
}

// Test Function16 with 61 point adaptive integration and absolute error
// bound. Should terminate due to max iterations reached
adaptive_test_error! {
        name: adaptive_gauss_kronrod_terminates_due_to_max_iterations_reached_61_point,
        function: util::Function16,
        alpha: i32 => 1,
        rule: GaussKronrod61,
        lower: -1.0,
        upper:  1.0,
        iterations: 3,
        exp_result:      9.565151449233894709,
        exp_error:       1.570369823891028460E+01,
        exp_iterations:  3,
        exp_evaluations: 305,
        tolerance_abs: 1e-14,
        test_abs_error_bound: 1e-15,
        test_rel_error_bound: 1e-6,
        kind: MaximumIterationsReached,
        "Function16 GaussKronrodAdaptive absolute bound 61-point (terminates due to max iterations reached)"
}

// Test Function16 with 51 point adaptive integration and absolute error
// bound. Should terminate due to singularity detected (singularity at x=-0.1).
#[test]
fn test_adaptive_singularity_51() -> Result<(), String> {
    let exp_iterations = 51;
    let exp_evaluations = 5151;

    let error_bound = ErrorBound::Absolute(1e-14);

    let rule = GaussKronrod51;
    let alpha = 2;

    let lower = -1.0;
    let upper = 1.0;

    let function = util::Function16 { alpha };

    let integral =
        GaussKronrodAdaptive::new(lower, upper, error_bound, rule, &function, 1000).unwrap();

    let integral_result = integral.integrate();

    if let Err(err) = integral_result {
        if let Kind::BadIntegrandBehaviour { lower, upper } = err.kind() {
            let iterations = err.iterations();
            let evaluations = err.function_evaluations();

            assert!(lower < -0.1f64);
            assert!(upper > -0.1f64);

            util::test_int(
                iterations,
                exp_iterations,
                "adaptive(f16,51) smooth iterations",
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                "adaptive(f16,51) smooth evaluations",
            )?;
        } else {
            return Err(String::from("wrong error kind"));
        }
    } else {
        return Err(String::from("should have been error"));
    }

    let error_bound = ErrorBound::Absolute(1e-14);
    let lower = 1.0;
    let upper = -1.0;

    let function = util::Function16 { alpha };

    let integral =
        GaussKronrodAdaptive::new(lower, upper, error_bound, rule, &function, 1000).unwrap();

    let integral_result = integral.integrate();

    if let Err(err) = integral_result {
        if let Kind::BadIntegrandBehaviour { lower, upper } = err.kind() {
            let iterations = err.iterations();
            let evaluations = err.function_evaluations();

            assert!(lower > -0.1f64);
            assert!(upper < -0.1f64);

            util::test_int(
                iterations,
                exp_iterations,
                "adaptive(f16,51) smooth iterations reverse",
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                "adaptive(f16,51) smooth evaluations reverse",
            )?;
        } else {
            return Err(String::from("wrong error kind"));
        }
    } else {
        return Err(String::from("should have been error"));
    }
    Ok(())
}
