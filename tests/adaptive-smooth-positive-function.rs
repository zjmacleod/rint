use gauss_kronrod_integration::integration::{ErrorBound, GaussKronrodAdaptive};
use gauss_kronrod_integration::rule::{GaussKronrod15, GaussKronrod21};

mod util;

// Test the smooth Function1 with 15 point adaptive integration and relative error
// bound.
#[test]
fn test_adaptive_smooth_relative_error_15() -> Result<(), String> {
    let exp_result = 7.716049382715854665E-02;
    let exp_abserr = 6.679384885865053037E-12;
    let exp_iterations = 6;

    let rel_error_bound = 1e-6;
    let abs_error_bound = 1e-15;
    let rule = GaussKronrod15;
    let error_bound = ErrorBound::Relative(1e-10);
    let alpha = 2.6;

    let lower = 0.0;
    let upper = 1.0;

    let function = util::Function1 { alpha };

    let integral = GaussKronrodAdaptive::new(
        lower,
        upper,
        error_bound,
        rule,
        &function,
        1000,
    )
    .unwrap();

    let integral_result = integral.integrate().unwrap();
    let result = integral_result.result();
    let error = integral_result.error();
    let iterations = integral_result.iterations();

    util::test_relative_error(
        result,
        exp_result,
        abs_error_bound,
        "qk15(f1) smooth result",
    )?;
    util::test_relative_error(
        error,
        exp_abserr,
        rel_error_bound,
        "qk15(f1) smooth abserr",
    )?;
    util::test_int(
        iterations,
        exp_iterations,
        "adaptive qk15(f1) smooth iterations",
    )?;

    let error_bound = ErrorBound::Relative(1e-10);
    let lower = 1.0;
    let upper = 0.0;

    let function = util::Function1 { alpha };

    let integral = GaussKronrodAdaptive::new(
        lower,
        upper,
        error_bound,
        rule,
        &function,
        1000,
    )
    .unwrap();

    let integral_result = integral.integrate().unwrap();
    let result = integral_result.result();
    let error = integral_result.error();
    let iterations = integral_result.iterations();

    util::test_relative_error(
        result,
        -exp_result,
        abs_error_bound,
        "qk15(f1) smooth result",
    )?;
    util::test_relative_error(
        error,
        exp_abserr,
        rel_error_bound,
        "qk15(f1) smooth abserr",
    )?;
    util::test_int(
        iterations,
        exp_iterations,
        "adaptive qk15(f1) smooth iterations",
    )?;

    Ok(())
}

// Test the smooth Function1 with 21 point adaptive integration and absolute error
// bound.
#[test]
fn test_adaptive_smooth_absolute_error_21() -> Result<(), String> {
    let exp_result = 7.716049382716050342E-02;
    let exp_abserr = 2.227969521869139532E-15;
    let exp_iterations = 8;

    let error_bound = ErrorBound::Absolute(1e-14);

    let abs_error_bound = 1e-15;
    let rel_error_bound = 1e-6;
    let rule = GaussKronrod21;
    let alpha = 2.6;

    let lower = 0.0;
    let upper = 1.0;

    let function = util::Function1 { alpha };

    let integral = GaussKronrodAdaptive::new(
        lower,
        upper,
        error_bound,
        rule,
        &function,
        1000,
    )
    .unwrap();

    let integral_result = integral.integrate().unwrap();
    let result = integral_result.result();
    let error = integral_result.error();
    let iterations = integral_result.iterations();

    util::test_relative_error(
        result,
        exp_result,
        abs_error_bound,
        "qk15(f1) smooth result",
    )?;
    util::test_relative_error(
        error,
        exp_abserr,
        rel_error_bound,
        "qk15(f1) smooth abserr",
    )?;
    util::test_int(
        iterations,
        exp_iterations,
        "adaptive qk15(f1) smooth iterations",
    )?;

    let error_bound = ErrorBound::Absolute(1e-14);
    let lower = 1.0;
    let upper = 0.0;

    let function = util::Function1 { alpha };

    let integral = GaussKronrodAdaptive::new(
        lower,
        upper,
        error_bound,
        rule,
        &function,
        1000,
    )
    .unwrap();

    let integral_result = integral.integrate().unwrap();
    let result = integral_result.result();
    let error = integral_result.error();
    let iterations = integral_result.iterations();

    util::test_relative_error(
        result,
        -exp_result,
        abs_error_bound,
        "qk15(f1) smooth result",
    )?;
    util::test_relative_error(
        error,
        exp_abserr,
        rel_error_bound,
        "qk15(f1) smooth abserr",
    )?;
    util::test_int(
        iterations,
        exp_iterations,
        "adaptive qk15(f1) smooth iterations",
    )?;

    Ok(())
}
