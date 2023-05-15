use gauss_kronrod_integration::integration::adaptive::Kind;
use gauss_kronrod_integration::integration::{ErrorBound, GaussKronrodAdaptiveSingularity};
use gauss_kronrod_integration::rule::{
    GaussKronrod15, GaussKronrod21, GaussKronrod31, GaussKronrod51, GaussKronrod61,
};

mod util;

// Test the smooth Function1 with 15 point adaptive integration and relative error
// bound.
#[test]
fn test_singularity_smooth_relative_error_21() -> Result<(), String> {
    let exp_result = 7.716049382715789440E-02;
    let exp_abserr = 2.216394961010438404E-12;
    let exp_iterations = 5;

    let rel_error_bound = 1e-6;
    let abs_error_bound = 1e-15;
    let rule = GaussKronrod21;
    let error_bound = ErrorBound::Relative(1e-10);
    let alpha = 2.6;

    let lower = 0.0;
    let upper = 1.0;

    let function = util::Function1 { alpha };

    let integral =
        GaussKronrodAdaptiveSingularity::new(lower, upper, error_bound, rule, &function, 1000)
            .unwrap();

    let integral_result = integral.integrate().unwrap();
    let result = integral_result.result();
    let error = integral_result.error();
    let iterations = integral_result.iterations();

    util::test_relative_error(
        result,
        exp_result,
        abs_error_bound,
        "adaptive(f1,21) smooth result",
    )?;
    util::test_relative_error(
        error,
        exp_abserr,
        rel_error_bound,
        "adaptive(f1,21) smooth abserr",
    )?;
    util::test_int(
        iterations,
        exp_iterations,
        "adaptive(f1,21) smooth iterations",
    )?;

    let error_bound = ErrorBound::Relative(1e-10);
    let lower = 1.0;
    let upper = 0.0;

    let function = util::Function1 { alpha };

    let integral =
        GaussKronrodAdaptiveSingularity::new(lower, upper, error_bound, rule, &function, 1000)
            .unwrap();

    let integral_result = integral.integrate().unwrap();
    let result = integral_result.result();
    let error = integral_result.error();
    let iterations = integral_result.iterations();

    util::test_relative_error(
        result,
        -exp_result,
        abs_error_bound,
        "adaptive(f1,21) smooth result reverse",
    )?;
    util::test_relative_error(
        error,
        exp_abserr,
        rel_error_bound,
        "adaptive(f1,21) smooth abserr reverse",
    )?;
    util::test_int(
        iterations,
        exp_iterations,
        "adaptive(f1,21) smooth iterations reverse",
    )?;

    Ok(())
}

#[test]
fn test_singularity_f11_absolute_error_21() -> Result<(), String> {
    // XXX Test passes but result for error estimate 2e-14 off from GSL.
    let exp_result = -5.908755278982136588E+03;
    let exp_abserr = 1.299646281053874554E-10;
    let exp_iterations = 9;

    let rel_error_bound = 1e-3;
    let abs_error_bound = 1e-15;
    let rule = GaussKronrod21;
    let error_bound = ErrorBound::Absolute(1e-7);
    let alpha = 2.0;

    let lower = 1.0;
    let upper = 1000.0;

    let function = util::Function11 { alpha };

    let integral =
        GaussKronrodAdaptiveSingularity::new(lower, upper, error_bound, rule, &function, 1000)
            .unwrap();

    let integral_result = integral.integrate().unwrap();
    let result = integral_result.result();
    let error = integral_result.error();
    let iterations = integral_result.iterations();

    util::test_relative_error(
        result,
        exp_result,
        abs_error_bound,
        "adaptive(f11,21) f11 result",
    )?;
    util::test_relative_error(
        error,
        exp_abserr,
        rel_error_bound,
        "adaptive(f11,21) f11 abserr",
    )?;
    util::test_int(
        iterations,
        exp_iterations,
        "adaptive(f11,21) f11 iterations",
    )?;

    let error_bound = ErrorBound::Absolute(1e-7);
    let alpha = 2.0;

    let lower = 1000.0;
    let upper = 1.0;

    let function = util::Function11 { alpha };

    let integral =
        GaussKronrodAdaptiveSingularity::new(lower, upper, error_bound, rule, &function, 1000)
            .unwrap();

    let integral_result = integral.integrate().unwrap();
    let result = integral_result.result();
    let error = integral_result.error();
    let iterations = integral_result.iterations();

    util::test_relative_error(
        result,
        -exp_result,
        abs_error_bound,
        "adaptive(f11,21) f11 result reverse",
    )?;
    util::test_relative_error(
        error,
        exp_abserr,
        rel_error_bound,
        "adaptive(f11,21) f11 abserr reverse",
    )?;
    util::test_int(
        iterations,
        exp_iterations,
        "adaptive(f11,21) f11 iterations reverse",
    )?;

    Ok(())
}
