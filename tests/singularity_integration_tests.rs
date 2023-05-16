//use rint::integration::adaptive::Kind;
use rint::integration::{ErrorBound, GaussKronrodAdaptiveSingularity};

mod util;

// Test the smooth Function1 with 15 point adaptive integration and relative error
// bound.
#[test]
fn test_singularity_smooth_relative_error_21() -> Result<(), String> {
    let exp_result = 7.716049382715789440E-02;
    let exp_abserr = 2.216394961010438404E-12;
    let exp_iterations = 5;
    let exp_evaluations = 189;

    let rel_error_bound = 1e-6;
    let abs_error_bound = 1e-15;
    let error_bound = ErrorBound::Relative(1e-10);
    let alpha = 2.6;

    let lower = 0.0;
    let upper = 1.0;

    let function = util::Function1 { alpha };

    let integral =
        GaussKronrodAdaptiveSingularity::general(lower, upper, error_bound, &function, 1000)
            .unwrap();

    let integral_result = integral.integrate().unwrap();
    let result = integral_result.result();
    let error = integral_result.error();
    let iterations = integral_result.iterations();
    let evaluations = integral_result.function_evaluations();

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
    util::test_int(
        evaluations,
        exp_evaluations,
        "adaptive(f1,21) smooth evaluations",
    )?;

    let lower = 1.0;
    let upper = 0.0;
    let error_bound = ErrorBound::Relative(1e-10);

    let integral =
        GaussKronrodAdaptiveSingularity::general(lower, upper, error_bound, &function, 1000)
            .unwrap();

    let integral_result = integral.integrate().unwrap();
    let result = integral_result.result();
    let error = integral_result.error();
    let iterations = integral_result.iterations();
    let evaluations = integral_result.function_evaluations();

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
    util::test_int(
        evaluations,
        exp_evaluations,
        "adaptive(f1,21) smooth evaluations reverse",
    )?;

    Ok(())
}

#[test]
fn test_singularity_f11_absolute_error_21() -> Result<(), String> {
    // XXX Test passes but result for error estimate 2e-14 off from GSL.
    let exp_result = -5.908755278982136588E+03;
    let exp_abserr = 1.299646281053874554E-10;
    let exp_iterations = 9;
    let exp_evaluations = 357;

    let rel_error_bound = 1e-3;
    let abs_error_bound = 1e-15;
    let error_bound = ErrorBound::Absolute(1e-7);
    let alpha = 2.0;

    let lower = 1.0;
    let upper = 1000.0;

    let function = util::Function11 { alpha };

    let integral =
        GaussKronrodAdaptiveSingularity::general(lower, upper, error_bound, &function, 1000)
            .unwrap();

    let integral_result = integral.integrate().unwrap();
    let result = integral_result.result();
    let error = integral_result.error();
    let iterations = integral_result.iterations();
    let evaluations = integral_result.function_evaluations();

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
    util::test_int(
        evaluations,
        exp_evaluations,
        "adaptive(f11,21) f11 evaluations",
    )?;

    let lower = 1000.0;
    let upper = 1.0;
    let error_bound = ErrorBound::Absolute(1e-7);

    let integral =
        GaussKronrodAdaptiveSingularity::general(lower, upper, error_bound, &function, 1000)
            .unwrap();

    let integral_result = integral.integrate().unwrap();
    let result = integral_result.result();
    let error = integral_result.error();
    let iterations = integral_result.iterations();
    let evaluations = integral_result.function_evaluations();

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
    util::test_int(
        evaluations,
        exp_evaluations,
        "adaptive(f11,21) f11 evaluations reverse",
    )?;

    Ok(())
}

/* Test infinite range integral f455 using a relative error bound */

#[test]
fn test_singularity_f455_relative_error_infinite() -> Result<(), String> {
    // XXX Test passes but result for error estimate 2e-14 off from GSL.
    let exp_result = -3.616892186127022568E-01;
    let exp_abserr = 3.016716913328831851E-06;
    let exp_iterations = 10;
    let exp_evaluations = 285;

    let rel_error_bound = 1e-5;
    let abs_error_bound = 1e-14;
    let error_bound = ErrorBound::Relative(1e-3);

    let function = util::Function455;

    let integral =
        GaussKronrodAdaptiveSingularity::semi_infinite_positive(0.0, error_bound, &function, 1000)
            .unwrap();

    let integral_result = integral.integrate().unwrap();
    let result = integral_result.result();
    let error = integral_result.error();
    let iterations = integral_result.iterations();
    let evaluations = integral_result.function_evaluations();

    util::test_relative_error(
        result,
        exp_result,
        abs_error_bound,
        "adaptive(f455,15) f455 result",
    )?;
    util::test_relative_error(
        error,
        exp_abserr,
        rel_error_bound,
        "adaptive(f455,15) f455 abserr",
    )?;
    util::test_int(
        iterations,
        exp_iterations,
        "adaptive(f455,15) f455 iterations",
    )?;
    util::test_int(
        evaluations,
        exp_evaluations,
        "adaptive(f455,15) f455 evaluations",
    )?;

    Ok(())
}

/* Test infinite range integral f15 using a relative error bound */
#[test]
fn test_singularity_f15_relative_error_infinite() -> Result<(), String> {
    // XXX Test passes but result for error estimate 2e-14 off from GSL.
    let exp_result = 6.553600000000024738E+04;
    let exp_abserr = 7.121667111456009280E-04;
    let exp_iterations = 10;
    let exp_evaluations = 285;

    let rel_error_bound = 1e-5;
    let abs_error_bound = 1e-14;
    let error_bound = ErrorBound::Relative(1e-7);

    let alpha = 5;

    let function = util::Function15 { alpha };

    let integral =
        GaussKronrodAdaptiveSingularity::semi_infinite_positive(0.0, error_bound, &function, 1000)
            .unwrap();

    let integral_result = integral.integrate().unwrap();
    let result = integral_result.result();
    let error = integral_result.error();
    let iterations = integral_result.iterations();
    let evaluations = integral_result.function_evaluations();

    util::test_relative_error(
        result,
        exp_result,
        abs_error_bound,
        "adaptive(f15,15) f15 result",
    )?;
    util::test_relative_error(
        error,
        exp_abserr,
        rel_error_bound,
        "adaptive(f15,15) f15 abserr",
    )?;
    util::test_int(
        iterations,
        exp_iterations,
        "adaptive(f15,15) f15 iterations",
    )?;
    util::test_int(
        evaluations,
        exp_evaluations,
        "adaptive(f15,15) f15 evaluations",
    )?;

    Ok(())
}

/* Test infinite range integral f16 using an absolute error bound */
#[test]
fn test_singularity_f16_absolute_error_infinite() -> Result<(), String> {
    let exp_result = 1.000000000006713292E-04;
    let exp_abserr = 3.084062020905636316E-09;
    let exp_iterations = 6;
    let exp_evaluations = 165;

    let rel_error_bound = 1e-5;
    let abs_error_bound = 1e-14;
    let error_bound = ErrorBound::Absolute(1e-7);

    let alpha = 1;

    let function = util::Function16 { alpha };

    let integral =
        GaussKronrodAdaptiveSingularity::semi_infinite_positive(99.9, error_bound, &function, 1000)
            .unwrap();

    let integral_result = integral.integrate().unwrap();
    let result = integral_result.result();
    let error = integral_result.error();
    let iterations = integral_result.iterations();
    let evaluations = integral_result.function_evaluations();

    util::test_relative_error(
        result,
        exp_result,
        abs_error_bound,
        "adaptive(f16,15) f16 result",
    )?;
    util::test_relative_error(
        error,
        exp_abserr,
        rel_error_bound,
        "adaptive(f16,15) f16 abserr",
    )?;
    util::test_int(
        iterations,
        exp_iterations,
        "adaptive(f16,15) f16 iterations",
    )?;
    util::test_int(
        evaluations,
        exp_evaluations,
        "adaptive(f16,15) f16 evaluations",
    )?;

    Ok(())
}

/* Test infinite range integral myfn1 using a absolute error bound */
#[test]
fn test_singularity_myfn1_absolute_error_infinite() -> Result<(), String> {
    let exp_result = 2.275875794468747770E+00;
    let exp_abserr = 7.436490118267390744E-09;
    let exp_iterations = 5;
    let exp_evaluations = 135;

    let rel_error_bound = 1e-5;
    let abs_error_bound = 1e-14;
    let error_bound = ErrorBound::Absolute(1e-7);

    let function = util::MyFunciton1;

    let integral = GaussKronrodAdaptiveSingularity::infinite(error_bound, &function, 1000).unwrap();

    let integral_result = integral.integrate().unwrap();
    let result = integral_result.result();
    let error = integral_result.error();
    let iterations = integral_result.iterations();
    let evaluations = integral_result.function_evaluations();

    util::test_relative_error(
        result,
        exp_result,
        abs_error_bound,
        "adaptive(myfn1,15) myfn1 result",
    )?;
    util::test_relative_error(
        error,
        exp_abserr,
        rel_error_bound,
        "adaptive(myfn1,15) myfn1 abserr",
    )?;
    util::test_int(
        iterations,
        exp_iterations,
        "adaptive(myfn1,15) myfn1 iterations",
    )?;
    util::test_int(
        evaluations,
        exp_evaluations,
        "adaptive(myfn1,15) myfn1 evaluations",
    )?;

    Ok(())
}

/* Test infinite range integral myfn1 using a absolute error bound lower */
#[test]
fn test_singularity_myfn1_absolute_error_semi_infinite() -> Result<(), String> {
    let exp_result = 2.718281828459044647E+00;
    let exp_abserr = 1.588185109253204805E-10;
    let exp_iterations = 5;
    let exp_evaluations = 135;

    let rel_error_bound = 1e-5;
    let abs_error_bound = 1e-14;
    let error_bound = ErrorBound::Absolute(1e-7);

    let alpha = 1.0;

    let function = util::MyFunciton2 { alpha };

    let integral =
        GaussKronrodAdaptiveSingularity::semi_infinite_negative(1.0, error_bound, &function, 1000)
            .unwrap();

    let integral_result = integral.integrate().unwrap();
    let result = integral_result.result();
    let error = integral_result.error();
    let iterations = integral_result.iterations();
    let evaluations = integral_result.function_evaluations();

    util::test_relative_error(
        result,
        exp_result,
        abs_error_bound,
        "adaptive(myfn1,15) myfn1 result",
    )?;
    util::test_relative_error(
        error,
        exp_abserr,
        rel_error_bound,
        "adaptive(myfn1,15) myfn1 abserr",
    )?;
    util::test_int(
        iterations,
        exp_iterations,
        "adaptive(myfn1,15) myfn1 iterations",
    )?;
    util::test_int(
        evaluations,
        exp_evaluations,
        "adaptive(myfn1,15) myfn1 evaluations",
    )?;

    Ok(())
}
