use gauss_kronrod_integration::integration::adaptive::Kind;
use gauss_kronrod_integration::integration::{ErrorBound, GaussKronrodAdaptive};
use gauss_kronrod_integration::rule::{
    GaussKronrod15, GaussKronrod21, GaussKronrod31, GaussKronrod51, GaussKronrod61,
};

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

    let integral =
        GaussKronrodAdaptive::new(lower, upper, error_bound, rule, &function, 1000).unwrap();

    let integral_result = integral.integrate().unwrap();
    let result = integral_result.result();
    let error = integral_result.error();
    let iterations = integral_result.iterations();

    util::test_relative_error(
        result,
        exp_result,
        abs_error_bound,
        "adaptive(f1,15) smooth result",
    )?;
    util::test_relative_error(
        error,
        exp_abserr,
        rel_error_bound,
        "adaptive(f1,15) smooth abserr",
    )?;
    util::test_int(
        iterations,
        exp_iterations,
        "adaptive(f1,15) smooth iterations",
    )?;

    let error_bound = ErrorBound::Relative(1e-10);
    let lower = 1.0;
    let upper = 0.0;

    let function = util::Function1 { alpha };

    let integral =
        GaussKronrodAdaptive::new(lower, upper, error_bound, rule, &function, 1000).unwrap();

    let integral_result = integral.integrate().unwrap();
    let result = integral_result.result();
    let error = integral_result.error();
    let iterations = integral_result.iterations();

    util::test_relative_error(
        result,
        -exp_result,
        abs_error_bound,
        "adaptive(f1,15) smooth result reverse",
    )?;
    util::test_relative_error(
        error,
        exp_abserr,
        rel_error_bound,
        "adaptive(f1,15) smooth abserr reverse",
    )?;
    util::test_int(
        iterations,
        exp_iterations,
        "adaptive(f1,15) smooth iterations reverse",
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

    let integral =
        GaussKronrodAdaptive::new(lower, upper, error_bound, rule, &function, 1000).unwrap();

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

    let error_bound = ErrorBound::Absolute(1e-14);
    let lower = 1.0;
    let upper = 0.0;

    let function = util::Function1 { alpha };

    let integral =
        GaussKronrodAdaptive::new(lower, upper, error_bound, rule, &function, 1000).unwrap();

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

// Test oscillating Function3 with 31 point adaptive integration and absolute error
// bound. Should terminate due to roundoff error.
#[test]
fn test_adaptive_oscillating_should_error_roundoff_31() -> Result<(), String> {
    let exp_result = -7.238969575482959717E-01;
    let exp_abserr = 1.285805464427459261E-14;
    let exp_iterations = 1;

    let error_bound = ErrorBound::Absolute(1e-14);

    let abs_error_bound = 1e-15;
    let rel_error_bound = 1e-6;
    let rule = GaussKronrod31;
    let alpha = 1.3;

    let lower = 0.3;
    let upper = 2.71;

    let function = util::Function3 { alpha };

    let integral =
        GaussKronrodAdaptive::new(lower, upper, error_bound, rule, &function, 1000).unwrap();

    let integral_result = integral.integrate();

    if let Err(err) = integral_result {
        if let Kind::FailedToReachToleranceRoundoff = err.kind() {
            let result = err.result();
            let error = err.error();
            let iterations = err.iterations();

            util::test_relative_error(
                result,
                exp_result,
                abs_error_bound,
                "adaptive(f3,31) smooth result",
            )?;
            util::test_relative_error(
                error,
                exp_abserr,
                rel_error_bound,
                "adaptive(f3,31) smooth abserr",
            )?;
            util::test_int(
                iterations,
                exp_iterations,
                "adaptive(f3,31) smooth iterations",
            )?;
        } else {
            return Err(String::from("wrong error kind"));
        }
    } else {
        return Err(String::from("should have been error"));
    }

    let error_bound = ErrorBound::Absolute(1e-14);
    let lower = 2.71;
    let upper = 0.3;

    let function = util::Function3 { alpha };

    let integral =
        GaussKronrodAdaptive::new(lower, upper, error_bound, rule, &function, 1000).unwrap();

    let integral_result = integral.integrate();

    if let Err(err) = integral_result {
        if let Kind::FailedToReachToleranceRoundoff = err.kind() {
            let result = err.result();
            let error = err.error();
            let iterations = err.iterations();

            util::test_relative_error(
                result,
                -exp_result,
                abs_error_bound,
                "adaptive(f3,31) smooth result reverse",
            )?;
            util::test_relative_error(
                error,
                exp_abserr,
                rel_error_bound,
                "adaptive(f3,31) smooth abserr reverse",
            )?;
            util::test_int(
                iterations,
                exp_iterations,
                "adaptive(f3,31) smooth iterations reverse",
            )?;
        } else {
            return Err(String::from("wrong error kind"));
        }
    } else {
        return Err(String::from("should have been error"));
    }

    Ok(())
}

// Test Function16 with 51 point adaptive integration and absolute error
// bound. Should terminate due to singularity detected (singularity at x=-0.1).
#[test]
fn test_adaptive_singularity_51() -> Result<(), String> {
    let exp_iterations = 51;

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
        if let Kind::PossibleSingularity { lower, upper } = err.kind() {
            let iterations = err.iterations();

            assert!(*lower < -0.1f64);
            assert!(*upper > -0.1f64);

            util::test_int(
                iterations,
                exp_iterations,
                "adaptive(f16,51) smooth iterations",
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
        if let Kind::PossibleSingularity { lower, upper } = err.kind() {
            let iterations = err.iterations();

            assert!(*lower > -0.1f64);
            assert!(*upper < -0.1f64);

            util::test_int(
                iterations,
                exp_iterations,
                "adaptive(f16,51) smooth iterations reverse",
            )?;
        } else {
            return Err(String::from("wrong error kind"));
        }
    } else {
        return Err(String::from("should have been error"));
    }
    Ok(())
}

// TODO: This should return an error, but we successfully integrate,
// Test Function16 with 61 point adaptive integration and absolute error
// bound. Should terminate due to max iterations reached
#[test]
fn test_adaptive_max_iterations_61() -> Result<(), String> {
    let exp_result = 9.565151449233894709;
    let exp_abserr = 1.570369823891028460E+01;
    let exp_iterations = 3;

    let abs_error_bound = 1e-15;
    let rel_error_bound = 1e-6;
    let error_bound = ErrorBound::Absolute(1e-14);

    let rule = GaussKronrod61;
    let alpha = 1;

    let lower = -1.0;
    let upper = 1.0;

    let function = util::Function16 { alpha };

    let integral =
        GaussKronrodAdaptive::new(lower, upper, error_bound, rule, &function, 3).unwrap();

    let integral_result = integral.integrate();

    if let Err(err) = integral_result {
        if let Kind::MaximumIterationsReached = err.kind() {
            let result = err.result();
            let error = err.error();
            let iterations = err.iterations();

            util::test_relative_error(
                result,
                exp_result,
                abs_error_bound,
                "adaptive(f16,61) smooth result",
            )?;
            util::test_relative_error(
                error,
                exp_abserr,
                rel_error_bound,
                "adaptive(f16,61) smooth abserr",
            )?;
            util::test_int(
                iterations,
                exp_iterations,
                "adaptive(f16,61) smooth iterations",
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
        GaussKronrodAdaptive::new(lower, upper, error_bound, rule, &function, 3).unwrap();

    let integral_result = integral.integrate();

    if let Err(err) = integral_result {
        if let Kind::MaximumIterationsReached = err.kind() {
            let result = err.result();
            let error = err.error();
            let iterations = err.iterations();

            util::test_relative_error(
                result,
                -exp_result,
                abs_error_bound,
                "adaptive(f16,61) smooth result negative reverse",
            )?;
            util::test_relative_error(
                error,
                exp_abserr,
                rel_error_bound,
                "adaptive(f16,61) smooth abserr negative reverse",
            )?;
            util::test_int(
                iterations,
                exp_iterations,
                "adaptive(f16,61) smooth iterations negative reverse",
            )?;
        } else {
            return Err(String::from("wrong error kind"));
        }
    } else {
        return Err(String::from("should have been error"));
    }

    Ok(())
}
