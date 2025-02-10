use crate::quadrature::Adaptive;
use crate::quadrature::Rule;
use crate::IntegrationErrorKind;
use crate::Limits;
use crate::Tolerance;

use crate::quadrature::tests::util;

#[test]
fn adaptive_gauss_kronrod_terminates_due_to_roundoff_error_31_point() -> Result<(), String> {
    let description =
        "Function3 Adaptive absolute bound 31-point (terminates due to roundoff error)";
    let exp_result = -7.238969575482959717E-01;
    let exp_error = 1.285805464427459261E-14;
    let exp_iterations = 1;
    let exp_evaluations = 31;

    let abs_error_bound = 1e-15;
    let rel_error_bound = 1e-6;
    let rule = Rule::gk31();
    let error_bound = Tolerance::Absolute(1e-14);
    let alpha: f64 = 1.3;

    let lower = 0.3;
    let upper = 2.71;

    let function = util::Function3::new(alpha);

    let integral = Adaptive::new(
        &function,
        &rule,
        Limits::new(lower, upper),
        error_bound,
        1000,
    )
    .unwrap();

    let integral_result = integral.integrate();

    if let Err(err) = integral_result {
        if let IntegrationErrorKind::RoundoffErrorDetected = err.kind() {
            let estimate = err.estimate();
            let result = estimate.result();
            let error = estimate.error();
            let iterations = estimate.iterations();
            let evaluations = estimate.evaluations();

            util::test_relative_error(
                result,
                exp_result,
                abs_error_bound,
                &format!("{} result", description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_error_bound,
                &format!("{} error", description),
            )?;
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations", description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations", description),
            )?;
        } else {
            return Err(format!(
                "Wrong error type. Expected {} but recieved {:?}",
                stringify!(RoundoffErrorDetected),
                err.kind()
            ));
        }
    } else {
        return Err(String::from(
            "Integration returned Ok, but should have been an error.",
        ));
    }

    let rule = Rule::gk31();
    let error_bound = Tolerance::Absolute(1e-14);
    let lower = 2.71;
    let upper = 0.3;

    let integral = Adaptive::new(
        &function,
        &rule,
        Limits::new(lower, upper),
        error_bound,
        1000,
    )
    .unwrap();

    let integral_result = integral.integrate();

    if let Err(err) = integral_result {
        if let IntegrationErrorKind::RoundoffErrorDetected = err.kind() {
            let estimate = err.estimate();
            let result = estimate.result();
            let error = estimate.error();
            let iterations = estimate.iterations();
            let evaluations = estimate.evaluations();

            util::test_relative_error(
                result,
                -exp_result,
                abs_error_bound,
                &format!("{} result reverse", description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_error_bound,
                &format!("{} error reverse", description),
            )?;
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations reverse", description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations reverse", description),
            )?;
        } else {
            return Err(format!(
                "Wrong error type. Expected {} but recieved {:?}",
                stringify!(RoundoffErrorDetected),
                err.kind()
            ));
        }
    } else {
        return Err(String::from(
            "Integration returned Ok, but should have been an error.",
        ));
    }

    Ok(())
}

#[test]
fn adaptive_gauss_kronrod_terminates_due_to_max_iterations_reached_61_point() -> Result<(), String>
{
    let description =
        "Function16 Adaptive absolute bound 61-point (terminates due to max iterations reached)";
    let exp_result = 9.565151449233894709;
    let exp_error = 1.570369823891028460E+01;
    let exp_iterations = 3;
    let exp_evaluations = 305;

    let abs_error_bound = 1e-15;
    let rel_error_bound = 1e-6;
    let rule = Rule::gk61();
    let error_bound = Tolerance::Absolute(1e-14);
    let alpha: i32 = 1;

    let lower = -1.0;
    let upper = 1.0;

    let function = util::Function16::new(alpha);

    let integral =
        Adaptive::new(&function, &rule, Limits::new(lower, upper), error_bound, 3).unwrap();

    let integral_result = integral.integrate();

    println!("{}", &format!("{}", description),);
    if let Err(err) = integral_result {
        let estimate = err.estimate();
        let result = estimate.result();
        let error = estimate.error();
        let iterations = estimate.iterations();
        let evaluations = estimate.evaluations();

        if let IntegrationErrorKind::MaximumIterationsReached(iter) = err.kind() {
            util::test_relative_error(
                result,
                exp_result,
                abs_error_bound,
                &format!("{} result", description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_error_bound,
                &format!("{} error", description),
            )?;
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations", description),
            )?;
            util::test_int(
                iterations,
                iter,
                &format!("{} iterations (2) reverse", description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations", description),
            )?;
        } else {
            return Err(format!(
                "Wrong error type. Expected {} but recieved {:?}",
                stringify!(MaximumIterationsReached(iter)),
                err.kind()
            ));
        }
    } else {
        return Err(String::from(
            "Integration returned Ok, but should have been an error.",
        ));
    }

    let rule = Rule::gk61();
    let error_bound = Tolerance::Absolute(1e-14);
    let lower = 1.0;
    let upper = -1.0;

    let integral =
        Adaptive::new(&function, &rule, Limits::new(lower, upper), error_bound, 3).unwrap();

    let integral_result = integral.integrate();

    if let Err(err) = integral_result {
        let estimate = err.estimate();
        let result = estimate.result();
        let error = estimate.error();
        let iterations = estimate.iterations();
        let evaluations = estimate.evaluations();

        if let IntegrationErrorKind::MaximumIterationsReached(iter) = err.kind() {
            util::test_relative_error(
                result,
                -exp_result,
                abs_error_bound,
                &format!("{} result reverse", description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_error_bound,
                &format!("{} error reverse", description),
            )?;
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations reverse", description),
            )?;
            util::test_int(
                iterations,
                iter,
                &format!("{} iterations (2) reverse", description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations reverse", description),
            )?;
        } else {
            return Err(format!(
                "Wrong error type. Expected {} but recieved {:?}",
                stringify!(MaximumIterationsReached(iter)),
                err.kind()
            ));
        }
    } else {
        return Err(String::from(
            "Integration returned Ok, but should have been an error.",
        ));
    }

    Ok(())
}

// Test Function16 with 51 point adaptive integration and absolute error
// bound. Should terminate due to singularity detected (singularity at x=-0.1).
#[test]
fn test_adaptive_singularity_51() -> Result<(), String> {
    let exp_iterations = 51;
    let exp_evaluations = 5151;

    let error_bound = Tolerance::Absolute(1e-14);

    let rule = Rule::gk51();
    let alpha = 2;

    let lower = -1.0;
    let upper = 1.0;

    let function = util::Function16 { alpha };

    let integral = Adaptive::new(
        &function,
        &rule,
        Limits::new(lower, upper),
        error_bound,
        1000,
    )
    .unwrap();

    let integral_result = integral.integrate();

    if let Err(err) = integral_result {
        if let IntegrationErrorKind::BadIntegrandBehaviour(limits) = err.kind() {
            let estimate = err.estimate();
            let iterations = estimate.iterations();
            let evaluations = estimate.evaluations();

            assert!(limits.lower() < -0.1f64);
            assert!(limits.upper() > -0.1f64);

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

    let rule = Rule::gk51();
    let error_bound = Tolerance::Absolute(1e-14);
    let lower = 1.0;
    let upper = -1.0;

    let function = util::Function16 { alpha };

    let integral = Adaptive::new(
        &function,
        &rule,
        Limits::new(lower, upper),
        error_bound,
        1000,
    )
    .unwrap();

    let integral_result = integral.integrate();

    if let Err(err) = integral_result {
        if let IntegrationErrorKind::BadIntegrandBehaviour(limits) = err.kind() {
            let estimate = err.estimate();
            let iterations = estimate.iterations();
            let evaluations = estimate.evaluations();

            assert!(limits.lower() > -0.1f64);
            assert!(limits.upper() < -0.1f64);

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
