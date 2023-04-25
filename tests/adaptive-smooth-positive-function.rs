use gauss_kronrod_integration::integration::{
    Adaptive, ErrorBound, GaussKronrodAdaptive, IntegrationError,
};
use gauss_kronrod_integration::rule::GaussKronrod15;

mod util;

#[test]
fn test_adaptive_smooth() -> Result<(), String> {
    let exp_result = 7.716049382715854665E-02;
    let exp_abserr = 6.679384885865053037E-12;
    let exp_last = 6;

    let alpha = 2.6;

    let lower = 0.0;
    let upper = 1.0;

    let function = util::Function1 { alpha };
    let rule = GaussKronrod15;

    let error_bound = ErrorBound::Relative(1e-10);
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

    let rel_error_bound = 1e-6;
    let abs_error_bound = 1e-15;

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

    println!("RESUL:\t{result}");
    println!("ERROR:\t{error}");
    println!("ITERA:\t{iterations}");

    Ok(())
}
