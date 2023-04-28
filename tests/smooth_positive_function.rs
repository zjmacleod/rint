// Test the basic Gauss-Kronrod integration rules with a smooth positive function.
// Ported from gsl-2.6/integration/test.c

use gauss_kronrod_integration::integration::GaussKronrodBasic;
use gauss_kronrod_integration::rule::{
    GaussKronrod15, GaussKronrod21, GaussKronrod31, GaussKronrod41, GaussKronrod51, GaussKronrod61,
};

mod util;

#[test]
fn gauss_kronrod_15() -> Result<(), String> {
    let exp_result = 7.716049357767090777E-02;
    let exp_abserr = 2.990224871000550874E-06;
    let exp_resabs = 7.716049357767090777E-02;
    let exp_resasc = 4.434273814139995384E-02;

    let abs_error_bound = 1.0e-15;
    let rel_error_bound = 1.0e-7;
    let alpha = 2.6;

    let lower = 0.0;
    let upper = 1.0;

    let function = util::Function1 { alpha };
    let rule = GaussKronrod15;

    let integral = GaussKronrodBasic::new(lower, upper, rule, &function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        exp_result,
        abs_error_bound,
        "basic(f1,15) smooth result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f1,15) smooth abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f1,15) smooth resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f1,15) smooth resasc",
    )?;

    let lower = 1.0;
    let upper = 0.0;

    let integral = GaussKronrodBasic::new(lower, upper, rule, &function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        -exp_result,
        abs_error_bound,
        "basic(f1,15) reverse result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f1,15) reverse abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f1,15) reverse resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f1,15) reverse resasc",
    )?;

    Ok(())
}

#[test]
fn gauss_kronrod_21() -> Result<(), String> {
    let exp_result = 7.716049379303084599E-02;
    let exp_abserr = 9.424302194248481445E-08;
    let exp_resabs = 7.716049379303084599E-02;
    let exp_resasc = 4.434311425038358484E-02;

    let abs_error_bound = 1.0e-15;
    let rel_error_bound = 1.0e-7;
    let alpha = 2.6;

    let lower = 0.0;
    let upper = 1.0;

    let function = util::Function1 { alpha };
    let rule = GaussKronrod21;

    let integral = GaussKronrodBasic::new(lower, upper, rule, &function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        exp_result,
        abs_error_bound,
        "basic(f1,21) smooth result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f1,21) smooth abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f1,21) smooth resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f1,21) smooth resasc",
    )?;

    let lower = 1.0;
    let upper = 0.0;

    let integral = GaussKronrodBasic::new(lower, upper, rule, &function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        -exp_result,
        abs_error_bound,
        "basic(f1,21) reverse result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f1,21) reverse abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f1,21) reverse resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f1,21) reverse resasc",
    )?;

    Ok(())
}

#[test]
fn gauss_kronrod_31() -> Result<(), String> {
    let exp_result = 7.716049382494900855E-02;
    let exp_abserr = 1.713503193600029893E-09;
    let exp_resabs = 7.716049382494900855E-02;
    let exp_resasc = 4.427995051868838933E-02;

    let abs_error_bound = 1.0e-15;
    let rel_error_bound = 1.0e-7;
    let alpha = 2.6;

    let lower = 0.0;
    let upper = 1.0;

    let function = util::Function1 { alpha };
    let rule = GaussKronrod31;

    let integral = GaussKronrodBasic::new(lower, upper, rule, &function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        exp_result,
        abs_error_bound,
        "basic(f1,31) smooth result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f1,31) smooth abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f1,31) smooth resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f1,31) smooth resasc",
    )?;

    let lower = 1.0;
    let upper = 0.0;

    let integral = GaussKronrodBasic::new(lower, upper, rule, &function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        -exp_result,
        abs_error_bound,
        "basic(f1,31) reverse result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f1,31) reverse abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f1,31) reverse resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f1,31) reverse resasc",
    )?;

    Ok(())
}

#[test]
fn gauss_kronrod_41() -> Result<(), String> {
    let exp_result = 7.716049382681375302E-02;
    let exp_abserr = 9.576386660975511224E-11;
    let exp_resabs = 7.716049382681375302E-02;
    let exp_resasc = 4.421521169637691873E-02;

    let abs_error_bound = 1.0e-15;
    let rel_error_bound = 1.0e-7;
    let alpha = 2.6;

    let lower = 0.0;
    let upper = 1.0;

    let function = util::Function1 { alpha };
    let rule = GaussKronrod41;

    let integral = GaussKronrodBasic::new(lower, upper, rule, &function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        exp_result,
        abs_error_bound,
        "basic(f1,41) smooth result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f1,41) smooth abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f1,41) smooth resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f1,41) smooth resasc",
    )?;

    let lower = 1.0;
    let upper = 0.0;

    let integral = GaussKronrodBasic::new(lower, upper, rule, &function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        -exp_result,
        abs_error_bound,
        "basic(f1,41) reverse result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f1,41) reverse abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f1,41) reverse resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f1,41) reverse resasc",
    )?;

    Ok(())
}

#[test]
fn gauss_kronrod_51() -> Result<(), String> {
    let exp_result = 7.716049382708510540E-02;
    let exp_abserr = 1.002079980317363772E-11;
    let exp_resabs = 7.716049382708510540E-02;
    let exp_resasc = 4.416474291216854892E-02;

    let abs_error_bound = 1.0e-15;
    let rel_error_bound = 1.0e-5;
    let alpha = 2.6;

    let lower = 0.0;
    let upper = 1.0;

    let function = util::Function1 { alpha };
    let rule = GaussKronrod51;

    let integral = GaussKronrodBasic::new(lower, upper, rule, &function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        exp_result,
        abs_error_bound,
        "basic(f1,51) smooth result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f1,51) smooth abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f1,51) smooth resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f1,51) smooth resasc",
    )?;

    let lower = 1.0;
    let upper = 0.0;

    let integral = GaussKronrodBasic::new(lower, upper, rule, &function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        -exp_result,
        abs_error_bound,
        "basic(f1,51) reverse result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f1,51) reverse abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f1,51) reverse resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f1,51) reverse resasc",
    )?;

    Ok(())
}

#[test]
fn gauss_kronrod_61() -> Result<(), String> {
    let exp_result = 7.716049382713800753E-02;
    let exp_abserr = 1.566060362296155616E-12;
    let exp_resabs = 7.716049382713800753E-02;
    let exp_resasc = 4.419287685934316506E-02;

    let abs_error_bound = 1.0e-15;
    let rel_error_bound = 1.0e-5;
    let alpha = 2.6;

    let lower = 0.0;
    let upper = 1.0;

    let function = util::Function1 { alpha };
    let rule = GaussKronrod61;

    let integral = GaussKronrodBasic::new(lower, upper, rule, &function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        exp_result,
        abs_error_bound,
        "basic(f1,61) smooth result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f1,61) smooth abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f1,61) smooth resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f1,61) smooth resasc",
    )?;

    let lower = 1.0;
    let upper = 0.0;

    let integral = GaussKronrodBasic::new(lower, upper, rule, &function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        -exp_result,
        abs_error_bound,
        "basic(f1,61) reverse result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f1,61) reverse abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f1,61) reverse resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f1,61) reverse resasc",
    )?;

    Ok(())
}
