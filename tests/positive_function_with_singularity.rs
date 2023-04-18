// Test the basic Gauss-Kronrod integration rules with a positive function that
// has a singularity. This should give large values of the absolute error which
// would find discrepancies in the absolute error calculation.
//
// Ported from gsl-2.6/integration/test.c

use gauss_kronrod_integration::integration::GaussKronrodBasic;
use gauss_kronrod_integration::rule::{
    GaussKronrod15, GaussKronrod21, GaussKronrod31, GaussKronrod41,
    GaussKronrod51, GaussKronrod61,
};

mod util;

#[test]
fn gauss_kronrod_15() -> Result<(), String> {
    let exp_result = 1.555688196612745777E+01;
    let exp_abserr = 2.350164577239293706E+01;
    let exp_resabs = 1.555688196612745777E+01;
    let exp_resasc = 2.350164577239293706E+01;

    let abs_error_bound = 1.0e-15;
    let rel_error_bound = 1.0e-7;
    let alpha = -0.9;

    let lower = 0.0;
    let upper = 1.0;

    let function = util::Function1 { alpha };
    let rule = GaussKronrod15;

    let integral = GaussKronrodBasic::new(lower, upper, rule, function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        exp_result,
        abs_error_bound,
        "qk15(f1) singular result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "qk15(f1) singular abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "qk15(f1) singular resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "qk15(f1) singular resasc",
    )?;

    let lower = 1.0;
    let upper = 0.0;

    let function = util::Function1 { alpha };
    let rule = GaussKronrod15;

    let integral = GaussKronrodBasic::new(lower, upper, rule, function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        -exp_result,
        abs_error_bound,
        "qk15(f1) reverse result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "qk15(f1) reverse abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "qk15(f1) reverse resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "qk15(f1) reverse resasc",
    )?;

    Ok(())
}

#[test]
fn gauss_kronrod_21() -> Result<(), String> {
    let exp_result = 1.799045317938126232E+01;
    let exp_abserr = 2.782360287710622515E+01;
    let exp_resabs = 1.799045317938126232E+01;
    let exp_resasc = 2.782360287710622515E+01;

    let abs_error_bound = 1.0e-15;
    let rel_error_bound = 1.0e-7;
    let alpha = -0.9;

    let lower = 0.0;
    let upper = 1.0;

    let function = util::Function1 { alpha };
    let rule = GaussKronrod21;

    let integral = GaussKronrodBasic::new(lower, upper, rule, function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        exp_result,
        abs_error_bound,
        "qk21(f1) singular result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "qk21(f1) singular abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "qk21(f1) singular resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "qk21(f1) singular resasc",
    )?;

    let lower = 1.0;
    let upper = 0.0;

    let function = util::Function1 { alpha };
    let rule = GaussKronrod21;

    let integral = GaussKronrodBasic::new(lower, upper, rule, function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        -exp_result,
        abs_error_bound,
        "qk21(f1) reverse result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "qk21(f1) reverse abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "qk21(f1) reverse resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "qk21(f1) reverse resasc",
    )?;

    Ok(())
}

#[test]
fn gauss_kronrod_31() -> Result<(), String> {
    let exp_result = 2.081873305159121657E+01;
    let exp_abserr = 3.296500137482590276E+01;
    let exp_resabs = 2.081873305159121301E+01;
    let exp_resasc = 3.296500137482590276E+01;

    let abs_error_bound = 1.0e-15;
    let rel_error_bound = 1.0e-7;
    let alpha = -0.9;

    let lower = 0.0;
    let upper = 1.0;

    let function = util::Function1 { alpha };
    let rule = GaussKronrod31;

    let integral = GaussKronrodBasic::new(lower, upper, rule, function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        exp_result,
        abs_error_bound,
        "qk31(f1) singular result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "qk31(f1) singular abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "qk31(f1) singular resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "qk31(f1) singular resasc",
    )?;

    let lower = 1.0;
    let upper = 0.0;

    let function = util::Function1 { alpha };
    let rule = GaussKronrod31;

    let integral = GaussKronrodBasic::new(lower, upper, rule, function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        -exp_result,
        abs_error_bound,
        "qk31(f1) reverse result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "qk31(f1) reverse abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "qk31(f1) reverse resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "qk31(f1) reverse resasc",
    )?;

    Ok(())
}

#[test]
fn gauss_kronrod_41() -> Result<(), String> {
    let exp_result = 2.288677623903126701E+01;
    let exp_abserr = 3.671538820274916048E+01;
    let exp_resabs = 2.288677623903126701E+01;
    let exp_resasc = 3.671538820274916048E+01;

    let abs_error_bound = 1.0e-15;
    let rel_error_bound = 1.0e-7;
    let alpha = -0.9;

    let lower = 0.0;
    let upper = 1.0;

    let function = util::Function1 { alpha };
    let rule = GaussKronrod41;

    let integral = GaussKronrodBasic::new(lower, upper, rule, function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        exp_result,
        abs_error_bound,
        "qk41(f1) singular result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "qk41(f1) singular abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "qk41(f1) singular resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "qk41(f1) singular resasc",
    )?;

    let lower = 1.0;
    let upper = 0.0;

    let function = util::Function1 { alpha };
    let rule = GaussKronrod41;

    let integral = GaussKronrodBasic::new(lower, upper, rule, function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        -exp_result,
        abs_error_bound,
        "qk41(f1) reverse result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "qk41(f1) reverse abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "qk41(f1) reverse resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "qk41(f1) reverse resasc",
    )?;

    Ok(())
}

#[test]
fn gauss_kronrod_51() -> Result<(), String> {
    let exp_result = 2.449953612016972215E+01;
    let exp_abserr = 3.967771249391228849E+01;
    let exp_resabs = 2.449953612016972215E+01;
    let exp_resasc = 3.967771249391228849E+01;

    let abs_error_bound = 1.0e-15;
    let rel_error_bound = 1.0e-7;
    let alpha = -0.9;

    let lower = 0.0;
    let upper = 1.0;

    let function = util::Function1 { alpha };
    let rule = GaussKronrod51;

    let integral = GaussKronrodBasic::new(lower, upper, rule, function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        exp_result,
        abs_error_bound,
        "qk51(f1) singular result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "qk51(f1) singular abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "qk51(f1) singular resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "qk51(f1) singular resasc",
    )?;

    let lower = 1.0;
    let upper = 0.0;

    let function = util::Function1 { alpha };
    let rule = GaussKronrod51;

    let integral = GaussKronrodBasic::new(lower, upper, rule, function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        -exp_result,
        abs_error_bound,
        "qk51(f1) reverse result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "qk51(f1) reverse abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "qk51(f1) reverse resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "qk51(f1) reverse resasc",
    )?;

    Ok(())
}

#[test]
fn gauss_kronrod_61() -> Result<(), String> {
    let exp_result = 2.583030240976628988E+01;
    let exp_abserr = 4.213750493076978643E+01;
    let exp_resabs = 2.583030240976628988E+01;
    let exp_resasc = 4.213750493076978643E+01;

    let abs_error_bound = 1.0e-15;
    let rel_error_bound = 1.0e-7;
    let alpha = -0.9;

    let lower = 0.0;
    let upper = 1.0;

    let function = util::Function1 { alpha };
    let rule = GaussKronrod61;

    let integral = GaussKronrodBasic::new(lower, upper, rule, function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        exp_result,
        abs_error_bound,
        "qk61(f1) singular result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "qk61(f1) singular abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "qk61(f1) singular resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "qk61(f1) singular resasc",
    )?;

    let lower = 1.0;
    let upper = 0.0;

    let function = util::Function1 { alpha };
    let rule = GaussKronrod61;

    let integral = GaussKronrodBasic::new(lower, upper, rule, function);

    let integral_result = integral.integrate();
    let result = integral_result.result();
    let abserr = integral_result.error();
    let resabs = integral_result.result_abs();
    let resasc = integral_result.result_asc();

    util::test_relative_error(
        result,
        -exp_result,
        abs_error_bound,
        "qk61(f1) reverse result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "qk61(f1) reverse abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "qk61(f1) reverse resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "qk61(f1) reverse resasc",
    )?;

    Ok(())
}
