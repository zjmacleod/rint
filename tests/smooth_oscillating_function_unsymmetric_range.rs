// Test the basic Gauss-Kronrod integration rules with a smooth oscillating
// function over an unsymmetric range. This should find any discrepancies in
// the abscissae.
//
// Ported from gsl-2.6/integration/test.c

use rint::integration::GaussKronrodBasic;
use rint::rule::{
    GaussKronrod15, GaussKronrod21, GaussKronrod31, GaussKronrod41, GaussKronrod51, GaussKronrod61,
};

mod util;

#[test]
fn gauss_kronrod_15() -> Result<(), String> {
    let exp_result = -7.238969575483799046E-01;
    let exp_abserr = 8.760080200939757174E-06;
    let exp_resabs = 1.165564172429140788E+00;
    let exp_resasc = 9.334560307787327371E-01;

    let abs_error_bound = 1.0e-15;
    let rel_error_bound = 1.0e-7;
    let alpha = 1.3;

    let lower = 0.3;
    let upper = 2.71;

    let function = util::Function3 { alpha };
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
        "basic(f3,15) oscillating result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f3,15) oscillating abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f3,15) oscillating resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f3,15) oscillating resasc",
    )?;

    let lower = 2.71;
    let upper = 0.3;

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
        "basic(f3,15) reverse result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f3,15) reverse abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f3,15) reverse resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f3,15) reverse resasc",
    )?;

    Ok(())
}

#[test]
fn gauss_kronrod_21() -> Result<(), String> {
    let exp_result = -7.238969575482959717E-01;
    let exp_abserr = 7.999213141433641888E-11;
    let exp_resabs = 1.150829032708484023E+00;
    let exp_resasc = 9.297591249133687619E-01;

    let abs_error_bound = 1.0e-15;
    let rel_error_bound = 1.0e-7;
    let alpha = 1.3;

    let lower = 0.3;
    let upper = 2.71;

    let function = util::Function3 { alpha };
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
        "basic(f3,21) oscillating result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f3,21) oscillating abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f3,21) oscillating resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f3,21) oscillating resasc",
    )?;

    let lower = 2.71;
    let upper = 0.3;

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
        "basic(f3,21) reverse result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f3,21) reverse abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f3,21) reverse resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f3,21) reverse resasc",
    )?;

    Ok(())
}

#[test]
fn gauss_kronrod_31() -> Result<(), String> {
    let exp_result = -7.238969575482959717E-01;
    let exp_abserr = 1.285805464427459261E-14;
    let exp_resabs = 1.158150602093290571E+00;
    let exp_resasc = 9.277828092501518853E-01;

    let abs_error_bound = 1.0e-15;
    let rel_error_bound = 1.0e-7;
    let alpha = 1.3;

    let lower = 0.3;
    let upper = 2.71;

    let function = util::Function3 { alpha };
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
        "basic(f3,31) oscillating result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f3,31) oscillating abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f3,31) oscillating resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f3,31) oscillating resasc",
    )?;

    let lower = 2.71;
    let upper = 0.3;

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
        "basic(f3,31) reverse result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f3,31) reverse abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f3,31) reverse resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f3,31) reverse resasc",
    )?;

    Ok(())
}

#[test]
fn gauss_kronrod_41() -> Result<(), String> {
    let exp_result = -7.238969575482959717E-01;
    let exp_abserr = 1.286535726271015626E-14;
    let exp_resabs = 1.158808363486595328E+00;
    let exp_resasc = 9.264382258645686985E-01;

    let abs_error_bound = 1.0e-15;
    let rel_error_bound = 1.0e-7;
    let alpha = 1.3;

    let lower = 0.3;
    let upper = 2.71;

    let function = util::Function3 { alpha };
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
        "basic(f3,41) oscillating result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f3,41) oscillating abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f3,41) oscillating resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f3,41) oscillating resasc",
    )?;

    let lower = 2.71;
    let upper = 0.3;

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
        "basic(f3,41) reverse result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f3,41) reverse abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f3,41) reverse resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f3,41) reverse resasc",
    )?;

    Ok(())
}

#[test]
fn gauss_kronrod_51() -> Result<(), String> {
    let exp_result = -7.238969575482961938E-01;
    let exp_abserr = 1.285290995039385778E-14;
    let exp_resabs = 1.157687209264406381E+00;
    let exp_resasc = 9.264666884071264263E-01;

    let abs_error_bound = 1.0e-15;
    let rel_error_bound = 1.0e-7;
    let alpha = 1.3;

    let lower = 0.3;
    let upper = 2.71;

    let function = util::Function3 { alpha };
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
        "basic(f3,51) oscillating result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f3,51) oscillating abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f3,51) oscillating resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f3,51) oscillating resasc",
    )?;

    let lower = 2.71;
    let upper = 0.3;

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
        "basic(f3,51) reverse result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f3,51) reverse abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f3,51) reverse resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f3,51) reverse resasc",
    )?;

    Ok(())
}

#[test]
fn gauss_kronrod_61() -> Result<(), String> {
    let exp_result = -7.238969575482959717E-01;
    let exp_abserr = 1.286438572027470736E-14;
    let exp_resabs = 1.158720854723590099E+00;
    let exp_resasc = 9.270469641771273972E-01;

    let abs_error_bound = 1.0e-15;
    let rel_error_bound = 1.0e-7;
    let alpha = 1.3;

    let lower = 0.3;
    let upper = 2.71;

    let function = util::Function3 { alpha };
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
        "basic(f3,61) oscillating result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f3,61) oscillating abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f3,61) oscillating resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f3,61) oscillating resasc",
    )?;

    let lower = 2.71;
    let upper = 0.3;

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
        "basic(f3,61) reverse result",
    )?;
    util::test_relative_error(
        abserr,
        exp_abserr,
        rel_error_bound,
        "basic(f3,61) reverse abserr",
    )?;
    util::test_relative_error(
        resabs,
        exp_resabs,
        abs_error_bound,
        "basic(f3,61) reverse resabs",
    )?;
    util::test_relative_error(
        resasc,
        exp_resasc,
        abs_error_bound,
        "basic(f3,61) reverse resasc",
    )?;

    Ok(())
}
