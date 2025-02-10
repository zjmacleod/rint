use rint::quadrature::Basic;
use rint::quadrature::Rule;

mod util;

#[test]
fn catalan_test_basic_15() -> Result<(), String> {
    const EXPECTED: f64 = util::CATALAN;

    let rule = Rule::gk15();
    let catalan = util::Catalan1::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk15();
    let catalan = util::Catalan2::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk15();
    let catalan = util::Catalan3::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk15();
    let catalan = util::Catalan4::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk15();
    let catalan = util::ComplexCatalan12::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk15();
    let catalan = util::ComplexCatalan13::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk15();
    let catalan = util::ComplexCatalan14::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk15();
    let catalan = util::ComplexCatalan23::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk15();
    let catalan = util::ComplexCatalan24::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk15();
    let catalan = util::ComplexCatalan34::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    Ok(())
}

#[test]
fn catalan_test_basic_21() -> Result<(), String> {
    const EXPECTED: f64 = util::CATALAN;

    let rule = Rule::gk21();
    let catalan = util::Catalan1::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk21();
    let catalan = util::Catalan2::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk21();
    let catalan = util::Catalan3::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk21();
    let catalan = util::Catalan4::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk21();
    let catalan = util::ComplexCatalan12::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk21();
    let catalan = util::ComplexCatalan13::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk21();
    let catalan = util::ComplexCatalan14::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk21();
    let catalan = util::ComplexCatalan23::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk21();
    let catalan = util::ComplexCatalan24::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk21();
    let catalan = util::ComplexCatalan34::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    Ok(())
}

#[test]
fn catalan_test_basic_31() -> Result<(), String> {
    const EXPECTED: f64 = util::CATALAN;

    let rule = Rule::gk31();
    let catalan = util::Catalan1::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk31();
    let catalan = util::Catalan2::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk31();
    let catalan = util::Catalan3::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk31();
    let catalan = util::Catalan4::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk31();
    let catalan = util::ComplexCatalan12::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk31();
    let catalan = util::ComplexCatalan13::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk31();
    let catalan = util::ComplexCatalan14::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk31();
    let catalan = util::ComplexCatalan23::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk31();
    let catalan = util::ComplexCatalan24::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk31();
    let catalan = util::ComplexCatalan34::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    Ok(())
}

#[test]
fn catalan_test_basic_41() -> Result<(), String> {
    const EXPECTED: f64 = util::CATALAN;

    let rule = Rule::gk41();
    let catalan = util::Catalan1::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk41();
    let catalan = util::Catalan2::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk41();
    let catalan = util::Catalan3::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk41();
    let catalan = util::Catalan4::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk41();
    let catalan = util::ComplexCatalan12::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk41();
    let catalan = util::ComplexCatalan13::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk41();
    let catalan = util::ComplexCatalan14::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk41();
    let catalan = util::ComplexCatalan23::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk41();
    let catalan = util::ComplexCatalan24::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk41();
    let catalan = util::ComplexCatalan34::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    Ok(())
}

#[test]
fn catalan_test_basic_51() -> Result<(), String> {
    const EXPECTED: f64 = util::CATALAN;

    let rule = Rule::gk51();
    let catalan = util::Catalan1::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk51();
    let catalan = util::Catalan2::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk51();
    let catalan = util::Catalan3::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk51();
    let catalan = util::Catalan4::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk51();
    let catalan = util::ComplexCatalan12::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk51();
    let catalan = util::ComplexCatalan13::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk51();
    let catalan = util::ComplexCatalan14::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk51();
    let catalan = util::ComplexCatalan23::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk51();
    let catalan = util::ComplexCatalan24::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk51();
    let catalan = util::ComplexCatalan34::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    Ok(())
}

#[test]
fn catalan_test_basic_61() -> Result<(), String> {
    const EXPECTED: f64 = util::CATALAN;

    let rule = Rule::gk61();
    let catalan = util::Catalan1::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk61();
    let catalan = util::Catalan2::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk61();
    let catalan = util::Catalan3::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk61();
    let catalan = util::Catalan4::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();

    assert!(abs_actual_error < error);
    //println!("res: {result:?}\nerr: {error}\nacre: {abs_actual_error}\n");

    let rule = Rule::gk61();
    let catalan = util::ComplexCatalan12::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk61();
    let catalan = util::ComplexCatalan13::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk61();
    let catalan = util::ComplexCatalan14::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk61();
    let catalan = util::ComplexCatalan23::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk61();
    let catalan = util::ComplexCatalan24::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    let rule = Rule::gk61();
    let catalan = util::ComplexCatalan34::new();
    let limits = catalan.limits();
    let integral = Basic::new(&catalan, &rule, limits).integrate();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();

    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    //println!(
    //    "res: {result:?}\nerr: {error}\nacre: {abs_actual_error_re}\nacin: {abs_actual_error_im}\n"
    //);

    Ok(())
}
