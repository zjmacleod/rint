use num_complex::ComplexFloat;
use rint::quadrature::Adaptive;
use rint::quadrature::Rule;
use rint::Tolerance;

mod util;

#[test]
fn catalan_test_adaptive_relative_tol_15() -> Result<(), String> {
    const EXPECTED: f64 = util::CATALAN;
    const TOL: f64 = 1.0e-13;

    let rule = Rule::gk15();
    let catalan = util::Catalan1::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 45);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk15();
    let catalan = util::Catalan2::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 2);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk15();
    let catalan = util::Catalan3::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 44);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk15();
    let catalan = util::Catalan4::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 27);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk15();
    let catalan = util::ComplexCatalan12::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 44);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk15();
    let catalan = util::ComplexCatalan13::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 45);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk15();
    let catalan = util::ComplexCatalan14::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 70);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk15();
    let catalan = util::ComplexCatalan23::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 44);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk15();
    let catalan = util::ComplexCatalan24::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 27);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk15();
    let catalan = util::ComplexCatalan34::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 69);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    Ok(())
}

#[test]
fn catalan_test_adaptive_absolute_tol_15() -> Result<(), String> {
    const EXPECTED: f64 = util::CATALAN;
    const TOL: f64 = 1.0e-7;

    let rule = Rule::gk15();
    let catalan = util::Catalan1::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 24);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk15();
    let catalan = util::Catalan2::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 1);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk15();
    let catalan = util::Catalan3::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 24);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk15();
    let catalan = util::Catalan4::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 13);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk15();
    let catalan = util::ComplexCatalan12::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 24);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk15();
    let catalan = util::ComplexCatalan13::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 25);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk15();
    let catalan = util::ComplexCatalan14::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 37);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk15();
    let catalan = util::ComplexCatalan23::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 24);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk15();
    let catalan = util::ComplexCatalan24::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 13);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk15();
    let catalan = util::ComplexCatalan34::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 36);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    Ok(())
}

#[test]
fn catalan_test_adaptive_relative_tol_21() -> Result<(), String> {
    const EXPECTED: f64 = util::CATALAN;
    const TOL: f64 = 1.0e-13;

    let rule = Rule::gk21();
    let catalan = util::Catalan1::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 45);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk21();
    let catalan = util::Catalan2::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 1);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk21();
    let catalan = util::Catalan3::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 44);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk21();
    let catalan = util::Catalan4::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 25);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk21();
    let catalan = util::ComplexCatalan12::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 44);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk21();
    let catalan = util::ComplexCatalan13::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 44);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk21();
    let catalan = util::ComplexCatalan14::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 68);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk21();
    let catalan = util::ComplexCatalan23::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 44);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk21();
    let catalan = util::ComplexCatalan24::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 25);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk21();
    let catalan = util::ComplexCatalan34::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 68);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    Ok(())
}

#[test]
fn catalan_test_adaptive_absolute_tol_21() -> Result<(), String> {
    const EXPECTED: f64 = util::CATALAN;
    const TOL: f64 = 1.0e-7;

    let rule = Rule::gk21();
    let catalan = util::Catalan1::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 24);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk21();
    let catalan = util::Catalan2::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 1);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk21();
    let catalan = util::Catalan3::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 24);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk21();
    let catalan = util::Catalan4::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 12);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk21();
    let catalan = util::ComplexCatalan12::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 24);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk21();
    let catalan = util::ComplexCatalan13::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 25);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk21();
    let catalan = util::ComplexCatalan14::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 35);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk21();
    let catalan = util::ComplexCatalan23::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 24);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk21();
    let catalan = util::ComplexCatalan24::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 12);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk21();
    let catalan = util::ComplexCatalan34::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 35);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    Ok(())
}

#[test]
fn catalan_test_adaptive_relative_tol_31() -> Result<(), String> {
    const EXPECTED: f64 = util::CATALAN;
    const TOL: f64 = 1.0e-13;

    let rule = Rule::gk31();
    let catalan = util::Catalan1::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 43);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk31();
    let catalan = util::Catalan2::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 1);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk31();
    let catalan = util::Catalan3::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 43);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk31();
    let catalan = util::Catalan4::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 24);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk31();
    let catalan = util::ComplexCatalan12::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 43);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk31();
    let catalan = util::ComplexCatalan13::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 43);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk31();
    let catalan = util::ComplexCatalan14::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 65);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk31();
    let catalan = util::ComplexCatalan23::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 43);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk31();
    let catalan = util::ComplexCatalan24::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 23);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk31();
    let catalan = util::ComplexCatalan34::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 65);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    Ok(())
}

#[test]
fn catalan_test_adaptive_absolute_tol_31() -> Result<(), String> {
    const EXPECTED: f64 = util::CATALAN;
    const TOL: f64 = 1.0e-7;

    let rule = Rule::gk31();
    let catalan = util::Catalan1::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 23);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk31();
    let catalan = util::Catalan2::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 1);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk31();
    let catalan = util::Catalan3::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 23);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk31();
    let catalan = util::Catalan4::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 10);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk31();
    let catalan = util::ComplexCatalan12::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 23);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk31();
    let catalan = util::ComplexCatalan13::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 24);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk31();
    let catalan = util::ComplexCatalan14::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 33);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk31();
    let catalan = util::ComplexCatalan23::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 23);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk31();
    let catalan = util::ComplexCatalan24::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 10);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk31();
    let catalan = util::ComplexCatalan34::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 32);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    Ok(())
}

#[test]
fn catalan_test_adaptive_relative_tol_41() -> Result<(), String> {
    const EXPECTED: f64 = util::CATALAN;
    const TOL: f64 = 1.0e-13;

    let rule = Rule::gk41();
    let catalan = util::Catalan1::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 42);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk41();
    let catalan = util::Catalan2::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 1);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk41();
    let catalan = util::Catalan3::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 42);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk41();
    let catalan = util::Catalan4::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 22);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk41();
    let catalan = util::ComplexCatalan12::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 42);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk41();
    let catalan = util::ComplexCatalan13::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 42);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk41();
    let catalan = util::ComplexCatalan14::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 63);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk41();
    let catalan = util::ComplexCatalan23::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 41);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk41();
    let catalan = util::ComplexCatalan24::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 22);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk41();
    let catalan = util::ComplexCatalan34::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 63);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    Ok(())
}

#[test]
fn catalan_test_adaptive_absolute_tol_41() -> Result<(), String> {
    const EXPECTED: f64 = util::CATALAN;
    const TOL: f64 = 1.0e-7;

    let rule = Rule::gk41();
    let catalan = util::Catalan1::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 22);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk41();
    let catalan = util::Catalan2::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 1);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk41();
    let catalan = util::Catalan3::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 22);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk41();
    let catalan = util::Catalan4::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 9);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk41();
    let catalan = util::ComplexCatalan12::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 22);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk41();
    let catalan = util::ComplexCatalan13::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 22);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk41();
    let catalan = util::ComplexCatalan14::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 30);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk41();
    let catalan = util::ComplexCatalan23::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 22);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk41();
    let catalan = util::ComplexCatalan24::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 9);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk41();
    let catalan = util::ComplexCatalan34::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 30);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    Ok(())
}

#[test]
fn catalan_test_adaptive_relative_tol_51() -> Result<(), String> {
    const EXPECTED: f64 = util::CATALAN;
    const TOL: f64 = 1.0e-13;

    let rule = Rule::gk51();
    let catalan = util::Catalan1::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 41);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk51();
    let catalan = util::Catalan2::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 1);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk51();
    let catalan = util::Catalan3::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 41);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk51();
    let catalan = util::Catalan4::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 22);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk51();
    let catalan = util::ComplexCatalan12::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 41);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk51();
    let catalan = util::ComplexCatalan13::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 41);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk51();
    let catalan = util::ComplexCatalan14::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 61);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk51();
    let catalan = util::ComplexCatalan23::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 40);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk51();
    let catalan = util::ComplexCatalan24::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 21);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk51();
    let catalan = util::ComplexCatalan34::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 61);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    Ok(())
}

#[test]
fn catalan_test_adaptive_absolute_tol_51() -> Result<(), String> {
    const EXPECTED: f64 = util::CATALAN;
    const TOL: f64 = 1.0e-7;

    let rule = Rule::gk51();
    let catalan = util::Catalan1::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 21);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk51();
    let catalan = util::Catalan2::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 1);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk51();
    let catalan = util::Catalan3::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 21);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk51();
    let catalan = util::Catalan4::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 8);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk51();
    let catalan = util::ComplexCatalan12::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 21);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk51();
    let catalan = util::ComplexCatalan13::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 21);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk51();
    let catalan = util::ComplexCatalan14::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 28);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk51();
    let catalan = util::ComplexCatalan23::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 21);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk51();
    let catalan = util::ComplexCatalan24::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 8);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk51();
    let catalan = util::ComplexCatalan34::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 28);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    Ok(())
}

#[test]
fn catalan_test_adaptive_relative_tol_61() -> Result<(), String> {
    const EXPECTED: f64 = util::CATALAN;
    const TOL: f64 = 1.0e-13;

    let rule = Rule::gk61();
    let catalan = util::Catalan1::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 41);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk61();
    let catalan = util::Catalan2::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 1);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk61();
    let catalan = util::Catalan3::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 40);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk61();
    let catalan = util::Catalan4::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 21);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk61();
    let catalan = util::ComplexCatalan12::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 40);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk61();
    let catalan = util::ComplexCatalan13::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 40);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk61();
    let catalan = util::ComplexCatalan14::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 60);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk61();
    let catalan = util::ComplexCatalan23::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 40);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk61();
    let catalan = util::ComplexCatalan24::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 20);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk61();
    let catalan = util::ComplexCatalan34::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Relative(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 59);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    Ok(())
}

#[test]
fn catalan_test_adaptive_absolute_tol_61() -> Result<(), String> {
    const EXPECTED: f64 = util::CATALAN;
    const TOL: f64 = 1.0e-7;

    let rule = Rule::gk61();
    let catalan = util::Catalan1::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 20);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk61();
    let catalan = util::Catalan2::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 1);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk61();
    let catalan = util::Catalan3::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 20);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk61();
    let catalan = util::Catalan4::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 7);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let rule = Rule::gk61();
    let catalan = util::ComplexCatalan12::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 20);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk61();
    let catalan = util::ComplexCatalan13::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 21);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk61();
    let catalan = util::ComplexCatalan14::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 27);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk61();
    let catalan = util::ComplexCatalan23::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 20);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk61();
    let catalan = util::ComplexCatalan24::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 7);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let rule = Rule::gk61();
    let catalan = util::ComplexCatalan34::new();
    let limits = catalan.limits();
    let tolerance = Tolerance::Absolute(TOL);
    let integral = Adaptive::new(&catalan, &rule, limits, tolerance, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let result_re = result.re;
    let result_im = result.im;
    let error = integral.error();
    let abs_actual_error_re = (EXPECTED - result_re).abs();
    let abs_actual_error_im = (EXPECTED - result_im).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 26);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    Ok(())
}
