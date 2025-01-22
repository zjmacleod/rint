use num_complex::ComplexFloat;
use rint::quadrature::AdaptiveSingularity;
use rint::Tolerance;

mod util;

#[test]
fn catalan_test_singularity_finite_relative_tol() -> Result<(), String> {
    const EXPECTED: f64 = util::CATALAN;
    const TOL: f64 = 1.0e-13;
    const TOL56: f64 = 1.0e-12; // Catalan5 and 6 have boundary singularities

    let catalan = util::Catalan1::new();
    let limits = catalan.limits();
    let error_bound = Tolerance::Relative(TOL);
    let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 8);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let catalan = util::Catalan2::new();
    let limits = catalan.limits();
    let error_bound = Tolerance::Relative(TOL);
    let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
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

    let catalan = util::Catalan3::new();
    let limits = catalan.limits();
    let error_bound = Tolerance::Relative(TOL);
    let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 6);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let catalan = util::Catalan4::new();
    let limits = catalan.limits();
    let error_bound = Tolerance::Relative(TOL);
    let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 8);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let catalan = util::Catalan5::new();
    let limits = catalan.limits();
    let error_bound = Tolerance::Relative(TOL56);
    let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL56 * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 10);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let catalan = util::Catalan6::new();
    let limits = catalan.limits();
    let error_bound = Tolerance::Relative(TOL56);
    let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL56 * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 17);
    assert!(abs_actual_error * 1e-1 < error);
    assert!(error < tol);

    let catalan = util::ComplexCatalan12::new();
    let limits = catalan.limits();
    let error_bound = Tolerance::Relative(TOL);
    let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
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
    assert_eq!(iters, 38);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let catalan = util::ComplexCatalan13::new();
    let limits = catalan.limits();
    let error_bound = Tolerance::Relative(TOL);
    let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
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
    assert_eq!(iters, 38);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let catalan = util::ComplexCatalan14::new();
    let limits = catalan.limits();
    let error_bound = Tolerance::Relative(TOL);
    let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
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

    let catalan = util::ComplexCatalan23::new();
    let limits = catalan.limits();
    let error_bound = Tolerance::Relative(TOL);
    let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
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
    assert_eq!(iters, 37);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let catalan = util::ComplexCatalan24::new();
    let limits = catalan.limits();
    let error_bound = Tolerance::Relative(TOL);
    let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
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

    let catalan = util::ComplexCatalan34::new();
    let limits = catalan.limits();
    let error_bound = Tolerance::Relative(TOL);
    let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
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

    Ok(())
}

#[test]
fn catalan_test_singularity_finite_absolute_tol() -> Result<(), String> {
    const EXPECTED: f64 = util::CATALAN;
    const TOL: f64 = 1.0e-7;

    let catalan = util::Catalan1::new();
    let limits = catalan.limits();
    let error_bound = Tolerance::Absolute(TOL);
    let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 6);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let catalan = util::Catalan2::new();
    let limits = catalan.limits();
    let error_bound = Tolerance::Absolute(TOL);
    let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
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

    let catalan = util::Catalan3::new();
    let limits = catalan.limits();
    let error_bound = Tolerance::Absolute(TOL);
    let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 6);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let catalan = util::Catalan4::new();
    let limits = catalan.limits();
    let error_bound = Tolerance::Absolute(TOL);
    let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL;
    let iters = integral.iterations();
    assert_eq!(iters, 6);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let catalan = util::ComplexCatalan12::new();
    let limits = catalan.limits();
    let error_bound = Tolerance::Absolute(TOL);
    let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
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
    assert_eq!(iters, 18);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let catalan = util::ComplexCatalan13::new();
    let limits = catalan.limits();
    let error_bound = Tolerance::Absolute(TOL);
    let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
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
    assert_eq!(iters, 18);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let catalan = util::ComplexCatalan14::new();
    let limits = catalan.limits();
    let error_bound = Tolerance::Absolute(TOL);
    let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
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

    let catalan = util::ComplexCatalan23::new();
    let limits = catalan.limits();
    let error_bound = Tolerance::Absolute(TOL);
    let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
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
    assert_eq!(iters, 18);
    assert!(abs_actual_error_re < error);
    assert!(abs_actual_error_im < error);
    assert!(error < tol);

    let catalan = util::ComplexCatalan24::new();
    let limits = catalan.limits();
    let error_bound = Tolerance::Absolute(TOL);
    let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
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

    let catalan = util::ComplexCatalan34::new();
    let limits = catalan.limits();
    let error_bound = Tolerance::Absolute(TOL);
    let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
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
fn catalan_test_singularity_semi_infinite_relative_tol() -> Result<(), String> {
    const EXPECTED: f64 = util::CATALAN;
    const TOL: f64 = 1.0e-13;

    let catalan = util::Catalan7::new();
    let limits = catalan.lower;
    let error_bound = Tolerance::Relative(TOL);
    let integral = AdaptiveSingularity::semi_infinite_upper(catalan, limits, error_bound, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 8);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let catalan = util::Catalan8::new();
    let limits = catalan.lower;
    let error_bound = Tolerance::Relative(TOL);
    let integral = AdaptiveSingularity::semi_infinite_upper(catalan, limits, error_bound, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 4);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let catalan = util::Catalan9::new();
    let limits = catalan.lower;
    let error_bound = Tolerance::Relative(TOL);
    let integral = AdaptiveSingularity::semi_infinite_upper(catalan, limits, error_bound, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 9);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let catalan = util::Catalan10::new();
    let limits = catalan.lower;
    let error_bound = Tolerance::Relative(TOL);
    let integral = AdaptiveSingularity::semi_infinite_upper(catalan, limits, error_bound, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 32);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let catalan = util::Catalan11::new();
    let limits = catalan.lower;
    let error_bound = Tolerance::Relative(TOL);
    let integral = AdaptiveSingularity::semi_infinite_upper(catalan, limits, error_bound, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 11);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    //let catalan = util::ComplexCatalan12::new();
    //let limits = catalan.limits();
    //let error_bound = Tolerance::Relative(TOL);
    //let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
    //    .unwrap()
    //    .integrate()
    //    .unwrap();
    //let result = integral.result();
    //let result_re = result.re;
    //let result_im = result.im;
    //let error = integral.error();
    //let abs_actual_error_re = (EXPECTED - result_re).abs();
    //let abs_actual_error_im = (EXPECTED - result_im).abs();
    //let tol = TOL * result.abs();
    //let iters = integral.iterations();
    //assert_eq!(iters, 38);
    //assert!(abs_actual_error_re < error);
    //assert!(abs_actual_error_im < error);
    //assert!(error < tol);

    Ok(())
}

#[test]
fn catalan_test_singularity_semi_infinite_absolute_tol() -> Result<(), String> {
    const EXPECTED: f64 = util::CATALAN;
    const TOL: f64 = 1.0e-10;

    let catalan = util::Catalan7::new();
    let limits = catalan.lower;
    let error_bound = Tolerance::Absolute(TOL);
    let integral = AdaptiveSingularity::semi_infinite_upper(catalan, limits, error_bound, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 7);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let catalan = util::Catalan8::new();
    let limits = catalan.lower;
    let error_bound = Tolerance::Absolute(TOL);
    let integral = AdaptiveSingularity::semi_infinite_upper(catalan, limits, error_bound, 1000)
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

    let catalan = util::Catalan9::new();
    let limits = catalan.lower;
    let error_bound = Tolerance::Absolute(TOL);
    let integral = AdaptiveSingularity::semi_infinite_upper(catalan, limits, error_bound, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 5);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    let catalan = util::Catalan10::new();
    let limits = catalan.lower;
    let error_bound = Tolerance::Absolute(TOL);
    let integral = AdaptiveSingularity::semi_infinite_upper(catalan, limits, error_bound, 1000)
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

    let catalan = util::Catalan11::new();
    let limits = catalan.lower;
    let error_bound = Tolerance::Absolute(TOL);
    let integral = AdaptiveSingularity::semi_infinite_upper(catalan, limits, error_bound, 1000)
        .unwrap()
        .integrate()
        .unwrap();
    let result = integral.result();
    let error = integral.error();
    let abs_actual_error = (EXPECTED - result).abs();
    let tol = TOL * result.abs();
    let iters = integral.iterations();
    assert_eq!(iters, 6);
    assert!(abs_actual_error < error);
    assert!(error < tol);

    //let catalan = util::ComplexCatalan12::new();
    //let limits = catalan.limits();
    //let error_bound = Tolerance::Absolute(TOL);
    //let integral = AdaptiveSingularity::finite(catalan, limits, error_bound, 1000)
    //    .unwrap()
    //    .integrate()
    //    .unwrap();
    //let result = integral.result();
    //let result_re = result.re;
    //let result_im = result.im;
    //let error = integral.error();
    //let abs_actual_error_re = (EXPECTED - result_re).abs();
    //let abs_actual_error_im = (EXPECTED - result_im).abs();
    //let tol = TOL * result.abs();
    //let iters = integral.iterations();
    //assert_eq!(iters, 38);
    //assert!(abs_actual_error_re < error);
    //assert!(abs_actual_error_im < error);
    //assert!(error < tol);

    Ok(())
}
