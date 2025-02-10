use rint::multi::Adaptive;
use rint::multi::{Rule07, Rule09, Rule11, Rule13};
use rint::Tolerance;

mod util;

#[test]
fn van_booren_f1() {
    use util::multi::{F1, F1_TARGET};

    let max_iterations = 10000;
    const NDIM: usize = 6;
    const TARGET: f64 = F1_TARGET;
    println!("target:\t{TARGET}");

    {
        const TOL: f64 = 1e-2;
        println!("REL TOL:\t{TOL:e}");
        let function = F1::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F1::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-3;
        println!("REL TOL:\t{TOL:e}");
        let function = F1::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F1::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-4;
        println!("REL TOL:\t{TOL:e}");
        let function = F1::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F1::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-5;
        println!("REL TOL:\t{TOL:e}");
        let function = F1::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F1::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);
        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-6;
        println!("REL TOL:\t{TOL:e}");
        let function = F1::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F1::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-7;
        println!("REL TOL:\t{TOL:e}");
        let function = F1::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F1::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }
}

#[test]
fn van_booren_f2() {
    use util::multi::{F2, F2_TARGET};

    let max_iterations = 10000;
    const NDIM: usize = 4;
    const TARGET: f64 = F2_TARGET;
    println!("target:\t{TARGET}");

    {
        const TOL: f64 = 1e-2;
        println!("REL TOL:\t{TOL:e}");
        let function = F2::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F2::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-3;
        println!("REL TOL:\t{TOL:e}");
        let function = F2::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F2::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-4;
        println!("REL TOL:\t{TOL:e}");
        let function = F2::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F2::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-5;
        println!("REL TOL:\t{TOL:e}");
        let function = F2::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F2::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);
        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-6;
        println!("REL TOL:\t{TOL:e}");
        let function = F2::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F2::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-7;
        println!("REL TOL:\t{TOL:e}");
        let function = F2::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F2::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }
}

#[test]
fn van_booren_f3() {
    use util::multi::{F3, F3_TARGET};

    let max_iterations = 10000;
    const NDIM: usize = 3;
    const TARGET: f64 = F3_TARGET;
    println!("target:\t{TARGET}");

    {
        const TOL: f64 = 1e-2;
        println!("REL TOL:\t{TOL:e}");
        let function = F3::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F3::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F3::new();

        let limits = function.limits;
        let rule = Rule11::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule11\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-3;
        println!("REL TOL:\t{TOL:e}");
        let function = F3::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F3::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F3::new();

        let limits = function.limits;
        let rule = Rule11::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule11\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-4;
        println!("REL TOL:\t{TOL:e}");
        let function = F3::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F3::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F3::new();

        let limits = function.limits;
        let rule = Rule11::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule11\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-5;
        println!("REL TOL:\t{TOL:e}");
        let function = F3::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F3::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);
        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F3::new();

        let limits = function.limits;
        let rule = Rule11::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule11\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-6;
        println!("REL TOL:\t{TOL:e}");
        let function = F3::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F3::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F3::new();

        let limits = function.limits;
        let rule = Rule11::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule11\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-7;
        println!("REL TOL:\t{TOL:e}");
        let function = F3::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F3::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F3::new();

        let limits = function.limits;
        let rule = Rule11::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule11\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }
}

#[test]
fn van_booren_f4() {
    use util::multi::{F4, F4_TARGET};

    let max_iterations = 1000000;
    const NDIM: usize = 5;
    const TARGET: f64 = F4_TARGET;
    println!("target:\t{TARGET}");

    {
        const TOL: f64 = 1e-2;
        println!("REL TOL:\t{TOL:e}");
        let function = F4::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F4::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-3;
        println!("REL TOL:\t{TOL:e}");
        let function = F4::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F4::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-4;
        println!("REL TOL:\t{TOL:e}");
        let function = F4::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F4::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-5;
        println!("REL TOL:\t{TOL:e}");
        let function = F4::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F4::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);
        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-6;
        println!("REL TOL:\t{TOL:e}");
        let function = F4::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F4::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-7;
        println!("REL TOL:\t{TOL:e}");
        let function = F4::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F4::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }
}

#[test]
fn van_booren_f5() {
    use util::multi::{F5, F5_TARGET};

    let max_iterations = 10000;
    const NDIM: usize = 4;
    const TARGET: f64 = F5_TARGET;
    println!("target:\t{TARGET}");

    {
        const TOL: f64 = 1e-2;
        println!("REL TOL:\t{TOL:e}");
        let function = F5::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F5::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-3;
        println!("REL TOL:\t{TOL:e}");
        let function = F5::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F5::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-4;
        println!("REL TOL:\t{TOL:e}");
        let function = F5::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F5::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-5;
        println!("REL TOL:\t{TOL:e}");
        let function = F5::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F5::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);
        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-6;
        println!("REL TOL:\t{TOL:e}");
        let function = F5::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F5::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-7;
        println!("REL TOL:\t{TOL:e}");
        let function = F5::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F5::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }
}

#[test]
fn van_booren_f6() {
    use util::multi::{F6, F6_TARGET};

    let max_iterations = 10000;
    const NDIM: usize = 2;
    const TARGET: f64 = F6_TARGET;
    println!("target:\t{TARGET}");

    {
        const TOL: f64 = 1e-2;
        println!("REL TOL:\t{TOL:e}");
        let function = F6::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F6::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-3;
        println!("REL TOL:\t{TOL:e}");
        let function = F6::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F6::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-4;
        println!("REL TOL:\t{TOL:e}");
        let function = F6::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F6::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-5;
        println!("REL TOL:\t{TOL:e}");
        let function = F6::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
        let function = F6::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-6;
        println!("REL TOL:\t{TOL:e}");
        let function = F6::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F6::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-7;
        println!("REL TOL:\t{TOL:e}");
        let function = F6::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F6::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }
}

#[test]
fn van_booren_f7() {
    use util::multi::{F7, F7_TARGET};

    let max_iterations = 10000;
    const NDIM: usize = 3;
    const TARGET: f64 = F7_TARGET;
    println!("target:\t{TARGET}");

    {
        const TOL: f64 = 1e-2;
        println!("REL TOL:\t{TOL:e}");
        let function = F7::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F7::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F7::new();

        let limits = function.limits;
        let rule = Rule11::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule11\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-3;
        println!("REL TOL:\t{TOL:e}");
        let function = F7::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F7::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F7::new();

        let limits = function.limits;
        let rule = Rule11::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule11\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-4;
        println!("REL TOL:\t{TOL:e}");
        let function = F7::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F7::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F7::new();

        let limits = function.limits;
        let rule = Rule11::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule11\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-5;
        println!("REL TOL:\t{TOL:e}");
        let function = F7::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F7::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);
        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F7::new();

        let limits = function.limits;
        let rule = Rule11::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule11\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-6;
        println!("REL TOL:\t{TOL:e}");
        let function = F7::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F7::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F7::new();

        let limits = function.limits;
        let rule = Rule11::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule11\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-7;
        println!("REL TOL:\t{TOL:e}");
        let function = F7::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F7::new();

        let limits = function.limits;
        let rule = Rule09::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule09\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F7::new();

        let limits = function.limits;
        let rule = Rule11::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule11\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }
}

#[test]
fn van_booren_f8() {
    use util::multi::{F8, F8_TARGET};

    let max_iterations = 10000;
    const NDIM: usize = 2;
    const TARGET: f64 = F8_TARGET;
    println!("target:\t{TARGET}");

    {
        const TOL: f64 = 1e-2;
        println!("REL TOL:\t{TOL:e}");
        let function = F8::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F8::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-3;
        println!("REL TOL:\t{TOL:e}");
        let function = F8::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F8::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-4;
        println!("REL TOL:\t{TOL:e}");
        let function = F8::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F8::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-5;
        println!("REL TOL:\t{TOL:e}");
        let function = F8::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
        let function = F8::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-6;
        println!("REL TOL:\t{TOL:e}");
        let function = F8::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F8::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-7;
        println!("REL TOL:\t{TOL:e}");
        let function = F8::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F8::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }
}

#[test]
fn van_booren_f9() {
    use util::multi::{F9, F9_TARGET};

    let max_iterations = 10000;
    const NDIM: usize = 2;
    const TARGET: f64 = F9_TARGET;
    println!("target:\t{TARGET}");

    {
        const TOL: f64 = 1e-2;
        println!("REL TOL:\t{TOL:e}");
        let function = F9::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F9::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-3;
        println!("REL TOL:\t{TOL:e}");
        let function = F9::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F9::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-4;
        println!("REL TOL:\t{TOL:e}");
        let function = F9::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F9::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-5;
        println!("REL TOL:\t{TOL:e}");
        let function = F9::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
        let function = F9::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-6;
        println!("REL TOL:\t{TOL:e}");
        let function = F9::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F9::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-7;
        println!("REL TOL:\t{TOL:e}");
        let function = F9::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F9::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }
}

#[test]
fn van_booren_f10() {
    use util::multi::{F10, F10_TARGET};

    let max_iterations = 10000;
    const NDIM: usize = 2;
    const TARGET: f64 = F10_TARGET;
    println!("target:\t{TARGET}");

    {
        const TOL: f64 = 1e-2;
        println!("REL TOL:\t{TOL:e}");
        let function = F10::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F10::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-3;
        println!("REL TOL:\t{TOL:e}");
        let function = F10::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F10::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-4;
        println!("REL TOL:\t{TOL:e}");
        let function = F10::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F10::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-5;
        println!("REL TOL:\t{TOL:e}");
        let function = F10::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
        let function = F10::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-6;
        println!("REL TOL:\t{TOL:e}");
        let function = F10::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F10::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }

    {
        const TOL: f64 = 1e-7;
        println!("REL TOL:\t{TOL:e}");
        let function = F10::new();

        let limits = function.limits;
        let rule = Rule07::<NDIM>::generate().unwrap();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule07\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");

        let function = F10::new();

        let limits = function.limits;
        let rule = Rule13::generate();
        let tolerance = Tolerance::Relative(TOL);

        let integrand = Adaptive::new(&function, &rule, limits, tolerance, max_iterations).unwrap();

        let integral = integrand.integrate().unwrap();
        let result = integral.result();
        let error = integral.error();
        let iterations = integral.iterations();
        let evaluations = integral.evaluations();

        let actual_error = (result - TARGET).abs();
        let requested_error = TOL * result.abs();

        assert!(actual_error < requested_error);

        println!("rule:\tRule13\tresult:\t{result:e}\terror:\t{error:e}\tactual:\t{actual_error:e}\titerations:\t{iterations}\tevaluations:\t{evaluations}");
    }
}
