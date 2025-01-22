//use rint::quadrature::adaptive::Kind;
use rint::quadrature::AdaptiveSingularity;
use rint::{Limits, Tolerance};

mod util;

// Test the smooth Function1 with 15 point adaptive integration and relative error
// bound.
singularity_test! {
        name: singularity_gauss_kronrod_smooth_function_relative_error_bound,
        function: util::Function1,
        alpha: f64 => 2.6,
        lower: 0.0,
        upper: 1.0,
        exp_result:      7.716049382715789440E-02,
        exp_error:       2.216394961010438404E-12,
        exp_iterations:  5,
        exp_evaluations: 189,
        tolerance_rel: 1e-10,
        test_abs_error_bound: 1e-15,
        test_rel_error_bound: 1e-6,
        "Function1 AdaptiveSingularity::general relative bound"
}

// XXX Test passes but result for error estimate 2e-14 off from GSL.
singularity_test! {
        name: singularity_gauss_kronrod_function11_absolute_error_bound,
        function: util::Function11,
        alpha: f64 => 2.0,
        lower: 1.0,
        upper: 1000.0,
        exp_result:     -5.908755278982136588E+03,
        exp_error:       1.299646281053874554E-10,
        exp_iterations:  9,
        exp_evaluations: 357,
        tolerance_abs: 1e-7,
        test_abs_error_bound: 1e-15,
        test_rel_error_bound: 1e-3,
        "Function11 AdaptiveSingularity::general absolute bound"
}

/*Test infinite range integral f455 using a relative error bound */
singularity_test! {
        name: singularity_gauss_kronrod_function455_relative_error_bound_semi_infinite_positive,
        function: util::Function455,
        lower: 0.0,
        exp_result:     -3.616892186127022568E-01,
        exp_error:       3.016716913328831851E-06,
        exp_iterations:  10,
        exp_evaluations: 570,
        tolerance_rel: 1e-3,
        test_abs_error_bound: 1e-14,
        test_rel_error_bound: 1e-5,
        "Function455 AdaptiveSingularity::semi_infinite_positive relative bound"
}

/*Test infinite range integral f15 using a relative error bound */
singularity_test! {
        name: singularity_gauss_kronrod_function15_relative_error_bound,
        function: util::Function15,
        alpha: i32 => 5,
        lower: 0.0,
        exp_result:      6.553600000000024738E+04,
        exp_error:       7.121667111456009280E-04,
        exp_iterations:  10,
        exp_evaluations: 570,
        tolerance_rel: 1e-7,
        test_abs_error_bound: 1e-14,
        test_rel_error_bound: 1e-5,
        "Function15 AdaptiveSingularity::semi_infinite_positive relative bound"
}

/*Test infinite range integral f16 using an absolute error bound */
singularity_test! {
        name: singularity_gauss_kronrod_function16_absolute_error_bound,
        function: util::Function16,
        alpha: i32 => 1,
        lower: 99.9,
        exp_result:      1.000000000006713292E-04,
        exp_error:       3.084062020905636316E-09,
        exp_iterations:  6,
        exp_evaluations: 330,
        tolerance_abs: 1e-7,
        test_abs_error_bound: 1e-14,
        test_rel_error_bound: 1e-5,
        "Function16 AdaptiveSingularity::semi_infinite_positive absolute bound"
}

/*Test infinite range integral myfn1 using a absolute error bound */
singularity_test! {
        name: singularity_gauss_kronrod_myfunction1_absolute_error_bound_semi_infinite_positive,
        function: util::MyFunciton1,
        exp_result:      2.275875794468747770E+00,
        exp_error:       7.436490118267390744E-09,
        exp_iterations:  5,
        exp_evaluations: 270,
        tolerance_abs: 1e-7,
        test_abs_error_bound: 1e-14,
        test_rel_error_bound: 1e-5,
        "MyFunction1 AdaptiveSingularity::infinite absolute bound"
}

/*Test infinite range integral myfn2 using a absolute error bound lower */
singularity_test! {
        name: singularity_gauss_kronrod_myfunction2_absolute_error_bound_semi_infinite_positive,
        function: util::MyFunciton2,
        alpha: f64 => 1.0,
        upper: 1.0,
        exp_result:      2.718281828459044647E+00,
        exp_error:       1.588185109253204805E-10,
        exp_iterations:  5,
        exp_evaluations: 270,
        tolerance_abs: 1e-7,
        test_abs_error_bound: 1e-14,
        test_rel_error_bound: 1e-5,
        "MyFunction2 AdaptiveSingularity::semi_infinite_negative absolute bound"
}
