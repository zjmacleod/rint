// Test the basic Gauss-Kronrod integration rules with a positive function that
// has a singularity. This should give large values of the absolute error which
// would find discrepancies in the absolute error calculation.
//
// Ported from gsl-2.6/integration/test.c

use crate::Limits;

basic_test! {
        name: basic_gauss_kronrod_positive_function_with_singularity_15_point,
        function: util::Function1,
        alpha: f64 => -0.9,
        rule: gk15,
        lower: 0.0,
        upper: 1.0,
        exp_result:     1.555688196612745777E+01,
        exp_error:      2.350164577239293706E+01,
        exp_result_abs: 1.555688196612745777E+01,
        exp_result_asc: 2.350164577239293706E+01,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-7,
        "Function1 with singularity Basic 15-point"
}

basic_test! {
        name: basic_gauss_kronrod_positive_function_with_singularity_21_point,
        function: util::Function1,
        alpha: f64 => -0.9,
        rule: gk21,
        lower: 0.0,
        upper: 1.0,
        exp_result:     1.799045317938126232E+01,
        exp_error:      2.782360287710622515E+01,
        exp_result_abs: 1.799045317938126232E+01,
        exp_result_asc: 2.782360287710622515E+01,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-7,
        "Function1 with singularity Basic 21-point"
}

basic_test! {
        name: basic_gauss_kronrod_positive_function_with_singularity_31_point,
        function: util::Function1,
        alpha: f64 => -0.9,
        rule: gk31,
        lower: 0.0,
        upper: 1.0,
        exp_result:     2.081873305159121657E+01,
        exp_error:      3.296500137482590276E+01,
        exp_result_abs: 2.081873305159121301E+01,
        exp_result_asc: 3.296500137482590276E+01,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-7,
        "Function1 with singularity Basic 31-point"
}

basic_test! {
        name: basic_gauss_kronrod_positive_function_with_singularity_41_point,
        function: util::Function1,
        alpha: f64 => -0.9,
        rule: gk41,
        lower: 0.0,
        upper: 1.0,
        exp_result:     2.288677623903126701E+01,
        exp_error:      3.671538820274916048E+01,
        exp_result_abs: 2.288677623903126701E+01,
        exp_result_asc: 3.671538820274916048E+01,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-7,
        "Function1 with singularity Basic 41-point"
}

basic_test! {
        name: basic_gauss_kronrod_positive_function_with_singularity_51_point,
        function: util::Function1,
        alpha: f64 => -0.9,
        rule: gk51,
        lower: 0.0,
        upper: 1.0,
        exp_result:     2.449953612016972215E+01,
        exp_error:      3.967771249391228849E+01,
        exp_result_abs: 2.449953612016972215E+01,
        exp_result_asc: 3.967771249391228849E+01,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-7,
        "Function1 with singularity Basic 51-point"
}

basic_test! {
        name: basic_gauss_kronrod_positive_function_with_singularity_61_point,
        function: util::Function1,
        alpha: f64 => -0.9,
        rule: gk61,
        lower: 0.0,
        upper: 1.0,
        exp_result:     2.583030240976628988E+01,
        exp_error:      4.213750493076978643E+01,
        exp_result_abs: 2.583030240976628988E+01,
        exp_result_asc: 4.213750493076978643E+01,
        abs_error_bound: 1.0e-15,
        rel_error_bound: 1.0e-7,
        "Function1 with singularity Basic 61-point"
}
