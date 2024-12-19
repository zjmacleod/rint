pub(crate) mod util;

#[macro_use]
mod macros {
    macro_rules! basic_test {
        (
            name: $func_name:ident,
            function: $type:ty,
            alpha: $alphaty:ty => $alpha:expr,
            rule: $rule:ident,
            lower: $lower:expr,
            upper: $upper:expr,
            exp_result: $exp_result:expr,
            exp_error: $exp_error:expr,
            exp_result_abs: $exp_result_abs:expr,
            exp_result_asc: $exp_result_asc:expr,
            abs_error_bound: $abs_error_bound:expr,
            rel_error_bound: $rel_error_bound:expr,
            $description:literal
        ) => {
            #[test]
            fn $func_name() -> Result<(), String> {
                use crate::quadrature::rule::Rule;
                use crate::quadrature::tests::util;
                use crate::quadrature::Integrator;
                let exp_result = $exp_result;
                let exp_error = $exp_error;
                let exp_result_abs = $exp_result_abs;
                let exp_result_asc = $exp_result_asc;

                let abs_error_bound = $abs_error_bound;
                let rel_error_bound = $rel_error_bound;
                let alpha: $alphaty = $alpha;

                let lower = $lower;
                let upper = $upper;

                let function = <$type>::new(alpha);
                let rule = Rule::$rule();

                let integral = Integrator::new(&function, &rule, Limits::new(lower, upper));

                let integral_result = integral.integrate();
                let result = integral_result.result();
                let error = integral_result.error();
                let result_abs = integral_result.result_abs();
                let result_asc = integral_result.result_asc();

                util::test_relative_error(
                    result,
                    exp_result,
                    abs_error_bound,
                    &format!("{} result", $description),
                )?;
                util::test_relative_error(
                    error,
                    exp_error,
                    rel_error_bound,
                    &format!("{} error", $description),
                )?;
                util::test_relative_error(
                    result_abs,
                    exp_result_abs,
                    abs_error_bound,
                    &format!("{} result_abs", $description),
                )?;
                util::test_relative_error(
                    result_asc,
                    exp_result_asc,
                    abs_error_bound,
                    &format!("{} result_asc", $description),
                )?;

                let lower = $upper;
                let upper = $lower;

                let integral = Integrator::new(&function, &rule, Limits::new(lower, upper));

                let integral_result = integral.integrate();
                let result = integral_result.result();
                let error = integral_result.error();
                let result_abs = integral_result.result_abs();
                let result_asc = integral_result.result_asc();

                util::test_relative_error(
                    result,
                    -exp_result,
                    abs_error_bound,
                    &format!("{} result reverse", $description),
                )?;
                util::test_relative_error(
                    error,
                    exp_error,
                    rel_error_bound,
                    &format!("{} error reverse", $description),
                )?;
                util::test_relative_error(
                    result_abs,
                    exp_result_abs,
                    abs_error_bound,
                    &format!("{} result_abs reverse", $description),
                )?;
                util::test_relative_error(
                    result_asc,
                    exp_result_asc,
                    abs_error_bound,
                    &format!("{} result_asc reverse", $description),
                )?;

                Ok(())
            }
        };
        (
            name: $func_name:ident,
            function: $type:ty,
            rule: $rule:ident,
            lower: $lower:expr,
            upper: $upper:expr,
            exp_result: $exp_result:expr,
            exp_error: $exp_error:expr,
            exp_result_abs: $exp_result_abs:expr,
            exp_result_asc: $exp_result_asc:expr,
            abs_error_bound: $abs_error_bound:expr,
            rel_error_bound: $rel_error_bound:expr
        ) => {
            #[test]
            fn $func_name() -> Result<(), String> {
                use crate::quadrature::rule::Rule;
                use crate::quadrature::tests::util;
                use crate::quadrature::Integrator;
                let exp_result = $exp_result;
                let exp_error = $exp_error;
                let exp_result_abs = $exp_result_abs;
                let exp_result_asc = $exp_result_asc;

                let abs_error_bound = $abs_error_bound;
                let rel_error_bound = $rel_error_bound;

                let lower = $lower;
                let upper = $upper;

                let function = <$type>::new();
                let rule = Rule::$rule();

                let integral = Integrator::new(Limits::new(lower, upper), &rule, &function);

                let integral_result = integral.integrate_internal();
                let result = integral_result.result();
                let error = integral_result.error();
                let result_abs = integral_result.result_abs();
                let result_asc = integral_result.result_asc();

                util::test_relative_error(
                    result,
                    exp_result,
                    abs_error_bound,
                    &format!("{} result", $description),
                )?;
                util::test_relative_error(
                    error,
                    exp_error,
                    rel_error_bound,
                    &format!("{} error", $description),
                )?;
                util::test_relative_error(
                    result_abs,
                    exp_result_abs,
                    abs_error_bound,
                    &format!("{} result_abs", $description),
                )?;
                util::test_relative_error(
                    result_asc,
                    exp_result_asc,
                    abs_error_bound,
                    &format!("{} result_asc", $description),
                )?;

                let lower = $upper;
                let upper = $lower;

                let integral = Integrator::new(Limits::new(lower, upper), &rule, &function);

                let integral_result = integral.integrate_internal();
                let result = integral_result.result();
                let error = integral_result.error();
                let result_abs = integral_result.result_abs();
                let result_asc = integral_result.result_asc();

                util::test_relative_error(
                    result,
                    -exp_result,
                    abs_error_bound,
                    &format!("{} result reverse", $description),
                )?;
                util::test_relative_error(
                    error,
                    exp_error,
                    rel_error_bound,
                    &format!("{} error reverse", $description),
                )?;
                util::test_relative_error(
                    result_abs,
                    exp_result_abs,
                    abs_error_bound,
                    &format!("{} result_abs reverse", $description),
                )?;
                util::test_relative_error(
                    result_asc,
                    exp_result_asc,
                    abs_error_bound,
                    &format!("{} result_asc reverse", $description),
                )?;

                Ok(())
            }
        };
    }
}

mod positive_function_with_singularity;
mod smooth_oscillating_function_unsymmetric_range;
mod smooth_positive_function;
