//#[macro_export]
//macro_rules! basic_test {
//    (
//        name: $func_name:ident,
//        function: $type:ty,
//        alpha: $alphaty:ty => $alpha:expr,
//        rule: $rule:ident,
//        lower: $lower:expr,
//        upper: $upper:expr,
//        exp_result: $exp_result:expr,
//        exp_error: $exp_error:expr,
//        exp_result_abs: $exp_result_abs:expr,
//        exp_result_asc: $exp_result_asc:expr,
//        abs_error_bound: $abs_error_bound:expr,
//        rel_error_bound: $rel_error_bound:expr,
//        $description:literal
//    ) => {
//        #[test]
//        fn $func_name() -> Result<(), String> {
//            let exp_result = $exp_result;
//            let exp_error = $exp_error;
//            let exp_result_abs = $exp_result_abs;
//            let exp_result_asc = $exp_result_asc;
//
//            let abs_error_bound = $abs_error_bound;
//            let rel_error_bound = $rel_error_bound;
//            let alpha: $alphaty = $alpha;
//
//            let lower = $lower;
//            let upper = $upper;
//
//            let function = <$type>::new(alpha);
//            let rule = Rule::$rule();
//
//            let integral = Basic::new(Limits::new(lower, upper), rule, &function);
//
//            let integral_result = integral.integrate();
//            let result = integral_result.result();
//            let error = integral_result.error();
//            let result_abs = integral_result.result_abs();
//            let result_asc = integral_result.result_asc();
//
//            util::test_relative_error(
//                result,
//                exp_result,
//                abs_error_bound,
//                &format!("{} result", $description),
//            )?;
//            util::test_relative_error(
//                error,
//                exp_error,
//                rel_error_bound,
//                &format!("{} error", $description),
//            )?;
//            util::test_relative_error(
//                result_abs,
//                exp_result_abs,
//                abs_error_bound,
//                &format!("{} result_abs", $description),
//            )?;
//            util::test_relative_error(
//                result_asc,
//                exp_result_asc,
//                abs_error_bound,
//                &format!("{} result_asc", $description),
//            )?;
//
//            let lower = $upper;
//            let upper = $lower;
//
//            let integral = Basic::new(Limits::new(lower, upper), rule, &function);
//
//            let integral_result = integral.integrate();
//            let result = integral_result.result();
//            let error = integral_result.error();
//            let result_abs = integral_result.result_abs();
//            let result_asc = integral_result.result_asc();
//
//            util::test_relative_error(
//                result,
//                -exp_result,
//                abs_error_bound,
//                &format!("{} result reverse", $description),
//            )?;
//            util::test_relative_error(
//                error,
//                exp_error,
//                rel_error_bound,
//                &format!("{} error reverse", $description),
//            )?;
//            util::test_relative_error(
//                result_abs,
//                exp_result_abs,
//                abs_error_bound,
//                &format!("{} result_abs reverse", $description),
//            )?;
//            util::test_relative_error(
//                result_asc,
//                exp_result_asc,
//                abs_error_bound,
//                &format!("{} result_asc reverse", $description),
//            )?;
//
//            Ok(())
//        }
//    };
//    (
//        name: $func_name:ident,
//        function: $type:ty,
//        rule: $rule:ident,
//        lower: $lower:expr,
//        upper: $upper:expr,
//        exp_result: $exp_result:expr,
//        exp_error: $exp_error:expr,
//        exp_result_abs: $exp_result_abs:expr,
//        exp_result_asc: $exp_result_asc:expr,
//        abs_error_bound: $abs_error_bound:expr,
//        rel_error_bound: $rel_error_bound:expr
//    ) => {
//        #[test]
//        fn $func_name() -> Result<(), String> {
//            let exp_result = $exp_result;
//            let exp_error = $exp_error;
//            let exp_result_abs = $exp_result_abs;
//            let exp_result_asc = $exp_result_asc;
//
//            let abs_error_bound = $abs_error_bound;
//            let rel_error_bound = $rel_error_bound;
//
//            let lower = $lower;
//            let upper = $upper;
//
//            let function = <$type>::new();
//            let rule = Rule::$rule();
//
//            let integral = Basic::new(Limits::new(lower, upper), rule, &function);
//
//            let integral_result = integral.integrate();
//            let result = integral_result.result();
//            let error = integral_result.error();
//            let result_abs = integral_result.result_abs();
//            let result_asc = integral_result.result_asc();
//
//            util::test_relative_error(
//                result,
//                exp_result,
//                abs_error_bound,
//                &format!("{} result", $description),
//            )?;
//            util::test_relative_error(
//                error,
//                exp_error,
//                rel_error_bound,
//                &format!("{} error", $description),
//            )?;
//            util::test_relative_error(
//                result_abs,
//                exp_result_abs,
//                abs_error_bound,
//                &format!("{} result_abs", $description),
//            )?;
//            util::test_relative_error(
//                result_asc,
//                exp_result_asc,
//                abs_error_bound,
//                &format!("{} result_asc", $description),
//            )?;
//
//            let lower = $upper;
//            let upper = $lower;
//
//            let integral = Basic::new(Limits::new(lower, upper), rule, &function);
//
//            let integral_result = integral.integrate();
//            let result = integral_result.result();
//            let error = integral_result.error();
//            let result_abs = integral_result.result_abs();
//            let result_asc = integral_result.result_asc();
//
//            util::test_relative_error(
//                result,
//                -exp_result,
//                abs_error_bound,
//                &format!("{} result reverse", $description),
//            )?;
//            util::test_relative_error(
//                error,
//                exp_error,
//                rel_error_bound,
//                &format!("{} error reverse", $description),
//            )?;
//            util::test_relative_error(
//                result_abs,
//                exp_result_abs,
//                abs_error_bound,
//                &format!("{} result_abs reverse", $description),
//            )?;
//            util::test_relative_error(
//                result_asc,
//                exp_result_asc,
//                abs_error_bound,
//                &format!("{} result_asc reverse", $description),
//            )?;
//
//            Ok(())
//        }
//    };
//}

#[macro_export]
macro_rules! adaptive_test_passing {
    (
        name: $func_name:ident,
        function: $type:ty,
        alpha: $alphaty:ty => $alpha:expr,
        rule: $rule:ident,
        lower: $lower:expr,
        upper: $upper:expr,
        exp_result:      $exp_result:expr,
        exp_error:       $exp_error:expr,
        exp_iterations:  $exp_iterations:expr,
        exp_evaluations: $exp_evaluations:expr,
        tolerance_rel: $tolerance:expr,
        test_abs_error_bound: $abs_error_bound:expr,
        test_rel_error_bound: $rel_error_bound:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_error_bound = $abs_error_bound;
            let rel_error_bound = $rel_error_bound;
            let rule = Rule::$rule();
            let error_bound = ErrorBound::Relative($tolerance);
            let alpha: $alphaty = $alpha;

            let lower = $lower;
            let upper = $upper;

            let function = <$type>::new(alpha);

            let integral = Adaptive::new(
                &function,
                rule,
                Limits::new(lower, upper),
                error_bound,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations", $description),
            )?;

            let rule = Rule::$rule();
            let error_bound = ErrorBound::Relative($tolerance);
            let lower = $upper;
            let upper = $lower;

            let integral = Adaptive::new(
                &function,
                rule,
                Limits::new(lower, upper),
                error_bound,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations reverse", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations reverse", $description),
            )?;

            Ok(())
        }
    };
    (
        name: $func_name:ident,
        function: $type:ty,
        alpha: $alphaty:ty => $alpha:expr,
        rule: $rule:ident,
        lower: $lower:expr,
        upper: $upper:expr,
        exp_result:      $exp_result:expr,
        exp_error:       $exp_error:expr,
        exp_iterations:  $exp_iterations:expr,
        exp_evaluations: $exp_evaluations:expr,
        tolerance_abs: $tolerance:expr,
        test_abs_error_bound: $abs_error_bound:expr,
        test_rel_error_bound: $rel_error_bound:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_error_bound = $abs_error_bound;
            let rel_error_bound = $rel_error_bound;
            let rule = Rule::$rule();
            let error_bound = ErrorBound::Absolute($tolerance);
            let alpha: $alphaty = $alpha;

            let lower = $lower;
            let upper = $upper;

            let function = <$type>::new(alpha);

            let integral = Adaptive::new(
                &function,
                rule,
                Limits::new(lower, upper),
                error_bound,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations", $description),
            )?;

            let rule = Rule::$rule();
            let error_bound = ErrorBound::Absolute($tolerance);
            let lower = $upper;
            let upper = $lower;

            let integral = Adaptive::new(
                &function,
                rule,
                Limits::new(lower, upper),
                error_bound,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations reverse", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations reverse", $description),
            )?;

            Ok(())
        }
    };
}

#[macro_export]
macro_rules! adaptive_test_error {
    (
        name: $func_name:ident,
        function: $type:ty,
        alpha: $alphaty:ty => $alpha:expr,
        rule: $rule:ident,
        lower: $lower:expr,
        upper: $upper:expr,
        exp_result:      $exp_result:expr,
        exp_error:       $exp_error:expr,
        exp_iterations:  $exp_iterations:expr,
        exp_evaluations: $exp_evaluations:expr,
        tolerance_abs: $tolerance:expr,
        test_abs_error_bound: $abs_error_bound:expr,
        test_rel_error_bound: $rel_error_bound:expr,
        kind: $kind:ident,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_error_bound = $abs_error_bound;
            let rel_error_bound = $rel_error_bound;
            let rule = Rule::$rule();
            let error_bound = ErrorBound::Absolute($tolerance);
            let alpha: $alphaty = $alpha;

            let lower = $lower;
            let upper = $upper;

            let function = <$type>::new(alpha);

            let integral = Adaptive::new(
                &function,
                rule,
                Limits::new(lower, upper),
                error_bound,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate();

            if let Err(err) = integral_result {
                if let Kind::$kind = err.kind() {
                    let result = err.result();
                    let error = err.error();
                    let iterations = err.iterations();
                    let evaluations = err.function_evaluations();

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
                    util::test_int(
                        iterations,
                        exp_iterations,
                        &format!("{} iterations", $description),
                    )?;
                    util::test_int(
                        evaluations,
                        exp_evaluations,
                        &format!("{} evaluations", $description),
                    )?;
                } else {
                    return Err(format!(
                        "Wrong error type. Expected {} but recieved {:?}",
                        stringify!($kind),
                        err.kind()
                    ));
                }
            } else {
                return Err(String::from(
                    "Integration returned Ok, but should have been an error.",
                ));
            }

            let rule = Rule::$rule();
            let error_bound = ErrorBound::Absolute($tolerance);
            let lower = $upper;
            let upper = $lower;

            let integral = Adaptive::new(
                &function,
                rule,
                Limits::new(lower, upper),
                error_bound,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate();

            if let Err(err) = integral_result {
                if let Kind::$kind = err.kind() {
                    let result = err.result();
                    let error = err.error();
                    let iterations = err.iterations();
                    let evaluations = err.function_evaluations();

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
                    util::test_int(
                        iterations,
                        exp_iterations,
                        &format!("{} iterations reverse", $description),
                    )?;
                    util::test_int(
                        evaluations,
                        exp_evaluations,
                        &format!("{} evaluations reverse", $description),
                    )?;
                } else {
                    return Err(format!(
                        "Wrong error type. Expected {} but recieved {:?}",
                        stringify!($kind),
                        err.kind()
                    ));
                }
            } else {
                return Err(String::from(
                    "Integration returned Ok, but should have been an error.",
                ));
            }

            Ok(())
        }
    };
    (
        name: $func_name:ident,
        function: $type:ty,
        alpha: $alphaty:ty => $alpha:expr,
        rule: $rule:ident,
        lower: $lower:expr,
        upper: $upper:expr,
        exp_result:      $exp_result:expr,
        exp_error:       $exp_error:expr,
        exp_iterations:  $exp_iterations:expr,
        exp_evaluations: $exp_evaluations:expr,
        tolerance_rel: $tolerance:expr,
        test_abs_error_bound: $abs_error_bound:expr,
        test_rel_error_bound: $rel_error_bound:expr,
        kind: $kind:ident,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_error_bound = $abs_error_bound;
            let rel_error_bound = $rel_error_bound;
            let rule = Rule::$rule();
            let error_bound = ErrorBound::Relative($tolerance);
            let alpha: $alphaty = $alpha;

            let lower = $lower;
            let upper = $upper;

            let function = <$type>::new(alpha);

            let integral = Adaptive::new(
                &function,
                rule,
                Limits::new(lower, upper),
                error_bound,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate();

            if let Err(err) = integral_result {
                if let Kind::$kind = err.kind() {
                    let result = err.result();
                    let error = err.error();
                    let iterations = err.iterations();
                    let evaluations = err.function_evaluations();

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
                    util::test_int(
                        iterations,
                        exp_iterations,
                        &format!("{} iterations", $description),
                    )?;
                    util::test_int(
                        evaluations,
                        exp_evaluations,
                        &format!("{} evaluations", $description),
                    )?;
                } else {
                    return Err(format!(
                        "Wrong error type. Expected {} but recieved {:?}",
                        stringify!($kind),
                        err.kind()
                    ));
                }
            } else {
                return Err(String::from(
                    "Integration returned Ok, but should have been an error.",
                ));
            }

            let error_bound = ErrorBound::Relative($tolerance);
            let lower = $upper;
            let upper = $lower;

            let integral = Adaptive::new(
                &function,
                rule,
                Limits::new(lower, upper),
                error_bound,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate();

            if let Err(err) = integral_result {
                if let Kind::$kind = err.kind() {
                    let result = err.result();
                    let error = err.error();
                    let iterations = err.iterations();
                    let evaluations = err.function_evaluations();

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
                    util::test_int(
                        iterations,
                        exp_iterations,
                        &format!("{} iterations reverse", $description),
                    )?;
                    util::test_int(
                        evaluations,
                        exp_evaluations,
                        &format!("{} evaluations reverse", $description),
                    )?;
                } else {
                    return Err(format!(
                        "Wrong error type. Expected {} but recieved {:?}",
                        stringify!($kind),
                        err.kind()
                    ));
                }
            } else {
                return Err(String::from(
                    "Integration returned Ok, but should have been an error.",
                ));
            }

            Ok(())
        }
    };
    (
        name: $func_name:ident,
        function: $type:ty,
        alpha: $alphaty:ty => $alpha:expr,
        rule: $rule:ident,
        lower: $lower:expr,
        upper: $upper:expr,
        iterations: $iter:expr,
        exp_result:      $exp_result:expr,
        exp_error:       $exp_error:expr,
        exp_iterations:  $exp_iterations:expr,
        exp_evaluations: $exp_evaluations:expr,
        tolerance_abs: $tolerance:expr,
        test_abs_error_bound: $abs_error_bound:expr,
        test_rel_error_bound: $rel_error_bound:expr,
        kind: $kind:ident,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_error_bound = $abs_error_bound;
            let rel_error_bound = $rel_error_bound;
            let rule = Rule::$rule();
            let error_bound = ErrorBound::Absolute($tolerance);
            let alpha: $alphaty = $alpha;

            let lower = $lower;
            let upper = $upper;

            let function = <$type>::new(alpha);

            let integral = Adaptive::new(
                &function,
                rule,
                Limits::new(lower, upper),
                error_bound,
                $iter,
            )
            .unwrap();

            let integral_result = integral.integrate();

            if let Err(err) = integral_result {
                if let Kind::$kind = err.kind() {
                    let result = err.result();
                    let error = err.error();
                    let iterations = err.iterations();
                    let evaluations = err.function_evaluations();

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
                    util::test_int(
                        iterations,
                        exp_iterations,
                        &format!("{} iterations", $description),
                    )?;
                    util::test_int(
                        evaluations,
                        exp_evaluations,
                        &format!("{} evaluations", $description),
                    )?;
                } else {
                    return Err(format!(
                        "Wrong error type. Expected {} but recieved {:?}",
                        stringify!($kind),
                        err.kind()
                    ));
                }
            } else {
                return Err(String::from(
                    "Integration returned Ok, but should have been an error.",
                ));
            }

            let rule = Rule::$rule();
            let error_bound = ErrorBound::Absolute($tolerance);
            let lower = $upper;
            let upper = $lower;

            let integral = Adaptive::new(
                &function,
                rule,
                Limits::new(lower, upper),
                error_bound,
                $iter,
            )
            .unwrap();

            let integral_result = integral.integrate();

            if let Err(err) = integral_result {
                if let Kind::$kind = err.kind() {
                    let result = err.result();
                    let error = err.error();
                    let iterations = err.iterations();
                    let evaluations = err.function_evaluations();

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
                    util::test_int(
                        iterations,
                        exp_iterations,
                        &format!("{} iterations reverse", $description),
                    )?;
                    util::test_int(
                        evaluations,
                        exp_evaluations,
                        &format!("{} evaluations reverse", $description),
                    )?;
                } else {
                    return Err(format!(
                        "Wrong error type. Expected {} but recieved {:?}",
                        stringify!($kind),
                        err.kind()
                    ));
                }
            } else {
                return Err(String::from(
                    "Integration returned Ok, but should have been an error.",
                ));
            }

            Ok(())
        }
    };
    (
        name: $func_name:ident,
        function: $type:ty,
        alpha: $alphaty:ty => $alpha:expr,
        rule: $rule:ident,
        lower: $lower:expr,
        upper: $upper:expr,
        iterations: $iter:expr,
        exp_result:      $exp_result:expr,
        exp_error:       $exp_error:expr,
        exp_iterations:  $exp_iterations:expr,
        exp_evaluations: $exp_evaluations:expr,
        tolerance_rel: $tolerance:expr,
        test_abs_error_bound: $abs_error_bound:expr,
        test_rel_error_bound: $rel_error_bound:expr,
        kind: $kind:ident,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_error_bound = $abs_error_bound;
            let rel_error_bound = $rel_error_bound;
            let rule = Rule::$rule();
            let error_bound = ErrorBound::Relative($tolerance);
            let alpha: $alphaty = $alpha;

            let lower = $lower;
            let upper = $upper;

            let function = <$type>::new(alpha);

            let integral = Adaptive::new(
                &function,
                rule,
                Limits::new(lower, upper),
                error_bound,
                $iter,
            )
            .unwrap();

            let integral_result = integral.integrate();

            if let Err(err) = integral_result {
                if let Kind::$kind = err.kind() {
                    let result = err.result();
                    let error = err.error();
                    let iterations = err.iterations();
                    let evaluations = err.function_evaluations();

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
                    util::test_int(
                        iterations,
                        exp_iterations,
                        &format!("{} iterations", $description),
                    )?;
                    util::test_int(
                        evaluations,
                        exp_evaluations,
                        &format!("{} evaluations", $description),
                    )?;
                } else {
                    return Err(format!(
                        "Wrong error type. Expected {} but recieved {:?}",
                        stringify!($kind),
                        err.kind()
                    ));
                }
            } else {
                return Err(String::from(
                    "Integration returned Ok, but should have been an error.",
                ));
            }

            let error_bound = ErrorBound::Relative($tolerance);
            let lower = $upper;
            let upper = $lower;

            let integral = Adaptive::new(
                &function,
                rule,
                Limits::new(lower, upper),
                error_bound,
                $iter,
            )
            .unwrap();

            let integral_result = integral.integrate();

            if let Err(err) = integral_result {
                if let Kind::$kind = err.kind() {
                    let result = err.result();
                    let error = err.error();
                    let iterations = err.iterations();
                    let evaluations = err.function_evaluations();

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
                    util::test_int(
                        iterations,
                        exp_iterations,
                        &format!("{} iterations reverse", $description),
                    )?;
                    util::test_int(
                        evaluations,
                        exp_evaluations,
                        &format!("{} evaluations reverse", $description),
                    )?;
                } else {
                    return Err(format!(
                        "Wrong error type. Expected {} but recieved {:?}",
                        stringify!($kind),
                        err.kind()
                    ));
                }
            } else {
                return Err(String::from(
                    "Integration returned Ok, but should have been an error.",
                ));
            }

            Ok(())
        }
    };
}

#[macro_export]
macro_rules! singularity_test {
    (
        name: $func_name:ident,
        function: $type:ty,
        alpha: $alphaty:ty => $alpha:expr,
        lower: $lower:expr,
        upper: $upper:expr,
        exp_result:      $exp_result:expr,
        exp_error:       $exp_error:expr,
        exp_iterations:  $exp_iterations:expr,
        exp_evaluations: $exp_evaluations:expr,
        tolerance_rel: $tolerance:expr,
        test_abs_error_bound: $abs_error_bound:expr,
        test_rel_error_bound: $rel_error_bound:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_error_bound = $abs_error_bound;
            let rel_error_bound = $rel_error_bound;
            let error_bound = ErrorBound::Relative($tolerance);
            let alpha: $alphaty = $alpha;

            let lower = $lower;
            let upper = $upper;

            let function = <$type>::new(alpha);

            let integral = AdaptiveSingularity::general(
                &function,
                Limits::new(lower, upper),
                error_bound,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations", $description),
            )?;

            let lower = $upper;
            let upper = $lower;
            let error_bound = ErrorBound::Relative($tolerance);

            let integral = AdaptiveSingularity::general(
                &function,
                Limits::new(lower, upper),
                error_bound,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations reverse", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations reverse", $description),
            )?;

            Ok(())
        }
    };
    (
        name: $func_name:ident,
        function: $type:ty,
        alpha: $alphaty:ty => $alpha:expr,
        lower: $lower:expr,
        upper: $upper:expr,
        exp_result:      $exp_result:expr,
        exp_error:       $exp_error:expr,
        exp_iterations:  $exp_iterations:expr,
        exp_evaluations: $exp_evaluations:expr,
        tolerance_abs: $tolerance:expr,
        test_abs_error_bound: $abs_error_bound:expr,
        test_rel_error_bound: $rel_error_bound:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_error_bound = $abs_error_bound;
            let rel_error_bound = $rel_error_bound;
            let error_bound = ErrorBound::Absolute($tolerance);
            let alpha: $alphaty = $alpha;

            let lower = $lower;
            let upper = $upper;

            let function = <$type>::new(alpha);

            let integral = AdaptiveSingularity::general(
                &function,
                Limits::new(lower, upper),
                error_bound,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations", $description),
            )?;

            let lower = $upper;
            let upper = $lower;
            let error_bound = ErrorBound::Absolute($tolerance);

            let integral = AdaptiveSingularity::general(
                &function,
                Limits::new(lower, upper),
                error_bound,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations reverse", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations reverse", $description),
            )?;

            Ok(())
        }
    };
    (
        name: $func_name:ident,
        function: $type:ty,
        alpha: $alphaty:ty => $alpha:expr,
        upper: $upper:expr,
        exp_result:      $exp_result:expr,
        exp_error:       $exp_error:expr,
        exp_iterations:  $exp_iterations:expr,
        exp_evaluations: $exp_evaluations:expr,
        tolerance_rel: $tolerance:expr,
        test_abs_error_bound: $abs_error_bound:expr,
        test_rel_error_bound: $rel_error_bound:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_error_bound = $abs_error_bound;
            let rel_error_bound = $rel_error_bound;
            let alpha: $alphaty = $alpha;
            let error_bound = ErrorBound::Relative($tolerance);

            let upper = $upper;

            let function = <$type>::new(alpha);

            let integral =
                AdaptiveSingularity::semi_infinite_negative(&function, upper, error_bound, 1000)
                    .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations", $description),
            )?;

            Ok(())
        }
    };
    (
        name: $func_name:ident,
        function: $type:ty,
        alpha: $alphaty:ty => $alpha:expr,
        upper: $upper:expr,
        exp_result:      $exp_result:expr,
        exp_error:       $exp_error:expr,
        exp_iterations:  $exp_iterations:expr,
        exp_evaluations: $exp_evaluations:expr,
        tolerance_abs: $tolerance:expr,
        test_abs_error_bound: $abs_error_bound:expr,
        test_rel_error_bound: $rel_error_bound:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_error_bound = $abs_error_bound;
            let rel_error_bound = $rel_error_bound;
            let alpha: $alphaty = $alpha;
            let error_bound = ErrorBound::Absolute($tolerance);

            let upper = $upper;

            let function = <$type>::new(alpha);

            let integral =
                AdaptiveSingularity::semi_infinite_negative(&function, upper, error_bound, 1000)
                    .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations", $description),
            )?;

            Ok(())
        }
    };
    (
        name: $func_name:ident,
        function: $type:ty,
        alpha: $alphaty:ty => $alpha:expr,
        lower: $lower:expr,
        exp_result:      $exp_result:expr,
        exp_error:       $exp_error:expr,
        exp_iterations:  $exp_iterations:expr,
        exp_evaluations: $exp_evaluations:expr,
        tolerance_rel: $tolerance:expr,
        test_abs_error_bound: $abs_error_bound:expr,
        test_rel_error_bound: $rel_error_bound:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_error_bound = $abs_error_bound;
            let rel_error_bound = $rel_error_bound;
            let alpha: $alphaty = $alpha;
            let error_bound = ErrorBound::Relative($tolerance);

            let lower = $lower;

            let function = <$type>::new(alpha);

            let integral =
                AdaptiveSingularity::semi_infinite_positive(&function, lower, error_bound, 1000)
                    .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations", $description),
            )?;

            Ok(())
        }
    };
    (
        name: $func_name:ident,
        function: $type:ty,
        alpha: $alphaty:ty => $alpha:expr,
        lower: $lower:expr,
        exp_result:      $exp_result:expr,
        exp_error:       $exp_error:expr,
        exp_iterations:  $exp_iterations:expr,
        exp_evaluations: $exp_evaluations:expr,
        tolerance_abs: $tolerance:expr,
        test_abs_error_bound: $abs_error_bound:expr,
        test_rel_error_bound: $rel_error_bound:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_error_bound = $abs_error_bound;
            let rel_error_bound = $rel_error_bound;
            let alpha: $alphaty = $alpha;
            let error_bound = ErrorBound::Absolute($tolerance);

            let lower = $lower;

            let function = <$type>::new(alpha);

            let integral =
                AdaptiveSingularity::semi_infinite_positive(&function, lower, error_bound, 1000)
                    .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations", $description),
            )?;

            Ok(())
        }
    };
    (
        name: $func_name:ident,
        function: $type:ty,
        alpha: $alphaty:ty => $alpha:expr,
        exp_result:      $exp_result:expr,
        exp_error:       $exp_error:expr,
        exp_iterations:  $exp_iterations:expr,
        exp_evaluations: $exp_evaluations:expr,
        tolerance_rel: $tolerance:expr,
        test_abs_error_bound: $abs_error_bound:expr,
        test_rel_error_bound: $rel_error_bound:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_error_bound = $abs_error_bound;
            let rel_error_bound = $rel_error_bound;
            let alpha: $alphaty = $alpha;
            let error_bound = ErrorBound::Relative($tolerance);

            let function = <$type>::new(alpha);

            let integral = AdaptiveSingularity::infinite(&function, error_bound, 1000).unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations", $description),
            )?;

            Ok(())
        }
    };
    (
        name: $func_name:ident,
        function: $type:ty,
        alpha: $alphaty:ty => $alpha:expr,
        exp_result:      $exp_result:expr,
        exp_error:       $exp_error:expr,
        exp_iterations:  $exp_iterations:expr,
        exp_evaluations: $exp_evaluations:expr,
        tolerance_abs: $tolerance:expr,
        test_abs_error_bound: $abs_error_bound:expr,
        test_rel_error_bound: $rel_error_bound:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_error_bound = $abs_error_bound;
            let rel_error_bound = $rel_error_bound;
            let alpha: $alphaty = $alpha;
            let error_bound = ErrorBound::Absolute($tolerance);

            let function = <$type>::new(alpha);

            let integral = AdaptiveSingularity::infinite(&function, error_bound, 1000).unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations", $description),
            )?;

            Ok(())
        }
    };
    (
        name: $func_name:ident,
        function: $type:ty,
        lower: $lower:expr,
        upper: $upper:expr,
        exp_result:      $exp_result:expr,
        exp_error:       $exp_error:expr,
        exp_iterations:  $exp_iterations:expr,
        exp_evaluations: $exp_evaluations:expr,
        tolerance_rel: $tolerance:expr,
        test_abs_error_bound: $abs_error_bound:expr,
        test_rel_error_bound: $rel_error_bound:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_error_bound = $abs_error_bound;
            let rel_error_bound = $rel_error_bound;
            let error_bound = ErrorBound::Relative($tolerance);

            let lower = $lower;
            let upper = $upper;

            let function = <$type>::new();

            let integral = AdaptiveSingularity::general(
                &function,
                Limits::new(lower, upper),
                error_bound,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations", $description),
            )?;

            let lower = $upper;
            let upper = $lower;
            let error_bound = ErrorBound::Relative(1e-10);

            let integral = AdaptiveSingularity::general(
                &function,
                Limits::new(lower, upper),
                error_bound,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations reverse", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations reverse", $description),
            )?;

            Ok(())
        }
    };
    (
        name: $func_name:ident,
        function: $type:ty,
        lower: $lower:expr,
        upper: $upper:expr,
        exp_result:      $exp_result:expr,
        exp_error:       $exp_error:expr,
        exp_iterations:  $exp_iterations:expr,
        exp_evaluations: $exp_evaluations:expr,
        tolerance_abs: $tolerance:expr,
        test_abs_error_bound: $abs_error_bound:expr,
        test_rel_error_bound: $rel_error_bound:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_error_bound = $abs_error_bound;
            let rel_error_bound = $rel_error_bound;
            let error_bound = ErrorBound::Absolute($tolerance);

            let lower = $lower;
            let upper = $upper;

            let function = <$type>::new();

            let integral = AdaptiveSingularity::general(
                &function,
                Limits::new(lower, upper),
                error_bound,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations", $description),
            )?;

            let lower = $upper;
            let upper = $lower;
            let error_bound = ErrorBound::Absolute(1e-10);

            let integral = AdaptiveSingularity::general(
                &function,
                Limits::new(lower, upper),
                error_bound,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations reverse", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations reverse", $description),
            )?;

            Ok(())
        }
    };
    (
        name: $func_name:ident,
        function: $type:ty,
        upper: $upper:expr,
        exp_result:      $exp_result:expr,
        exp_error:       $exp_error:expr,
        exp_iterations:  $exp_iterations:expr,
        exp_evaluations: $exp_evaluations:expr,
        tolerance_rel: $tolerance:expr,
        test_abs_error_bound: $abs_error_bound:expr,
        test_rel_error_bound: $rel_error_bound:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_error_bound = $abs_error_bound;
            let rel_error_bound = $rel_error_bound;
            let error_bound = ErrorBound::Relative($tolerance);

            let upper = $upper;

            let function = <$type>::new();

            let integral =
                AdaptiveSingularity::semi_infinite_negative(&function, upper, error_bound, 1000)
                    .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations", $description),
            )?;

            Ok(())
        }
    };
    (
        name: $func_name:ident,
        function: $type:ty,
        upper: $upper:expr,
        exp_result:      $exp_result:expr,
        exp_error:       $exp_error:expr,
        exp_iterations:  $exp_iterations:expr,
        exp_evaluations: $exp_evaluations:expr,
        tolerance_abs: $tolerance:expr,
        test_abs_error_bound: $abs_error_bound:expr,
        test_rel_error_bound: $rel_error_bound:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_error_bound = $abs_error_bound;
            let rel_error_bound = $rel_error_bound;
            let error_bound = ErrorBound::Absolute($tolerance);

            let upper = $upper;

            let function = <$type>::new();

            let integral =
                AdaptiveSingularity::semi_infinite_negative(&function, upper, error_bound, 1000)
                    .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations", $description),
            )?;

            Ok(())
        }
    };
    (
        name: $func_name:ident,
        function: $type:ty,
        lower: $lower:expr,
        exp_result:      $exp_result:expr,
        exp_error:       $exp_error:expr,
        exp_iterations:  $exp_iterations:expr,
        exp_evaluations: $exp_evaluations:expr,
        tolerance_rel: $tolerance:expr,
        test_abs_error_bound: $abs_error_bound:expr,
        test_rel_error_bound: $rel_error_bound:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_error_bound = $abs_error_bound;
            let rel_error_bound = $rel_error_bound;
            let error_bound = ErrorBound::Relative($tolerance);

            let lower = $lower;

            let function = <$type>::new();

            let integral =
                AdaptiveSingularity::semi_infinite_positive(&function, lower, error_bound, 1000)
                    .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations", $description),
            )?;

            Ok(())
        }
    };
    (
        name: $func_name:ident,
        function: $type:ty,
        lower: $lower:expr,
        exp_result:      $exp_result:expr,
        exp_error:       $exp_error:expr,
        exp_iterations:  $exp_iterations:expr,
        exp_evaluations: $exp_evaluations:expr,
        tolerance_abs: $tolerance:expr,
        test_abs_error_bound: $abs_error_bound:expr,
        test_rel_error_bound: $rel_error_bound:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_error_bound = $abs_error_bound;
            let rel_error_bound = $rel_error_bound;
            let error_bound = ErrorBound::Absolute($tolerance);

            let lower = $lower;

            let function = <$type>::new();

            let integral =
                AdaptiveSingularity::semi_infinite_positive(&function, lower, error_bound, 1000)
                    .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations", $description),
            )?;

            Ok(())
        }
    };
    (
        name: $func_name:ident,
        function: $type:ty,
        exp_result:      $exp_result:expr,
        exp_error:       $exp_error:expr,
        exp_iterations:  $exp_iterations:expr,
        exp_evaluations: $exp_evaluations:expr,
        tolerance_rel: $tolerance:expr,
        test_abs_error_bound: $abs_error_bound:expr,
        test_rel_error_bound: $rel_error_bound:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_error_bound = $abs_error_bound;
            let rel_error_bound = $rel_error_bound;
            let error_bound = ErrorBound::Relative($tolerance);

            let function = <$type>::new();

            let integral = AdaptiveSingularity::infinite(&function, error_bound, 1000).unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations", $description),
            )?;

            Ok(())
        }
    };
    (
        name: $func_name:ident,
        function: $type:ty,
        exp_result:      $exp_result:expr,
        exp_error:       $exp_error:expr,
        exp_iterations:  $exp_iterations:expr,
        exp_evaluations: $exp_evaluations:expr,
        tolerance_abs: $tolerance:expr,
        test_abs_error_bound: $abs_error_bound:expr,
        test_rel_error_bound: $rel_error_bound:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_error_bound = $abs_error_bound;
            let rel_error_bound = $rel_error_bound;
            let error_bound = ErrorBound::Absolute($tolerance);

            let function = <$type>::new();

            let integral = AdaptiveSingularity::infinite(&function, error_bound, 1000).unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.function_evaluations();

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
            util::test_int(
                iterations,
                exp_iterations,
                &format!("{} iterations", $description),
            )?;
            util::test_int(
                evaluations,
                exp_evaluations,
                &format!("{} evaluations", $description),
            )?;

            Ok(())
        }
    };
}
