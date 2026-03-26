#[macro_export]
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
        abs_tolerance: $abs_tolerance:expr,
        rel_tolerance: $rel_tolerance:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;

            let abs_tolerance = $abs_tolerance;
            let rel_tolerance = $rel_tolerance;
            let alpha: $alphaty = $alpha;

            let lower = $lower;
            let upper = $upper;

            let function = <$type>::new(alpha);
            let rule = Rule::$rule();
            let integral = Basic::new(&function, &rule, Limits::new(lower, upper).unwrap());

            let integral_result = integral.integrate();
            let result = integral_result.result();
            let error = integral_result.error();

            util::test_relative_error(
                result,
                exp_result,
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
                &format!("{} error", $description),
            )?;

            let lower = $upper;
            let upper = $lower;

            let rule = Rule::$rule();
            let integral = Basic::new(&function, &rule, Limits::new(lower, upper).unwrap());

            let integral_result = integral.integrate();
            let result = integral_result.result();
            let error = integral_result.error();

            util::test_relative_error(
                result,
                -exp_result,
                abs_tolerance,
                &format!("{} result reverse", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
                &format!("{} error reverse", $description),
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
        abs_tolerance: $abs_tolerance:expr,
        rel_tolerance: $rel_tolerance:expr
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;

            let abs_tolerance = $abs_tolerance;
            let rel_tolerance = $rel_tolerance;

            let lower = $lower;
            let upper = $upper;

            let function = <$type>::new();
            let rule = Rule::$rule();
            let integral = Basic::new(&function, &rule, Limits::new(lower, upper).unwrap());

            let integral_result = integral.integrate();
            let result = integral_result.result();
            let error = integral_result.error();

            util::test_relative_error(
                result,
                exp_result,
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
                &format!("{} error", $description),
            )?;

            let lower = $upper;
            let upper = $lower;

            let rule = Rule::$rule();
            let integral = Basic::new(&function, &rule, Limits::new(lower, upper).unwrap());

            let integral_result = integral.integrate();
            let result = integral_result.result();
            let error = integral_result.error();

            util::test_relative_error(
                result,
                -exp_result,
                abs_tolerance,
                &format!("{} result reverse", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
                &format!("{} error reverse", $description),
            )?;

            Ok(())
        }
    };
}

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
        test_abs_tolerance: $abs_tolerance:expr,
        test_rel_tolerance: $rel_tolerance:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_tolerance = $abs_tolerance;
            let rel_tolerance = $rel_tolerance;
            let rule = Rule::$rule();
            let tolerance = Tolerance::Relative($tolerance);
            let alpha: $alphaty = $alpha;

            let lower = $lower;
            let upper = $upper;

            let function = <$type>::new(alpha);

            let integral = Adaptive::new(
                &function,
                &rule,
                Limits::new(lower, upper).unwrap(),
                tolerance,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                exp_result,
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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
            let tolerance = Tolerance::Relative($tolerance);
            let lower = $upper;
            let upper = $lower;

            let integral = Adaptive::new(
                &function,
                &rule,
                Limits::new(lower, upper).unwrap(),
                tolerance,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                -exp_result,
                abs_tolerance,
                &format!("{} result reverse", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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
        test_abs_tolerance: $abs_tolerance:expr,
        test_rel_tolerance: $rel_tolerance:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_tolerance = $abs_tolerance;
            let rel_tolerance = $rel_tolerance;
            let rule = Rule::$rule();
            let tolerance = Tolerance::Absolute($tolerance);
            let alpha: $alphaty = $alpha;

            let lower = $lower;
            let upper = $upper;

            let function = <$type>::new(alpha);

            let integral = Adaptive::new(
                &function,
                &rule,
                Limits::new(lower, upper).unwrap(),
                tolerance,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                exp_result,
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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
            let tolerance = Tolerance::Absolute($tolerance);
            let lower = $upper;
            let upper = $lower;

            let integral = Adaptive::new(
                &function,
                &rule,
                Limits::new(lower, upper).unwrap(),
                tolerance,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                -exp_result,
                abs_tolerance,
                &format!("{} result reverse", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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
        test_abs_tolerance: $abs_tolerance:expr,
        test_rel_tolerance: $rel_tolerance:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_tolerance = $abs_tolerance;
            let rel_tolerance = $rel_tolerance;
            let tolerance = Tolerance::Relative($tolerance);
            let alpha: $alphaty = $alpha;

            let lower = $lower;
            let upper = $upper;

            let function = <$type>::new(alpha);

            let integral = AdaptiveSingularity::finite(
                &function,
                Limits::new(lower, upper).unwrap(),
                tolerance,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                exp_result,
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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
            let tolerance = Tolerance::Relative($tolerance);

            let integral = AdaptiveSingularity::finite(
                &function,
                Limits::new(lower, upper).unwrap(),
                tolerance,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                -exp_result,
                abs_tolerance,
                &format!("{} result reverse", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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
        test_abs_tolerance: $abs_tolerance:expr,
        test_rel_tolerance: $rel_tolerance:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_tolerance = $abs_tolerance;
            let rel_tolerance = $rel_tolerance;
            let tolerance = Tolerance::Absolute($tolerance);
            let alpha: $alphaty = $alpha;

            let lower = $lower;
            let upper = $upper;

            let function = <$type>::new(alpha);

            let integral = AdaptiveSingularity::finite(
                &function,
                Limits::new(lower, upper).unwrap(),
                tolerance,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                exp_result,
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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
            let tolerance = Tolerance::Absolute($tolerance);

            let integral = AdaptiveSingularity::finite(
                &function,
                Limits::new(lower, upper).unwrap(),
                tolerance,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                -exp_result,
                abs_tolerance,
                &format!("{} result reverse", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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
        test_abs_tolerance: $abs_tolerance:expr,
        test_rel_tolerance: $rel_tolerance:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_tolerance = $abs_tolerance;
            let rel_tolerance = $rel_tolerance;
            let alpha: $alphaty = $alpha;
            let tolerance = Tolerance::Relative($tolerance);

            let upper = $upper;

            let function = <$type>::new(alpha);

            let integral =
                AdaptiveSingularity::semi_infinite_lower(&function, upper, tolerance, 1000)
                    .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                exp_result,
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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
        test_abs_tolerance: $abs_tolerance:expr,
        test_rel_tolerance: $rel_tolerance:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_tolerance = $abs_tolerance;
            let rel_tolerance = $rel_tolerance;
            let alpha: $alphaty = $alpha;
            let tolerance = Tolerance::Absolute($tolerance);

            let upper = $upper;

            let function = <$type>::new(alpha);

            let integral =
                AdaptiveSingularity::semi_infinite_lower(&function, upper, tolerance, 1000)
                    .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                exp_result,
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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
        test_abs_tolerance: $abs_tolerance:expr,
        test_rel_tolerance: $rel_tolerance:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_tolerance = $abs_tolerance;
            let rel_tolerance = $rel_tolerance;
            let alpha: $alphaty = $alpha;
            let tolerance = Tolerance::Relative($tolerance);

            let lower = $lower;

            let function = <$type>::new(alpha);

            let integral =
                AdaptiveSingularity::semi_infinite_upper(&function, lower, tolerance, 1000)
                    .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                exp_result,
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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
        test_abs_tolerance: $abs_tolerance:expr,
        test_rel_tolerance: $rel_tolerance:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_tolerance = $abs_tolerance;
            let rel_tolerance = $rel_tolerance;
            let alpha: $alphaty = $alpha;
            let tolerance = Tolerance::Absolute($tolerance);

            let lower = $lower;

            let function = <$type>::new(alpha);

            let integral =
                AdaptiveSingularity::semi_infinite_upper(&function, lower, tolerance, 1000)
                    .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                exp_result,
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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
        test_abs_tolerance: $abs_tolerance:expr,
        test_rel_tolerance: $rel_tolerance:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_tolerance = $abs_tolerance;
            let rel_tolerance = $rel_tolerance;
            let alpha: $alphaty = $alpha;
            let tolerance = Tolerance::Relative($tolerance);

            let function = <$type>::new(alpha);

            let integral = AdaptiveSingularity::infinite(&function, tolerance, 1000).unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                exp_result,
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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
        test_abs_tolerance: $abs_tolerance:expr,
        test_rel_tolerance: $rel_tolerance:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_tolerance = $abs_tolerance;
            let rel_tolerance = $rel_tolerance;
            let alpha: $alphaty = $alpha;
            let tolerance = Tolerance::Absolute($tolerance);

            let function = <$type>::new(alpha);

            let integral = AdaptiveSingularity::infinite(&function, tolerance, 1000).unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                exp_result,
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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
        test_abs_tolerance: $abs_tolerance:expr,
        test_rel_tolerance: $rel_tolerance:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_tolerance = $abs_tolerance;
            let rel_tolerance = $rel_tolerance;
            let tolerance = Tolerance::Relative($tolerance);

            let lower = $lower;
            let upper = $upper;

            let function = <$type>::new();

            let integral = AdaptiveSingularity::finite(
                &function,
                Limits::new(lower, upper).unwrap(),
                tolerance,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                exp_result,
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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
            let tolerance = Tolerance::Relative(1e-10);

            let integral = AdaptiveSingularity::finite(
                &function,
                Limits::new(lower, upper).unwrap(),
                tolerance,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                -exp_result,
                abs_tolerance,
                &format!("{} result reverse", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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
        test_abs_tolerance: $abs_tolerance:expr,
        test_rel_tolerance: $rel_tolerance:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_tolerance = $abs_tolerance;
            let rel_tolerance = $rel_tolerance;
            let tolerance = Tolerance::Absolute($tolerance);

            let lower = $lower;
            let upper = $upper;

            let function = <$type>::new();

            let integral = AdaptiveSingularity::finite(
                &function,
                Limits::new(lower, upper).unwrap(),
                tolerance,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                exp_result,
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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
            let tolerance = Tolerance::Absolute(1e-10);

            let integral = AdaptiveSingularity::finite(
                &function,
                Limits::new(lower, upper).unwrap(),
                tolerance,
                1000,
            )
            .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                -exp_result,
                abs_tolerance,
                &format!("{} result reverse", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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
        test_abs_tolerance: $abs_tolerance:expr,
        test_rel_tolerance: $rel_tolerance:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_tolerance = $abs_tolerance;
            let rel_tolerance = $rel_tolerance;
            let tolerance = Tolerance::Relative($tolerance);

            let upper = $upper;

            let function = <$type>::new();

            let integral =
                AdaptiveSingularity::semi_infinite_lower(&function, upper, tolerance, 1000)
                    .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                exp_result,
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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
        test_abs_tolerance: $abs_tolerance:expr,
        test_rel_tolerance: $rel_tolerance:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_tolerance = $abs_tolerance;
            let rel_tolerance = $rel_tolerance;
            let tolerance = Tolerance::Absolute($tolerance);

            let upper = $upper;

            let function = <$type>::new();

            let integral =
                AdaptiveSingularity::semi_infinite_lower(&function, upper, tolerance, 1000)
                    .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                exp_result,
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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
        test_abs_tolerance: $abs_tolerance:expr,
        test_rel_tolerance: $rel_tolerance:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_tolerance = $abs_tolerance;
            let rel_tolerance = $rel_tolerance;
            let tolerance = Tolerance::Relative($tolerance);

            let lower = $lower;

            let function = <$type>::new();

            let integral =
                AdaptiveSingularity::semi_infinite_upper(&function, lower, tolerance, 1000)
                    .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                exp_result,
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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
        test_abs_tolerance: $abs_tolerance:expr,
        test_rel_tolerance: $rel_tolerance:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_tolerance = $abs_tolerance;
            let rel_tolerance = $rel_tolerance;
            let tolerance = Tolerance::Absolute($tolerance);

            let lower = $lower;

            let function = <$type>::new();

            let integral =
                AdaptiveSingularity::semi_infinite_upper(&function, lower, tolerance, 1000)
                    .unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                exp_result,
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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
        test_abs_tolerance: $abs_tolerance:expr,
        test_rel_tolerance: $rel_tolerance:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_tolerance = $abs_tolerance;
            let rel_tolerance = $rel_tolerance;
            let tolerance = Tolerance::Relative($tolerance);

            let function = <$type>::new();

            let integral = AdaptiveSingularity::infinite(&function, tolerance, 1000).unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                exp_result,
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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
        test_abs_tolerance: $abs_tolerance:expr,
        test_rel_tolerance: $rel_tolerance:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result = $exp_result;
            let exp_error = $exp_error;
            let exp_iterations = $exp_iterations;
            let exp_evaluations = $exp_evaluations;

            let abs_tolerance = $abs_tolerance;
            let rel_tolerance = $rel_tolerance;
            let tolerance = Tolerance::Absolute($tolerance);

            let function = <$type>::new();

            let integral = AdaptiveSingularity::infinite(&function, tolerance, 1000).unwrap();

            let integral_result = integral.integrate().unwrap();
            let result = integral_result.result();
            let error = integral_result.error();
            let iterations = integral_result.iterations();
            let evaluations = integral_result.evaluations();

            util::test_relative_error(
                result,
                exp_result,
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
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

#[macro_export]
macro_rules! basic_test_complex_equal {
    (
        name: $func_name:ident,
        function: $type:ty,
        alpha: $alphaty:ty => $alpha:expr,
        rule: $rule:ident,
        lower: $lower:expr,
        upper: $upper:expr,
        exp_result: $exp_result:expr,
        exp_error: $exp_error:expr,
        abs_tolerance: $abs_tolerance:expr,
        rel_tolerance: $rel_tolerance:expr,
        $description:literal
    ) => {
        #[test]
        fn $func_name() -> Result<(), String> {
            let exp_result_re = $exp_result;
            let exp_result_im = $exp_result;
            let exp_result = Complex::new(exp_result_re, exp_result_im);
            let exp_error_re = $exp_error;
            let exp_error_im = $exp_error;
            let exp_error_complex = Complex::new(exp_error_re, exp_error_im);
            let exp_error = exp_error_complex.abs();

            let abs_tolerance = $abs_tolerance;
            let rel_tolerance = $rel_tolerance;
            let alpha1: $alphaty = $alpha;
            let alpha2: $alphaty = $alpha;

            let lower = $lower;
            let upper = $upper;

            let function = <$type>::new(alpha1, alpha2);
            let rule = Rule::$rule();
            let integral = Basic::new(&function, &rule, Limits::new(lower, upper).unwrap());

            let integral_result = integral.integrate();
            let result = integral_result.result();
            let result_re = result.re();
            let result_im = result.im();
            let error = integral_result.error();

            util::test_relative_error(
                result_re,
                exp_result.re(),
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                result_im,
                exp_result.im(),
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
                &format!("{} error", $description),
            )?;

            let lower = $upper;
            let upper = $lower;

            let rule = Rule::$rule();
            let integral = Basic::new(&function, &rule, Limits::new(lower, upper).unwrap());

            let integral_result = integral.integrate();
            let result = integral_result.result();
            let result_re = result.re();
            let result_im = result.im();
            let error = integral_result.error();

            util::test_relative_error(
                result_re,
                -exp_result_re,
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                result_im,
                -exp_result_im,
                abs_tolerance,
                &format!("{} result", $description),
            )?;
            util::test_relative_error(
                error,
                exp_error,
                rel_tolerance,
                &format!("{} error", $description),
            )?;

            Ok(())
        }
    };
}
