use num_complex::ComplexFloat;
use num_traits::Zero;

use crate::multi::geometry::Geometry;
use crate::multi::region::Region;
use crate::multi::rule::Rule;
use crate::Limits;
use crate::MultiDimensionalIntegrand;
use crate::ScalarF64;

pub(crate) struct Integrator<'a, I, const NDIM: usize, const FINAL: usize, const TOTAL: usize>
where
    I: MultiDimensionalIntegrand<NDIM>,
{
    function: &'a I,
    rule: &'a Rule<NDIM, FINAL, TOTAL>,
    limits: [Limits; NDIM],
}

impl<'a, I, const NDIM: usize, const FINAL: usize, const TOTAL: usize>
    Integrator<'a, I, NDIM, FINAL, TOTAL>
where
    I: MultiDimensionalIntegrand<NDIM>,
{
    pub(crate) const fn new(
        function: &'a I,
        rule: &'a Rule<NDIM, FINAL, TOTAL>,
        limits: [Limits; NDIM],
    ) -> Self {
        Self {
            function,
            rule,
            limits,
        }
    }

    const fn geometry(&self) -> Geometry<NDIM> {
        Geometry::new(&self.limits)
    }

    pub(crate) fn integrate(&self) -> Region<NDIM, I::Scalar> {
        let evaluations = self.rule.evaluations();
        let basic_error_coeff = self.rule.basic_error_coeff();
        let ratio = self.rule.ratio();

        let Geometry {
            centre,
            half_widths,
            volume,
            largest_axis,
        } = self.geometry();

        let mut bisection_index = largest_axis;

        let [data0, data1, data2] = self.rule.initial_data();

        let point0 = data0.generator().point(&centre, &half_widths);
        let fval_0 = self.function.evaluate(&point0);

        // NOTE see note below
        //let fval_0_abs = fval_0.abs();

        let mut intermediate_result = IntermediateResult::new(
            fval_0 * data0.weight(),
            fval_0 * data0.null1(),
            fval_0 * data0.null2(),
            fval_0 * data0.null3(),
            fval_0 * data0.null4(),
        );

        let gen1 = data1.get_first_value();
        let gen2 = data2.get_first_value();

        let mut point = centre;

        let mut difmax = 0.0;

        for i in 0..NDIM {
            point[i] = centre[i] - half_widths[i] * gen1;
            let fval_1_minus = self.function.evaluate(&point);

            point[i] = centre[i] + half_widths[i] * gen1;
            let fval_1_plus = self.function.evaluate(&point);

            let fval_1 = fval_1_minus + fval_1_plus;

            point[i] = centre[i] - half_widths[i] * gen2;
            let fval_2_minus = self.function.evaluate(&point);

            point[i] = centre[i] + half_widths[i] * gen2;
            let fval_2_plus = self.function.evaluate(&point);

            let fval_2 = fval_2_minus + fval_2_plus;

            point[i] = centre[i];

            let fourth_diff = fval_0 * 2.0 * (1.0 - ratio) - fval_2 + fval_1 * ratio;
            let fourth_diff_abs = fourth_diff.abs();

            // NOTE: dcuhre uses this, we have removed it for now, as comparison
            // tests with dcuhre currently pass
            // Ignore differences below roundoff
            //let difsum = if fval_0_abs + fourth_diff_abs / 4.0 != fval_0_abs {
            //    fourth_diff_abs
            //} else {
            //    0.0
            //};

            let difsum = fourth_diff_abs;

            intermediate_result.add(
                fval_1 * data1.weight() + fval_2 * data2.weight(),
                fval_1 * data1.null1() + fval_2 * data2.null1(),
                fval_1 * data1.null2() + fval_2 * data2.null2(),
                fval_1 * data1.null3() + fval_2 * data2.null3(),
                fval_1 * data1.null4() + fval_2 * data2.null4(),
            );

            if difsum > difmax {
                difmax = difsum;
                bisection_index = i;
            }
        }

        let final_data = self.rule.final_data();

        for data in final_data {
            let zero = <I::Scalar as Zero>::zero();
            let fval = data
                .generator()
                .point_permutations(&centre, &half_widths)
                .map(|p| self.function.evaluate(&p))
                .fold(zero, |a, v| a + v);

            intermediate_result.add(
                fval * data.weight(),
                fval * data.null1(),
                fval * data.null2(),
                fval * data.null3(),
                fval * data.null4(),
            );
        }

        let mut nulls = [
            intermediate_result.null1,
            intermediate_result.null2,
            intermediate_result.null3,
            intermediate_result.null4,
        ];

        let mut nulls_abs = [
            intermediate_result.null1.abs(),
            intermediate_result.null2.abs(),
            intermediate_result.null3.abs(),
            intermediate_result.null4.abs(),
        ];

        let scales = self.rule.scales().0;
        let norms = self.rule.norms().0;

        for i in 0..3 {
            let mut search = <I::Scalar as Zero>::zero();
            for k in 0..(Rule::<NDIM, FINAL, TOTAL>::total()) {
                let other = (nulls[i + 1] + nulls[i] * scales[i][k]) * norms[i][k];
                search = if search.abs() < other.abs() {
                    other
                } else {
                    search
                };
            }
            nulls[i] = search;
            nulls_abs[i] = search.abs();
        }

        let [n1, n2, n3, _] = nulls;
        let [n1, n2, n3] = [n1.abs(), n2.abs(), n3.abs()];

        let c1 = basic_error_coeff.c1();
        let c2 = basic_error_coeff.c2();
        let c3 = basic_error_coeff.c3();
        let c4 = basic_error_coeff.c4();

        let error = if c1 * n1 <= n2 && c2 * n2 <= n3 {
            c3 * n1 * volume
        } else {
            c4 * n1.max(n2.max(n3)) * volume
        };
        let error = error.abs();
        let result = intermediate_result.result * volume;

        Region::unevaluated()
            .with_error(error)
            .with_result(result)
            .with_limits(self.limits)
            .with_function_evaluations(evaluations)
            .with_largest_error_axis(bisection_index)
    }
}

struct IntermediateResult<T> {
    result: T,
    null1: T,
    null2: T,
    null3: T,
    null4: T,
}

impl<T: ScalarF64> IntermediateResult<T> {
    fn new(result: T, null1: T, null2: T, null3: T, null4: T) -> Self {
        Self {
            result,
            null1,
            null2,
            null3,
            null4,
        }
    }

    fn add(&mut self, result: T, null1: T, null2: T, null3: T, null4: T) -> &mut Self {
        self.result += result;
        self.null1 += null1;
        self.null2 += null2;
        self.null3 += null3;
        self.null4 += null4;
        self
    }
}
