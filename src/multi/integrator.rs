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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::multi::{Rule07, Rule09, Rule11, Rule13};

    #[test]
    fn compare_basic_7point_with_dcuhre_output_ndim_2() {
        const NDIM: usize = 2;

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let exponent = -x * y;
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.79659970249839818;
            let dcuhre_error = 2.7359932817440052E-005;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-13);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-7);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let exponent = -(x * x + y * y);
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.55774638300971069;
            let dcuhre_error = 1.3049001762496107E-005;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-12);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-7);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let exponent = -((x * x).cos().powi(2) * (y * y).cos().powi(2));
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.51478351071879791;
            let dcuhre_error = 0.37449223525594855;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-10);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-9);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    f64::cos((1.0 + x * x).ln() * (1.0 + y * y).ln())
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.99332124356102158;
            let dcuhre_error = 9.3939393118662768E-005;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-10);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-7);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    (x * x / (2.0 - x.cos())) + (y * y / (2.0 - y.cos()))
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limits = [Limits::new(-2.0, 3.0), Limits::new(1.0, 10.0)];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 887.97458362328393;
            let dcuhre_error = 611.03739938210026;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-7);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-8);
        }
    }

    #[test]
    fn compare_basic_7point_with_dcuhre_output_ndim_3() {
        const NDIM: usize = 3;

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent = -x * y * z;
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.89120896753764178;
            let dcuhre_error = 1.3869851434747299E-005;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-12);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-7);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent = -(x * x + y * y + z * z);
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.41653847897138208;
            let dcuhre_error = 5.8625669977822451E-006;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-12);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 2e-7);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent =
                        -((x * x).cos().powi(2) * (y * y).cos().powi(2) * (z * z).cos().powi(2));
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.57845060677495597;
            let dcuhre_error = 0.47067683202127575;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-11);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-8);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    f64::cos((1.0 + x * x).ln() * (1.0 + y * y).ln() * (1.0 + z * z).ln())
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.99921595669136032;
            let dcuhre_error = 2.5792992219862503E-003;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-12);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-8);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    (x * x / (2.0 - x.cos()))
                        + (y * y / (2.0 - y.cos()))
                        + (z * z / (2.0 - z.cos()))
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limits = [
                Limits::new(-2.0, 3.0),
                Limits::new(1.0, 10.0),
                Limits::new(0.0, -1.0),
            ];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = -899.77650881070997;
            let dcuhre_error = 543.83282190394414;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-7);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-8);
        }
    }

    #[test]
    fn compare_basic_9point_with_dcuhre_output_ndim_3() {
        const NDIM: usize = 3;

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent = -x * y * z;
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule09::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.89121284475576390;
            let dcuhre_error = 2.8431093456305666E-003;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-13);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent = -(x * x + y * y + z * z);
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule09::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.41653839142919569;
            let dcuhre_error = 7.2117833072453187E-008;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-8);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent =
                        -((x * x).cos().powi(2) * (y * y).cos().powi(2) * (z * z).cos().powi(2));
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule09::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.57854494486523600;
            let dcuhre_error = 4.2597555388569463E-005;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-14);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-10);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    f64::cos((1.0 + x * x).ln() * (1.0 + y * y).ln() * (1.0 + z * z).ln())
                }
            }

            let function = Function;
            let rule = Rule09::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.99922153302627292;
            let dcuhre_error = 7.5347754069994284E-004;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-12);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    (x * x / (2.0 - x.cos()))
                        + (y * y / (2.0 - y.cos()))
                        + (z * z / (2.0 - z.cos()))
                }
            }

            let function = Function;
            let rule = Rule09::generate().unwrap();
            let limits = [
                Limits::new(-2.0, 3.0),
                Limits::new(1.0, 10.0),
                Limits::new(0.0, -1.0),
            ];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = -905.91915061432690;
            let dcuhre_error = 144.00666640485426;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-14);
        }
    }

    #[test]
    fn compare_basic_11point_with_dcuhre_output_ndim_3() {
        const NDIM: usize = 3;

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent = -x * y * z;
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule11::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.89121279793446351;
            let dcuhre_error = 1.5894785143767116E-008;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-20);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-20);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent = -(x * x + y * y + z * z);
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule11::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.41653838596655557;
            let dcuhre_error = 3.2348830851556911E-006;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-20);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-20);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent =
                        -((x * x).cos().powi(2) * (y * y).cos().powi(2) * (z * z).cos().powi(2));
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule11::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.57853848358343507;
            let dcuhre_error = 2.3874271216029219E-003;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-20);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-13);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    f64::cos((1.0 + x * x).ln() * (1.0 + y * y).ln() * (1.0 + z * z).ln())
                }
            }

            let function = Function;
            let rule = Rule11::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.99922120779864210;
            let dcuhre_error = 1.3732744589107246E-007;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-20);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-9);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    (x * x / (2.0 - x.cos()))
                        + (y * y / (2.0 - y.cos()))
                        + (z * z / (2.0 - z.cos()))
                }
            }

            let function = Function;
            let rule = Rule11::generate();
            let limits = [
                Limits::new(-2.0, 3.0),
                Limits::new(1.0, 10.0),
                Limits::new(0.0, -1.0),
            ];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = -904.73349759612938;
            let dcuhre_error = 36.531147257907918;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-20);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-14);
        }
    }

    #[test]
    fn compare_basic_13point_with_dcuhre_output_ndim_2() {
        const NDIM: usize = 2;

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let exponent = -x * y;
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule13::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.79659959929708202;
            let dcuhre_error = 3.2673674642535072E-013;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-20);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-20);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let exponent = -(x * x + y * y);
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule13::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.55774628535131576;
            let dcuhre_error = 9.7456124932917485E-012;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-20);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-20);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let exponent = -((x * x).cos().powi(2) * (y * y).cos().powi(2));
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule13::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.51478983417297963;
            let dcuhre_error = 4.4473459718302440E-004;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-20);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-12);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    f64::cos((1.0 + x * x).ln() * (1.0 + y * y).ln())
                }
            }

            let function = Function;
            let rule = Rule13::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.99331725120688252;
            let dcuhre_error = 8.0707600531214359E-010;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-20);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-20);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    (x * x / (2.0 - x.cos())) + (y * y / (2.0 - y.cos()))
                }
            }

            let function = Function;
            let rule = Rule13::generate();
            let limits = [Limits::new(-2.0, 3.0), Limits::new(1.0, 10.0)];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 911.85973569354837;
            let dcuhre_error = 129.79778763945998;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-20);
        }
    }
}
