use num_complex::ComplexFloat;
use num_traits::Zero;

use crate::multi::geometry::Geometry;
use crate::multi::region::Region;
use crate::multi::rule::Data;
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

    fn initial_integration(
        &self,
        centre: &[f64; NDIM],
        half_widths: &[f64; NDIM],
        largest_axis: usize,
    ) -> IntermediateResult<I::Scalar> {
        let ratio = self.rule.ratio();
        let [data0, data1, data2] = self.rule.initial_data();

        let point0 = data0.generator().point(centre, half_widths);
        let fval_0 = self.function.evaluate(&point0);

        let mut intermediate_result = IntermediateResult::new(fval_0, data0, largest_axis);

        let gen1 = data1.get_first_value();
        let gen2 = data2.get_first_value();

        let mut point = *centre;

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

            let difsum = if (fval_0 + fourth_diff / 4.0).abs().to_bits() != fval_0.abs().to_bits() {
                fourth_diff_abs
            } else {
                0.0
            };

            intermediate_result.add(fval_1, data1).add(fval_2, data2);

            if difsum > difmax {
                difmax = difsum;
                intermediate_result.update_bisection_axis(i);
            }
        }

        intermediate_result
    }

    fn complete_integration(
        &self,
        centre: &[f64; NDIM],
        half_widths: &[f64; NDIM],
        volume: f64,
        initial: IntermediateResult<I::Scalar>,
    ) -> (I::Scalar, f64) {
        let mut intermediate_result = initial;
        let final_data = self.rule.final_data();

        for data in final_data {
            let zero = I::Scalar::zero();
            let fval = data
                .generator()
                .point_permutations(centre, half_widths)
                .map(|p| self.function.evaluate(&p))
                .fold(zero, |a, v| a + v);

            intermediate_result.add(fval, data);
        }

        let result = intermediate_result.result * volume;
        let error = self.calculate_error(&intermediate_result, volume);

        (result, error)
    }

    fn calculate_error(
        &self,
        intermediate_result: &IntermediateResult<I::Scalar>,
        volume: f64,
    ) -> f64 {
        let consecutive_nulls = [
            (intermediate_result.null1, intermediate_result.null2),
            (intermediate_result.null2, intermediate_result.null3),
            (intermediate_result.null3, intermediate_result.null4),
        ];

        let scales_norms = self.rule.scales_norms();

        let mut max_nulls_abs = [0f64; 3];

        for (i, ((n1, n2), sn)) in consecutive_nulls
            .iter()
            .zip(scales_norms.iter())
            .enumerate()
        {
            let max = sn.max_consecutive_null(*n1, *n2).abs();
            max_nulls_abs[i] = max;
        }

        let [n1, n2, n3] = max_nulls_abs;

        let basic_error_coeff = self.rule.basic_error_coeff();
        let c1 = basic_error_coeff.c1();
        let c2 = basic_error_coeff.c2();
        let c3 = basic_error_coeff.c3();
        let c4 = basic_error_coeff.c4();

        let error = if c1 * n1 <= n2 && c2 * n2 <= n3 {
            c3 * n1 * volume
        } else {
            c4 * n1.max(n2.max(n3)) * volume
        };
        error.abs()
    }

    pub(crate) fn integrate(&self) -> Region<I::Scalar, NDIM> {
        let Geometry {
            centre,
            half_widths,
            volume,
            largest_axis,
        } = self.geometry();

        let initial_result = self.initial_integration(&centre, &half_widths, largest_axis);

        let bisection_axis = initial_result.bisection_axis();

        let (result, error) =
            self.complete_integration(&centre, &half_widths, volume, initial_result);

        Region::unevaluated()
            .with_result(result)
            .with_error(error)
            .with_bisection_axis(bisection_axis)
            .with_limits(self.limits)
            .with_volume(volume)
    }
}

struct IntermediateResult<T> {
    result: T,
    null1: T,
    null2: T,
    null3: T,
    null4: T,
    bisection_axis: usize,
}

impl<T: ScalarF64> IntermediateResult<T> {
    fn new<const NDIM: usize>(fval: T, data: &Data<NDIM>, bisection_axis: usize) -> Self {
        Self {
            result: fval * data.weight(),
            null1: fval * data.null1(),
            null2: fval * data.null2(),
            null3: fval * data.null3(),
            null4: fval * data.null4(),
            bisection_axis,
        }
    }

    fn add<const NDIM: usize>(&mut self, fval: T, data: &Data<NDIM>) -> &mut Self {
        self.result += fval * data.weight();
        self.null1 += fval * data.null1();
        self.null2 += fval * data.null2();
        self.null3 += fval * data.null3();
        self.null4 += fval * data.null4();
        self
    }

    fn update_bisection_axis(&mut self, bisection_axis: usize) -> &mut Self {
        self.bisection_axis = bisection_axis;
        self
    }

    fn bisection_axis(&self) -> usize {
        self.bisection_axis
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::multi::{Rule07, Rule09, Rule09N2, Rule11, Rule13};

    #[test]
    #[allow(clippy::too_many_lines)]
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 0.796_599_702_498_398_18;
            let dcuhre_error = 2.735_993_281_744_005_2E-005;
            let dcuhre_bisection_axis = 0;

            assert_eq!(axis, dcuhre_bisection_axis);
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 0.557_746_383_009_710_69;
            let dcuhre_error = 1.304_900_176_249_610_7E-005;
            let dcuhre_bisection_axis = 0;

            assert_eq!(axis, dcuhre_bisection_axis);
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 0.514_783_510_718_797_91;
            let dcuhre_error = 0.374_492_235_255_948_55;
            let dcuhre_bisection_axis = 0;

            assert_eq!(axis, dcuhre_bisection_axis);
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 0.993_321_243_561_021_58;
            let dcuhre_error = 9.393_939_311_866_276_8E-005;
            let dcuhre_bisection_axis = 0;

            assert_eq!(axis, dcuhre_bisection_axis);
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 887.974_583_623_283_93;
            let dcuhre_error = 611.037_399_382_100_26;
            let dcuhre_bisection_axis = 1;

            assert_eq!(axis, dcuhre_bisection_axis);
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-7);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-8);
        }
    }

    #[test]
    #[allow(clippy::too_many_lines)]
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 0.891_208_967_537_641_78;
            let dcuhre_error = 1.386_985_143_474_729_9E-005;
            let dcuhre_bisection_axis = 0;

            assert_eq!(axis, dcuhre_bisection_axis);
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 0.578_450_606_774_955_97;
            let dcuhre_error = 0.470_676_832_021_275_75;
            let dcuhre_bisection_axis = 0;

            assert_eq!(axis, dcuhre_bisection_axis);
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 0.999_215_956_691_360_32;
            let dcuhre_error = 2.579_299_221_986_250_3E-003;
            let dcuhre_bisection_axis = 0;

            assert_eq!(axis, dcuhre_bisection_axis);
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = -899.776_508_810_709_97;
            let dcuhre_error = 543.832_821_903_944_14;
            let dcuhre_bisection_axis = 1;

            assert_eq!(axis, dcuhre_bisection_axis);
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-7);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-8);
        }
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn compare_basic_9point_with_dcuhre_output_ndim_2() {
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
            let rule = Rule09N2::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 0.796_599_604_742_765_92;
            let dcuhre_error = 1.179_045_419_797_914_0E-007;
            let dcuhre_bisection_axis = 0;

            assert_eq!(axis, dcuhre_bisection_axis);
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-9);
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
            let rule = Rule09N2::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 0.557_746_293_208_851_61;
            let dcuhre_error = 1.099_467_361_768_036_7E-007;
            let dcuhre_bisection_axis = 0;

            assert_eq!(axis, dcuhre_bisection_axis);
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
                    let exponent = -((x * x).cos().powi(2) * (y * y).cos().powi(2));
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule09N2::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 0.514_797_676_048_813_88;
            let dcuhre_error = 1.718_363_434_173_172_1E-005;
            let dcuhre_bisection_axis = 0;

            assert_eq!(axis, dcuhre_bisection_axis);
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
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
            let rule = Rule09N2::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 0.993_317_491_794_954_10;
            let dcuhre_error = 4.525_423_632_376_451_2E-006;
            let dcuhre_bisection_axis = 0;

            assert_eq!(axis, dcuhre_bisection_axis);
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-9);
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
            let rule = Rule09N2::generate();
            let limits = [Limits::new(-2.0, 3.0), Limits::new(1.0, 10.0)];

            let integral = Integrator::new(&function, &rule, limits);

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 894.117_253_674_850_64;
            let dcuhre_error = 205.107_021_121_106_53;
            let dcuhre_bisection_axis = 1;

            assert_eq!(axis, dcuhre_bisection_axis);
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-14);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-9);
        }
    }

    #[test]
    #[allow(clippy::too_many_lines)]
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 0.891_212_844_755_763_90;
            let dcuhre_error = 2.843_109_345_630_566_6E-003;
            let dcuhre_bisection_axis = 0;

            assert_eq!(axis, dcuhre_bisection_axis);
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 0.416_538_391_429_195_69;
            let dcuhre_error = 7.211_783_307_245_318_7E-008;
            let dcuhre_bisection_axis = 0;

            assert_eq!(axis, dcuhre_bisection_axis);
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 0.578_544_944_865_236_00;
            let dcuhre_error = 4.259_755_538_856_946_3E-005;
            let dcuhre_bisection_axis = 0;

            assert_eq!(axis, dcuhre_bisection_axis);
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 0.999_221_533_026_272_92;
            let dcuhre_error = 7.534_775_406_999_428_4E-004;
            let dcuhre_bisection_axis = 0;

            assert_eq!(axis, dcuhre_bisection_axis);
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = -905.919_150_614_326_90;
            let dcuhre_error = 144.006_666_404_854_26;
            let dcuhre_bisection_axis = 1;

            assert_eq!(axis, dcuhre_bisection_axis);
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-14);
        }
    }

    #[test]
    #[allow(clippy::too_many_lines)]
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 0.891_212_797_934_463_51;
            let dcuhre_error = 1.589_478_514_376_711_6E-008;
            let dcuhre_bisection_axis = 0;

            assert_eq!(axis, dcuhre_bisection_axis);
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 0.416_538_385_966_555_57;
            let dcuhre_error = 3.234_883_085_155_691_1E-006;
            let dcuhre_bisection_axis = 0;

            assert_eq!(axis, dcuhre_bisection_axis);
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 0.578_538_483_583_435_07;
            let dcuhre_error = 2.387_427_121_602_921_9E-003;
            let dcuhre_bisection_axis = 0;

            assert_eq!(axis, dcuhre_bisection_axis);
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 0.999_221_207_798_642_10;
            let dcuhre_error = 1.373_274_458_910_724_6E-007;
            let dcuhre_bisection_axis = 0;

            assert_eq!(axis, dcuhre_bisection_axis);
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = -904.733_497_596_129_38;
            let dcuhre_error = 36.531_147_257_907_918;
            let dcuhre_bisection_axis = 1;

            assert_eq!(axis, dcuhre_bisection_axis);
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-20);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-14);
        }
    }

    #[test]
    #[allow(clippy::too_many_lines)]
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 0.796_599_599_297_082_02;
            let dcuhre_error = 3.267_367_464_253_507_2E-013;
            let dcuhre_bisection_axis = 0;

            assert_eq!(axis, dcuhre_bisection_axis);
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 0.557_746_285_351_315_76;
            let dcuhre_error = 9.745_612_493_291_748_5E-012;
            let dcuhre_bisection_axis = 0;

            assert_eq!(axis, dcuhre_bisection_axis);
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 0.514_789_834_172_979_63;
            let dcuhre_error = 4.447_345_971_830_244_0E-004;
            let dcuhre_bisection_axis = 0;

            assert_eq!(axis, dcuhre_bisection_axis);
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 0.993_317_251_206_882_52;
            let dcuhre_error = 8.070_760_053_121_435_9E-010;
            let dcuhre_bisection_axis = 0;

            assert_eq!(axis, dcuhre_bisection_axis);
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
            let axis = integral_result.bisection_axis();
            let dcuhre_result = 911.859_735_693_548_37;
            let dcuhre_error = 129.797_787_639_459_98;
            let dcuhre_bisection_axis = 1;

            assert_eq!(axis, dcuhre_bisection_axis);
            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-20);
        }
    }
}
