#![allow(unused)]
use crate::multi::{Error, Kind};
use crate::Limits;
use crate::MultiDimensionalIntegrand;
use crate::ScalarF64;
use num_complex::ComplexFloat;
use num_traits::Zero;

pub(crate) struct MultiDimensionalRegion<const NDIM: usize, T: ScalarF64> {
    pub(crate) error: f64,
    pub(crate) result: T,
    pub(crate) result_abs: f64,
    pub(crate) error_abs: f64,
    pub(crate) limits: [Limits; NDIM],
    pub(crate) largest_error_axis: usize,
    pub(crate) function_evaluations: usize,
}

impl<const NDIM: usize, T: ScalarF64> MultiDimensionalRegion<NDIM, T> {
    fn unevaluated() -> Self {
        let zero = <T as Zero>::zero();
        Self {
            error: 0.0,
            result: zero,
            result_abs: 0.0,
            error_abs: 0.0,
            limits: [Limits::new(0.0, 0.0); NDIM],
            largest_error_axis: 0,
            function_evaluations: 0,
        }
    }

    fn with_error(mut self, error: f64) -> Self {
        self.error = error;
        self
    }

    fn with_result(mut self, result: T) -> Self {
        self.result = result;
        self
    }

    fn with_result_abs(mut self, abs: f64) -> Self {
        self.result_abs = abs;
        self
    }

    fn with_error_abs(mut self, error_abs: f64) -> Self {
        self.error_abs = error_abs;
        self
    }

    fn with_limits(mut self, limits: [Limits; NDIM]) -> Self {
        self.limits = limits;
        self
    }

    fn with_largest_error_axis(mut self, largest_error_axis: usize) -> Self {
        self.largest_error_axis = largest_error_axis;
        self
    }

    fn with_function_evaluations(mut self, function_evaluations: usize) -> Self {
        self.function_evaluations = function_evaluations;
        self
    }

    fn error(&self) -> f64 {
        self.error
    }

    fn result(&self) -> T {
        self.result
    }

    fn result_abs(&self) -> f64 {
        self.result_abs
    }

    fn error_abs(&self) -> f64 {
        self.error_abs
    }

    fn limits(&self) -> &[Limits; NDIM] {
        &self.limits
    }

    fn largest_error_axis(&self) -> usize {
        self.largest_error_axis
    }

    fn function_evaluations(&self) -> usize {
        self.function_evaluations
    }

    fn bisect<I: MultiDimensionalIntegrand<NDIM>>(
        &self,
        function: &I,
    ) -> [MultiDimensionalRegion<NDIM, I::Scalar>; 2] {
        let axis_to_bisect = self.largest_error_axis;
        let previous_limits = self.limits();

        let [lower, upper] = previous_limits[axis_to_bisect].bisect();

        let mut lower_limits = *previous_limits;
        lower_limits[axis_to_bisect] = lower;

        let mut upper_limits = *previous_limits;
        upper_limits[axis_to_bisect] = upper;

        let lower_integral =
            MultiDimensionalBasic::new_unchecked(lower_limits, function).integrate_internal();
        let upper_integral =
            MultiDimensionalBasic::new_unchecked(upper_limits, function).integrate_internal();

        [lower_integral, upper_integral]
    }
}

pub struct MultiDimensionalBasic<I, const NDIM: usize>
where
    I: MultiDimensionalIntegrand<NDIM>,
{
    limits: [Limits; NDIM],
    function: I,
}

impl<I, const NDIM: usize> MultiDimensionalBasic<I, NDIM>
where
    I: MultiDimensionalIntegrand<NDIM>,
{
    /// # Errors
    /// TODO
    #[allow(clippy::cast_possible_truncation)]
    pub fn new(limits: [Limits; NDIM], function: I) -> Result<Self, Error<I::Scalar>> {
        if NDIM < 2 || NDIM > 15 {
            Err(Error::unevaluated(Kind::WrongDimensionality))
        } else {
            Ok(Self { limits, function })
        }
    }

    fn new_unchecked(limits: [Limits; NDIM], function: I) -> Self {
        Self { limits, function }
    }

    fn geometry(&self) -> Geometry<NDIM> {
        let mut centres = [0.0; NDIM];
        let mut widths = [0.0; NDIM];

        let mut volume = 1.0;

        for j in 0..NDIM {
            centres[j] = self.limits[j].centre();
            widths[j] = self.limits[j].half_width();
            volume *= 2.0 * widths[j];
        }

        Geometry {
            centres,
            widths,
            volume,
        }
    }

    const fn minimum_function_calls() -> usize {
        minimum_function_calls(NDIM)
    }

    #[allow(clippy::too_many_lines)]
    fn integrate_internal(&self) -> MultiDimensionalRegion<NDIM, I::Scalar> {
        let Geometry {
            centres,
            widths,
            volume,
        } = self.geometry();

        let sum1 = self.function.evaluate(&centres);
        let sum1_abs = sum1.abs();

        let zero = <I::Scalar as Zero>::zero();

        let mut difmax = zero;
        let mut sum2 = zero;
        let mut sum3 = zero;

        let mut sum2_abs = 0.0;
        let mut sum3_abs = 0.0;

        let mut coordinate = centres;
        let mut abscissa = [0.0; NDIM];

        let mut largest_error_axis = 0;

        for j in 0..NDIM {
            abscissa[j] = Self::lambda_2() * widths[j];

            coordinate[j] = centres[j] - abscissa[j];
            let fval_1 = self.function.evaluate(&coordinate);

            coordinate[j] = centres[j] + abscissa[j];
            let fval_2 = self.function.evaluate(&coordinate);

            sum2 += fval_1 + fval_2;
            sum2_abs += fval_1.abs() + fval_2.abs();

            abscissa[j] = Self::lambda_4() * widths[j];

            coordinate[j] = centres[j] - abscissa[j];
            let fval_3 = self.function.evaluate(&coordinate);

            coordinate[j] = centres[j] + abscissa[j];
            let fval_4 = self.function.evaluate(&coordinate);

            sum3 += fval_3 + fval_4;
            sum3_abs += fval_3.abs() + fval_4.abs();

            let dif = (fval_1 + fval_2) * 7.0 - (fval_3 + fval_4) - sum1 * 12.0;

            // added .abs() for testing
            if dif.abs() >= difmax.abs() {
                difmax = dif;
                // adaptive_integrator_multidin.cpp has j + 1;
                largest_error_axis = j;
            }

            coordinate[j] = centres[j];
        }

        let mut sum4 = zero;
        let mut sum4_abs = 0.0;
        for j in 1..NDIM {
            let j1 = j - 1;
            for k in j..NDIM {
                for _l in 0..2 {
                    abscissa[j1] = -abscissa[j1];
                    coordinate[j1] = centres[j1] + abscissa[j1];
                    for _m in 0..2 {
                        abscissa[k] = -abscissa[k];
                        coordinate[k] = centres[k] + abscissa[k];
                        let fval_4 = self.function.evaluate(&coordinate);
                        sum4 += fval_4;
                        sum4_abs += fval_4.abs();
                    }
                }
                coordinate[k] = centres[k];
            }
            coordinate[j1] = centres[j1];
        }

        let mut sum5 = zero;
        let mut sum5_abs = 0.0;

        for j in 0..NDIM {
            abscissa[j] = -Self::lambda_5() * widths[j];
            coordinate[j] = centres[j] + abscissa[j];
        }

        let mut continue_loop = true;
        'outer: while continue_loop {
            let fval_5 = self.function.evaluate(&coordinate);
            sum5 += fval_5;
            sum5_abs += fval_5.abs();
            for j in 0..NDIM {
                abscissa[j] = -abscissa[j];
                coordinate[j] = centres[j] + abscissa[j];
                if abscissa[j].is_sign_positive() {
                    continue 'outer;
                }
            }
            continue_loop = false;
        }

        let compare = (sum1 * Self::prime_weight_1()
            + sum2 * Self::prime_weight_2()
            + sum3 * Self::prime_weight_3()
            + sum4 * Self::prime_weight_4())
            * volume;

        let result = (sum1 * Self::weight_1()
            + sum2 * Self::weight_2()
            + sum3 * Self::weight_3()
            + sum4 * Self::weight_4()
            + sum5 * Self::weight_5())
            * volume;

        let compare_abs = volume
            * (Self::prime_weight_1() * sum1_abs
                + Self::prime_weight_2() * sum2_abs
                + Self::prime_weight_3() * sum3_abs
                + Self::prime_weight_4() * sum4_abs);

        let result_abs = volume
            * (Self::weight_1() * sum1_abs
                + Self::weight_2() * sum2_abs
                + Self::weight_3() * sum3_abs
                + Self::weight_4() * sum4_abs
                + Self::weight_5() * sum5_abs);

        let error = (result - compare).abs();
        let error_abs = (result_abs - compare_abs).abs();

        let function_evaluations = Self::minimum_function_calls();

        MultiDimensionalRegion::unevaluated()
            .with_error(error)
            .with_result(result)
            .with_result_abs(result_abs)
            .with_error_abs(error_abs)
            .with_limits(self.limits)
            .with_largest_error_axis(largest_error_axis)
            .with_function_evaluations(function_evaluations)
    }

    const LAMBDA_2: f64 = 0.358_568_582_800_318_073;
    const LAMBDA_4: f64 = 0.948_683_298_050_513_796;
    const LAMBDA_5: f64 = 0.688_247_201_611_685_289;

    const fn lambda_2() -> f64 {
        Self::LAMBDA_2
    }

    const fn lambda_4() -> f64 {
        Self::LAMBDA_4
    }

    const fn lambda_5() -> f64 {
        Self::LAMBDA_5
    }

    const fn weight_1() -> f64 {
        Self::WEIGHT_1[NDIM - 2]
    }

    const fn weight_2() -> f64 {
        Self::WEIGHT_2
    }

    const fn weight_3() -> f64 {
        Self::WEIGHT_3[NDIM - 2]
    }

    const fn weight_4() -> f64 {
        Self::WEIGHT_4
    }

    const fn weight_5() -> f64 {
        Self::WEIGHT_5[NDIM - 2]
    }

    // The weights w1 for fixed dimension NDIM, where w1_n = WEIGHT_1_DIM_MINUS_2[NDIM-2]
    // w1 / 2^n
    const WEIGHT_1: [f64; 14] = [
        -0.193_872_885_230_909_911,
        -0.555_606_360_818_980_835,
        -0.876_695_625_666_819_078,
        -1.157_140_679_774_424_59,
        -1.396_941_523_141_797_43,
        -1.596_098_155_768_937_54,
        -1.754_610_577_655_844_94,
        -1.872_478_788_802_519_83,
        -1.949_702_789_208_962_01,
        -1.986_282_578_875_171_46,
        -1.982_218_157_801_148_18,
        -1.937_509_525_986_892_19,
        -1.852_156_683_432_403_47,
        -1.726_159_630_137_682_25,
    ];

    //w2 / 2^n
    const WEIGHT_2: f64 = 980.0 / 6561.0;

    // The weights w3 for fixed dimension NDIM, where w3_n = WEIGHT_1_DIM_MINUS_2[NDIM-2]
    // w3 / 2^n
    const WEIGHT_3: [f64; 14] = [
        0.051_821_368_693_796_676_8,
        0.031_499_263_323_680_333_0,
        0.011_177_157_953_563_989_1,
        -0.009_144_947_416_552_354_73,
        -0.029_467_052_786_668_698_6,
        -0.049_789_158_156_785_042_4,
        -0.070_111_263_526_901_376_8,
        -0.090_433_368_897_017_724_1,
        -0.110_755_474_267_134_071,
        -0.131_077_579_637_250_419,
        -0.151_399_685_007_366_752,
        -0.171_721_790_377_483_099,
        -0.192_043_895_747_599_447,
        -0.212_366_001_117_715_794,
    ];

    //w4 / 2^n
    const WEIGHT_4: f64 = 200.0 / 19683.0;

    // The weights w5 for fixed dimension NDIM, where w5_n = WEIGHT_1_DIM_MINUS_2[NDIM-2]
    // w5 / 2^n
    const WEIGHT_5: [f64; 14] = [
        0.871_183_254_585_174_982e-01,
        0.435_591_627_292_587_508e-01,
        0.217_795_813_646_293_754e-01,
        0.108_897_906_823_146_873e-01,
        0.544_489_534_115_734_364e-02,
        0.272_244_767_057_867_193e-02,
        0.136_122_383_528_933_596e-02,
        0.680_611_917_644_667_955e-03,
        0.340_305_958_822_333_977e-03,
        0.170_152_979_411_166_995e-03,
        0.850_764_897_055_834_977e-04,
        0.425_382_448_527_917_472e-04,
        0.212_691_224_263_958_736e-04,
        0.106_345_612_131_979_372e-04,
    ];

    const fn prime_weight_1() -> f64 {
        Self::PRIME_WEIGHT_1[NDIM - 2]
    }

    const fn prime_weight_2() -> f64 {
        Self::PRIME_WEIGHT_2
    }

    const fn prime_weight_3() -> f64 {
        Self::PRIME_WEIGHT_3[NDIM - 2]
    }

    const fn prime_weight_4() -> f64 {
        Self::PRIME_WEIGHT_4
    }

    // w1' / 2^n
    const PRIME_WEIGHT_1: [f64; 14] = [
        -1.331_961_591_220_850_45,
        -2.292_181_069_958_847_63,
        -3.115_226_337_448_559_59,
        -3.801_097_393_689_986_11,
        -4.349_794_238_683_127_42,
        -4.761_316_872_427_983_52,
        -5.035_665_294_924_554_17,
        -5.172_839_506_172_839_39,
        -5.172_839_506_172_839_39,
        -5.035_665_294_924_554_17,
        -4.761_316_872_427_983_52,
        -4.349_794_238_683_127_42,
        -3.801_097_393_689_986_11,
        -3.115_226_337_448_559_59,
    ];

    //w2' / 2^n
    const PRIME_WEIGHT_2: f64 = 245.0 / 486.0;

    // w3' / 2^n
    const PRIME_WEIGHT_3: [f64; 14] = [
        0.044_581_618_655_692_729_2,
        -0.024_005_486_968_449_930_9,
        -0.092_592_592_592_592_587_5,
        -0.161_179_698_216_735_251,
        -0.229_766_803_840_877_915,
        -0.298_353_909_465_020_564,
        -0.366_941_015_089_163_228,
        -0.435_528_120_713_305_891,
        -0.504_115_226_337_448_555,
        -0.572_702_331_961_591_218,
        -0.641_289_437_585_733_882,
        -0.709_876_543_209_876_532,
        -0.778_463_648_834_019_195,
        -0.847_050_754_458_161_859,
    ];

    //w4' / 2^n
    const PRIME_WEIGHT_4: f64 = 25.0 / 729.0;
}

struct Geometry<const NDIM: usize> {
    centres: [f64; NDIM],
    widths: [f64; NDIM],
    volume: f64,
}

const fn minimum_function_calls(n: usize) -> usize {
    let two_pow_n = two_pow_n(n);

    two_pow_n + 2 * n * (n + 1) + 1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_geometry() {
        let lim1 = Limits::new(0.0, 1.0);
        let lim2 = Limits::new(0.0, 2.0);
        let lim3 = Limits::new(1.0, 2.0);

        let limits = [lim1, lim2, lim3];

        struct Foo {
            alpha: f64,
        }

        impl MultiDimensionalIntegrand<3> for Foo {
            type Scalar = f64;
            fn evaluate(&self, coordinates: &[f64; 3]) -> Self::Scalar {
                let x = coordinates[0];
                let y = coordinates[1];
                let z = coordinates[2];
                self.alpha * (x.powi(2) + y.powi(2) + z.powi(2))
            }
        }

        let foo = Foo { alpha: 1.0 };

        let integral = MultiDimensionalBasic::new(limits, &foo).unwrap();

        let Geometry {
            centres,
            widths,
            volume,
        } = integral.geometry();

        let centres_should_be = [0.5, 1.0, 1.5];
        let widths_should_be = [0.5, 1.0, 0.5];
        let volume_should_be = 2.0;

        assert_eq!(centres, centres_should_be);
        assert_eq!(widths, widths_should_be);
        assert_eq!(volume, volume_should_be);
    }

    #[test]
    fn test_error_wrong_dimensionality() {
        {
            let lim1 = Limits::new(0.0, 1.0);

            let limits = [lim1];

            struct Foo {
                alpha: f64,
            }

            impl MultiDimensionalIntegrand<1> for Foo {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; 1]) -> Self::Scalar {
                    let x = coordinates[0];
                    self.alpha * x.powi(2)
                }
            }

            let foo = Foo { alpha: 1.0 };

            let integral = MultiDimensionalBasic::new(limits, &foo);

            assert!(integral.is_err());
        }

        {
            let lim1 = Limits::new(0.0, 1.0);

            let limits = [lim1; 16];

            struct Bar {
                alpha: f64,
            }

            impl MultiDimensionalIntegrand<16> for Bar {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; 16]) -> Self::Scalar {
                    let x = coordinates[0];
                    self.alpha * x.powi(2)
                }
            }

            let bar = Bar { alpha: 1.0 };

            let integral = MultiDimensionalBasic::new(limits, &bar);

            assert!(integral.is_err());
        }
    }

    #[test]
    fn simple_polynomial_multidimensional_integrals() {
        {
            struct Simple;

            impl MultiDimensionalIntegrand<2> for Simple {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; 2]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    x * x + y * y + 2.0 * x * y
                }
            }

            let limit = Limits::new(-1.0, 1.0);
            let limits = [limit; 2];

            let function = Simple;
            let integral = MultiDimensionalBasic::new(limits, function).unwrap();

            let integral_result = integral.integrate_internal();

            let result = integral_result.result();
            let should_be = 8.0 / 3.0;

            assert!((result - should_be).abs() / should_be.abs() < 1e-10);
        }

        {
            struct Simple;

            impl MultiDimensionalIntegrand<2> for Simple {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; 2]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    x * x + y * y + x * y
                }
            }

            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; 2];

            let function = Simple;
            let integral = MultiDimensionalBasic::new(limits, function).unwrap();

            let integral_result = integral.integrate_internal();

            let result = integral_result.result();
            let should_be = 11.0 / 12.0;

            assert!((result - should_be).abs() / should_be.abs() < 1e-10);
        }

        {
            struct Simple;

            impl MultiDimensionalIntegrand<2> for Simple {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; 2]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    x.powi(4) + y.powi(2) + x.powi(2) * y.powi(2)
                }
            }

            let limits = [Limits::new(-2.0, 2.0), Limits::new(-3.0, 3.0)];

            let function = Simple;
            let integral = MultiDimensionalBasic::new(limits, function).unwrap();

            let integral_result = integral.integrate_internal();

            let result = integral_result.result();
            let should_be = 1224.0 / 5.0;

            assert!((result - should_be).abs() / should_be.abs() < 1e-10);
        }
    }

    #[test]
    fn test_polynomial_with_exponentials() {
        {
            struct Simple;

            impl MultiDimensionalIntegrand<2> for Simple {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; 2]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    x * f64::exp(-x * y)
                }
            }

            let limits = [Limits::new(0.0, 2.0), Limits::new(0.0, 3.0)];

            let function = Simple;
            let integral = MultiDimensionalBasic::new(limits, function).unwrap();

            let integral_result = integral.integrate_internal();

            let result = integral_result.result();
            let should_be = 1.667_492_917_392_222;

            assert!((result - should_be).abs() / should_be.abs() < 1e-3);
        }
    }

    #[test]
    fn min_function_calls() {
        let n = 2;
        assert!(minimum_function_calls(n) == 17);
        let n = 3;
        assert!(minimum_function_calls(n) == 33);
        let n = 4;
        assert!(minimum_function_calls(n) == 57);
        let n = 5;
        assert!(minimum_function_calls(n) == 93);
    }

    #[test]
    fn polynomial_fraction_multidimensional_integrals() {
        {
            struct Simple;

            impl MultiDimensionalIntegrand<3> for Simple {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; 3]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let numerator = x.powi(2) + y.powi(2) + z.powi(2);
                    let denominator = 1.0 + x * y * z;
                    numerator / denominator
                }
            }

            let limits = [
                Limits::new(0.0, 2.0),
                Limits::new(0.0, 3.0),
                Limits::new(0.0, 1.0),
            ];

            let function = Simple;
            let integral = MultiDimensionalBasic::new(limits, function).unwrap();

            let integral_result = integral.integrate_internal();

            let result = integral_result.result();
            let should_be = 16.174_761_084;

            assert!((result - should_be).abs() / should_be.abs() < 1e-3);
        }

        {
            struct Simple;

            impl MultiDimensionalIntegrand<3> for Simple {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; 3]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let numerator = x.powi(2) + y.powi(2) + z.powi(2);
                    let denominator = 1.0 + x.powi(2);
                    numerator / denominator
                }
            }

            let limits = [
                Limits::new(0.0, 2.0),
                Limits::new(0.0, 3.0),
                Limits::new(0.0, 1.0),
            ];

            let function = Simple;
            let integral = MultiDimensionalBasic::new(limits, function).unwrap();

            let integral_result = integral.integrate_internal();

            let result = integral_result.result();
            let should_be = 13.750_041_016;

            assert!((result - should_be).abs() / should_be.abs() < 1e-4);
        }
        {
            struct Simple;

            impl MultiDimensionalIntegrand<3> for Simple {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; 3]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let numerator = x.powi(2) + x * y * z + z.powi(2);
                    let denominator = 1.0 + y.powi(2);
                    numerator / denominator
                }
            }

            let limits = [
                Limits::new(0.0, 2.0),
                Limits::new(0.0, 3.0),
                Limits::new(0.0, 1.0),
            ];

            let function = Simple;
            let integral = MultiDimensionalBasic::new(limits, function).unwrap();

            let integral_result = integral.integrate_internal();

            let result = integral_result.result();
            let should_be = 5.314_778_459;

            assert!((result - should_be).abs() / should_be.abs() < 1e-2);
        }
    }
}
