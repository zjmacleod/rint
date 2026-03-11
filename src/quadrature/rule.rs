/// A one-dimensional Gauss-Kronrod integration rule.
///
/// A Gaussian numerical integration rule approximates an integral of a function by performing a
/// weighted sum of the function evaluated at defined points/abscissae. The order of an integration
/// rule, $n$, denotes the number of abscissae, $x_{i}$, at which the function is evaluated and the
/// number of weights $w_{i}$ for the weighted sum, such that the approximation is,
/// $$
/// I = \int_{b}^{a} f(x) dx \approx \sum_{i = 1}^{n} W_{i} f(X_{i}) = I_{n}
/// $$
/// where the $X_{i}$ and $W_{i}$ are the rescaled abscissae and weights,
/// $$
/// X_{i} = \frac{b + a + (a - b) x_{i}}{2} ~~~~~~~~ W_{i} = \frac{(a - b) w_{i}}{2}
/// $$
/// A Gauss-Kronrod integration rule combines two rules of different order for efficient estimation
/// of the numerical error. The rules for an $n$-point Gauss-Kronrod rule contain $m = (n - 1) / 2$
/// abscissae _shared_ by the Gaussian and Kronrod rules and an extended set of $n - m$ Kronrod
/// abscissae. The weighted sum of the full set of $n$ Kronrod function evaluations are used to
/// approximate the result of the integration, while the weighted sum of the lower order set of $m$
/// Gaussian points are used to calculate the numerical error in the routine,
/// $$
/// E = |I_{n} - I_{m}|
/// $$
/// This approach is efficient, as only $n$ total function evaluations are required to obtain the
/// result approximation and error estimate.
///
/// The [`Rule`] struct defines a Gauss-Kronrod quadrature rule for use in the one-dimensional
/// numerical integration routines [`Basic`], [`Adaptive`], and [`AdaptiveSingularity`]. Rules of varying
/// order $n$ are generated through dedicated constructor functions [`Rule::gk*`].
///
/// [`Basic`]: crate::quadrature::Basic
/// [`Adaptive`]: crate::quadrature::Adaptive
/// [`AdaptiveSingularity`]: crate::quadrature::AdaptiveSingularity
/// [`Rule::gk*`]: struct.Rule.html#impl-Rule
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Rule {
    evaluations: usize,
    kronrod_centre: f64,
    gauss_centre: Option<f64>,
    shared_data: &'static [SharedData],
    extended_data: &'static [ExtendedData],
}

/// Constructors
impl Rule {
    /// 15-point Gauss-Kronrod rule
    ///
    /// Generates the evaluation points/abscissae $x_{i}$ and weights $w_{i}$ for a 15-point Gauss-Kronrod
    /// integration rule.
    #[must_use]
    pub const fn gk15() -> Self {
        Self {
            evaluations: gk15::EVALUATIONS,
            kronrod_centre: gk15::KRONROD_CENTRE,
            gauss_centre: gk15::GAUSS_CENTRE,
            shared_data: &gk15::SHARED_DATA,
            extended_data: &gk15::EXTENDED_DATA,
        }
    }

    /// 21-point Gauss-Kronrod rule
    ///
    /// Generates the evaluation points/abscissae $x_{i}$ and weights $w_{i}$ for a 21-point Gauss-Kronrod
    /// integration rule.
    #[must_use]
    pub const fn gk21() -> Self {
        Self {
            evaluations: gk21::EVALUATIONS,
            kronrod_centre: gk21::KRONROD_CENTRE,
            gauss_centre: gk21::GAUSS_CENTRE,
            shared_data: &gk21::SHARED_DATA,
            extended_data: &gk21::EXTENDED_DATA,
        }
    }

    /// 31-point Gauss-Kronrod rule
    ///
    /// Generates the evaluation points/abscissae $x_{i}$ and weights $w_{i}$ for a 31-point Gauss-Kronrod
    /// integration rule.
    #[must_use]
    pub const fn gk31() -> Self {
        Self {
            evaluations: gk31::EVALUATIONS,
            kronrod_centre: gk31::KRONROD_CENTRE,
            gauss_centre: gk31::GAUSS_CENTRE,
            shared_data: &gk31::SHARED_DATA,
            extended_data: &gk31::EXTENDED_DATA,
        }
    }

    /// 41-point Gauss-Kronrod rule
    ///
    /// Generates the evaluation points/abscissae $x_{i}$ and weights $w_{i}$ for a 41-point Gauss-Kronrod
    /// integration rule.
    #[must_use]
    pub const fn gk41() -> Self {
        Self {
            evaluations: gk41::EVALUATIONS,
            kronrod_centre: gk41::KRONROD_CENTRE,
            gauss_centre: gk41::GAUSS_CENTRE,
            shared_data: &gk41::SHARED_DATA,
            extended_data: &gk41::EXTENDED_DATA,
        }
    }

    /// 51-point Gauss-Kronrod rule
    ///
    /// Generates the evaluation points/abscissae $x_{i}$ and weights $w_{i}$ for a 51-point Gauss-Kronrod
    /// integration rule.
    #[must_use]
    pub const fn gk51() -> Self {
        Self {
            evaluations: gk51::EVALUATIONS,
            kronrod_centre: gk51::KRONROD_CENTRE,
            gauss_centre: gk51::GAUSS_CENTRE,
            shared_data: &gk51::SHARED_DATA,
            extended_data: &gk51::EXTENDED_DATA,
        }
    }

    /// 61-point Gauss-Kronrod rule
    ///
    /// Generates the evaluation points/abscissae $x_{i}$ and weights $w_{i}$ for a 61-point Gauss-Kronrod
    /// integration rule.
    #[must_use]
    pub const fn gk61() -> Self {
        Self {
            evaluations: gk61::EVALUATIONS,
            kronrod_centre: gk61::KRONROD_CENTRE,
            gauss_centre: gk61::GAUSS_CENTRE,
            shared_data: &gk61::SHARED_DATA,
            extended_data: &gk61::EXTENDED_DATA,
        }
    }
}

impl Rule {
    /// The number of function evaluations required by the rule.
    pub(crate) const fn evaluations(&self) -> usize {
        self.evaluations
    }

    /// The Kronrod rules are all of odd order, and so have an abscissa/weight at the centre.
    pub(crate) const fn kronrod_centre(&self) -> f64 {
        self.kronrod_centre
    }

    /// The Gaussian rules can be of even _or_ odd order, and so conditionally have an
    /// abscissa/weight at the centre.
    pub(crate) const fn gauss_centre(&self) -> Option<f64> {
        self.gauss_centre
    }

    /// Return a slice corresponding to the abscissae shared between the Gaussian and Kronrod
    /// integration rules.
    pub(crate) const fn shared_data(&self) -> &[SharedData] {
        self.shared_data
    }

    /// Return a slice corresponding to only the extended points of the Kronrad rule.
    pub(crate) const fn extended_data(&self) -> &[ExtendedData] {
        self.extended_data
    }
}

/// Data shared by the Gaussian and Kronrod rules.
///
/// Each [`SharedData`] contains a field `point` which corresponds to the abscissa at which the
/// function is to be evaluated, a field `gauss` which is the Gaussian rule weight associated with
/// the abscissa, and a field `kronrod` which is the Kronrod rule weight associated with the
/// abscissa. The data are contained in the [`Rule`] in this form for most efficient use of Rust's
/// [`Iterator`] trait methods in calculations.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub(crate) struct SharedData {
    /// The abscissa at which the function is to be evaluated $x_{i}$.
    point: f64,
    /// The weight associated with the Gaussian integration rule $w^{G}_{i}$ at the point $x_{i}$.
    gauss: f64,
    /// The weight associated with the Kronrod integration rule $w^{K}_{i}$ at the point $x_{i}$.
    kronrod: f64,
}

impl SharedData {
    /// Generate a new instance of [`SharedData`].
    pub(crate) const fn new(point: f64, gauss: f64, kronrod: f64) -> Self {
        Self {
            point,
            gauss,
            kronrod,
        }
    }

    /// Return the abscissa/point at which the function is to be evaluated.
    pub(crate) const fn point(&self) -> f64 {
        self.point
    }

    /// Return the Gaussian weight corresponding to the shared data absicssa.
    pub(crate) const fn gauss(&self) -> f64 {
        self.gauss
    }

    /// Return the Kronrod weight corresponding to the shared data absicssa.
    pub(crate) const fn kronrod(&self) -> f64 {
        self.kronrod
    }
}

/// Data unique to the extended Kronrod rule.
///
/// Each [`ExtendedData`] contains a field `point` which corresponds to the abscissa at which the
/// function is to be evaluated and a field `kronrod` which is the Kronrod rule weight associated
/// with the abscissa. The data are contained in the [`Rule`] in this form for most efficient use
/// of Rust's [`Iterator`] trait methods in calculations.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub(crate) struct ExtendedData {
    point: f64,
    kronrod: f64,
}

impl ExtendedData {
    /// Generate a new instance of [`ExtendedData`].
    pub(crate) const fn new(point: f64, kronrod: f64) -> Self {
        Self { point, kronrod }
    }

    /// Return the abscissa/point at which the function is to be evaluated.
    pub(crate) const fn point(&self) -> f64 {
        self.point
    }

    /// Return the Kronrod weight corresponding to the extended data absicssa.
    pub(crate) const fn kronrod(&self) -> f64 {
        self.kronrod
    }
}

/// Data for the 15-point Gauss-Kronrod integration rule.
mod gk15 {
    use super::{ExtendedData, SharedData};
    pub(super) const KRONROD_CENTRE: f64 = 0.209_482_141_084_727_828_012_999_174_891_714;

    pub(super) const GAUSS_CENTRE: Option<f64> =
        Some(0.417_959_183_673_469_387_755_102_040_816_327);

    pub(super) const EVALUATIONS: usize = 15;

    pub(super) const SHARED_DATA: [SharedData; 3] = [
        SharedData::new(
            0.949_107_912_342_758_524_526_189_684_047_851,
            0.129_484_966_168_869_693_270_611_432_679_082,
            0.063_092_092_629_978_553_290_700_663_189_204,
        ),
        SharedData::new(
            0.741_531_185_599_394_439_863_864_773_280_788,
            0.279_705_391_489_276_667_901_467_771_423_780,
            0.140_653_259_715_525_918_745_189_590_510_238,
        ),
        SharedData::new(
            0.405_845_151_377_397_166_906_606_412_076_961,
            0.381_830_050_505_118_944_950_369_775_488_975,
            0.190_350_578_064_785_409_913_256_402_421_014,
        ),
    ];

    pub(super) const EXTENDED_DATA: [ExtendedData; 4] = [
        ExtendedData::new(
            0.991_455_371_120_812_639_206_854_697_526_329,
            0.022_935_322_010_529_224_963_732_008_058_970,
        ),
        ExtendedData::new(
            0.864_864_423_359_769_072_789_712_788_640_926,
            0.104_790_010_322_250_183_839_876_322_541_518,
        ),
        ExtendedData::new(
            0.586_087_235_467_691_130_294_144_838_258_730,
            0.169_004_726_639_267_902_826_583_426_598_550,
        ),
        ExtendedData::new(
            0.207_784_955_007_898_467_600_689_403_773_245,
            0.204_432_940_075_298_892_414_161_999_234_649,
        ),
    ];
}

/// Data for the 21-point Gauss-Kronrod integration rule.
mod gk21 {
    use super::{ExtendedData, SharedData};
    pub(super) const KRONROD_CENTRE: f64 = 0.149_445_554_002_916_905_664_936_468_389_821;

    pub(super) const GAUSS_CENTRE: Option<f64> = None;

    pub(super) const EVALUATIONS: usize = 21;

    pub(super) const SHARED_DATA: [SharedData; 5] = [
        SharedData::new(
            0.973_906_528_517_171_720_077_964_012_084_452,
            0.066_671_344_308_688_137_593_568_809_893_332,
            0.032_558_162_307_964_727_478_818_972_459_390,
        ),
        SharedData::new(
            0.865_063_366_688_984_510_732_096_688_423_493,
            0.149_451_349_150_580_593_145_776_339_657_697,
            0.075_039_674_810_919_952_767_043_140_916_190,
        ),
        SharedData::new(
            0.679_409_568_299_024_406_234_327_365_114_874,
            0.219_086_362_515_982_043_995_534_934_228_163,
            0.109_387_158_802_297_641_899_210_590_325_805,
        ),
        SharedData::new(
            0.433_395_394_129_247_190_799_265_943_165_784,
            0.269_266_719_309_996_355_091_226_921_569_469,
            0.134_709_217_311_473_325_928_054_001_771_707,
        ),
        SharedData::new(
            0.148_874_338_981_631_210_884_826_001_129_720,
            0.295_524_224_714_752_870_173_892_994_651_338,
            0.147_739_104_901_338_491_374_841_515_972_068,
        ),
    ];

    pub(super) const EXTENDED_DATA: [ExtendedData; 5] = [
        ExtendedData::new(
            0.995_657_163_025_808_080_735_527_280_689_003,
            0.011_694_638_867_371_874_278_064_396_062_192,
        ),
        ExtendedData::new(
            0.930_157_491_355_708_226_001_207_180_059_508,
            0.054_755_896_574_351_996_031_381_300_244_580,
        ),
        ExtendedData::new(
            0.780_817_726_586_416_897_063_717_578_345_042,
            0.093_125_454_583_697_605_535_065_465_083_366,
        ),
        ExtendedData::new(
            0.562_757_134_668_604_683_339_000_099_272_694,
            0.123_491_976_262_065_851_077_958_109_831_074,
        ),
        ExtendedData::new(
            0.294_392_862_701_460_198_131_126_603_103_866,
            0.142_775_938_577_060_080_797_094_273_138_717,
        ),
    ];
}

/// Data for the 31-point Gauss-Kronrod integration rule.
mod gk31 {
    use super::{ExtendedData, SharedData};
    pub(super) const KRONROD_CENTRE: f64 = 0.101_330_007_014_791_549_017_374_792_767_493;

    pub(super) const GAUSS_CENTRE: Option<f64> =
        Some(0.202_578_241_925_561_272_880_620_199_967_519);

    pub(super) const EVALUATIONS: usize = 31;

    pub(super) const SHARED_DATA: [SharedData; 7] = [
        SharedData::new(
            0.987_992_518_020_485_428_489_565_718_586_613,
            0.030_753_241_996_117_268_354_628_393_577_204,
            0.015_007_947_329_316_122_538_374_763_075_807,
        ),
        SharedData::new(
            0.937_273_392_400_705_904_307_758_947_710_209,
            0.070_366_047_488_108_124_709_267_416_450_667,
            0.035_346_360_791_375_846_222_037_948_478_360,
        ),
        SharedData::new(
            0.848_206_583_410_427_216_200_648_320_774_217,
            0.107_159_220_467_171_935_011_869_546_685_869,
            0.053_481_524_690_928_087_265_343_147_239_430,
        ),
        SharedData::new(
            0.724_417_731_360_170_047_416_186_054_613_938,
            0.139_570_677_926_154_314_447_804_794_511_028,
            0.069_854_121_318_728_258_709_520_077_099_147,
        ),
        SharedData::new(
            0.570_972_172_608_538_847_537_226_737_253_911,
            0.166_269_205_816_993_933_553_200_860_481_209,
            0.083_080_502_823_133_021_038_289_247_286_104,
        ),
        SharedData::new(
            0.394_151_347_077_563_369_897_207_370_981_045,
            0.186_161_000_015_562_211_026_800_561_866_423,
            0.093_126_598_170_825_321_225_486_872_747_346,
        ),
        SharedData::new(
            0.201_194_093_997_434_522_300_628_303_394_596,
            0.198_431_485_327_111_576_456_118_326_443_839,
            0.099_173_598_721_791_959_332_393_173_484_603,
        ),
    ];

    pub(super) const EXTENDED_DATA: [ExtendedData; 8] = [
        ExtendedData::new(
            0.998_002_298_693_397_060_285_172_840_152_271,
            0.005_377_479_872_923_348_987_792_051_430_128,
        ),
        ExtendedData::new(
            0.967_739_075_679_139_134_257_347_978_784_337,
            0.025_460_847_326_715_320_186_874_001_019_653,
        ),
        ExtendedData::new(
            0.897_264_532_344_081_900_882_509_656_454_496,
            0.044_589_751_324_764_876_608_227_299_373_280,
        ),
        ExtendedData::new(
            0.790_418_501_442_465_932_967_649_294_817_947,
            0.062_009_567_800_670_640_285_139_230_960_803,
        ),
        ExtendedData::new(
            0.650_996_741_297_416_970_533_735_895_313_275,
            0.076_849_680_757_720_378_894_432_777_482_659,
        ),
        ExtendedData::new(
            0.485_081_863_640_239_680_693_655_740_232_351,
            0.088_564_443_056_211_770_647_275_443_693_774,
        ),
        ExtendedData::new(
            0.299_180_007_153_168_812_166_780_024_266_389,
            0.096_642_726_983_623_678_505_179_907_627_589,
        ),
        ExtendedData::new(
            0.101_142_066_918_717_499_027_074_231_447_392,
            0.100_769_845_523_875_595_044_946_662_617_570,
        ),
    ];
}

/// Data for the 41-point Gauss-Kronrod integration rule.
mod gk41 {
    use super::{ExtendedData, SharedData};
    pub(super) const KRONROD_CENTRE: f64 = 0.076_600_711_917_999_656_445_049_901_530_102;

    pub(super) const GAUSS_CENTRE: Option<f64> = None;

    pub(super) const EVALUATIONS: usize = 41;

    pub(super) const SHARED_DATA: [SharedData; 10] = [
        SharedData::new(
            0.993_128_599_185_094_924_786_122_388_471_320,
            0.017_614_007_139_152_118_311_861_962_351_853,
            0.008_600_269_855_642_942_198_661_787_950_102,
        ),
        SharedData::new(
            0.963_971_927_277_913_791_267_666_131_197_277,
            0.040_601_429_800_386_941_331_039_952_274_932,
            0.020_388_373_461_266_523_598_010_231_432_755,
        ),
        SharedData::new(
            0.912_234_428_251_325_905_867_752_441_203_298,
            0.062_672_048_334_109_063_569_506_535_187_042,
            0.031_287_306_777_032_798_958_543_119_323_801,
        ),
        SharedData::new(
            0.839_116_971_822_218_823_394_529_061_701_521,
            0.083_276_741_576_704_748_724_758_143_222_046,
            0.041_668_873_327_973_686_263_788_305_936_895,
        ),
        SharedData::new(
            0.746_331_906_460_150_792_614_305_070_355_642,
            0.101_930_119_817_240_435_036_750_135_480_350,
            0.050_944_573_923_728_691_932_707_670_050_345,
        ),
        SharedData::new(
            0.636_053_680_726_515_025_452_836_696_226_286,
            0.118_194_531_961_518_417_312_377_377_711_382,
            0.059_111_400_880_639_572_374_967_220_648_594,
        ),
        SharedData::new(
            0.510_867_001_950_827_098_004_364_050_955_251,
            0.131_688_638_449_176_626_898_494_499_748_163,
            0.065_834_597_133_618_422_111_563_556_969_398,
        ),
        SharedData::new(
            0.373_706_088_715_419_560_672_548_177_024_927,
            0.142_096_109_318_382_051_329_298_325_067_165,
            0.071_054_423_553_444_068_305_790_361_723_210,
        ),
        SharedData::new(
            0.227_785_851_141_645_078_080_496_195_368_575,
            0.149_172_986_472_603_746_787_828_737_001_969,
            0.074_582_875_400_499_188_986_581_418_362_488,
        ),
        SharedData::new(
            0.076_526_521_133_497_333_754_640_409_398_838,
            0.152_753_387_130_725_850_698_084_331_955_098,
            0.076_377_867_672_080_736_705_502_835_038_061,
        ),
    ];

    pub(super) const EXTENDED_DATA: [ExtendedData; 10] = [
        ExtendedData::new(
            0.998_859_031_588_277_663_838_315_576_545_863,
            0.003_073_583_718_520_531_501_218_293_246_031,
        ),
        ExtendedData::new(
            0.981_507_877_450_250_259_193_342_994_720_217,
            0.014_626_169_256_971_252_983_787_960_308_868,
        ),
        ExtendedData::new(
            0.940_822_633_831_754_753_519_982_722_212_443,
            0.025_882_133_604_951_158_834_505_067_096_153,
        ),
        ExtendedData::new(
            0.878_276_811_252_281_976_077_442_995_113_078,
            0.036_600_169_758_200_798_030_557_240_707_211,
        ),
        ExtendedData::new(
            0.795_041_428_837_551_198_350_638_833_272_788,
            0.046_434_821_867_497_674_720_231_880_926_108,
        ),
        ExtendedData::new(
            0.693_237_656_334_751_384_805_490_711_845_932,
            0.055_195_105_348_285_994_744_832_372_419_777,
        ),
        ExtendedData::new(
            0.575_140_446_819_710_315_342_946_036_586_425,
            0.062_653_237_554_781_168_025_870_122_174_255,
        ),
        ExtendedData::new(
            0.443_593_175_238_725_103_199_992_213_492_640,
            0.068_648_672_928_521_619_345_623_411_885_368,
        ),
        ExtendedData::new(
            0.301_627_868_114_913_004_320_555_356_858_592,
            0.073_030_690_332_786_667_495_189_417_658_913,
        ),
        ExtendedData::new(
            0.152_605_465_240_922_675_505_220_241_022_678,
            0.075_704_497_684_556_674_659_542_775_376_617,
        ),
    ];
}

/// Data for the 51-point Gauss-Kronrod integration rule.
mod gk51 {
    use super::{ExtendedData, SharedData};
    pub(super) const KRONROD_CENTRE: f64 = 0.061_580_818_067_832_935_078_759_824_240_066;

    pub(super) const GAUSS_CENTRE: Option<f64> =
        Some(0.123_176_053_726_715_451_203_902_873_079_050);

    pub(super) const EVALUATIONS: usize = 51;

    pub(super) const SHARED_DATA: [SharedData; 12] = [
        SharedData::new(
            0.995_556_969_790_498_097_908_784_946_893_902,
            0.011_393_798_501_026_287_947_902_964_113_235,
            0.005_561_932_135_356_713_758_040_236_901_066,
        ),
        SharedData::new(
            0.976_663_921_459_517_511_498_315_386_479_594,
            0.026_354_986_615_032_137_261_901_815_295_299,
            0.013_236_229_195_571_674_813_656_405_846_976,
        ),
        SharedData::new(
            0.942_974_571_228_974_339_414_011_169_658_471,
            0.040_939_156_701_306_312_655_623_487_711_646,
            0.020_435_371_145_882_835_456_568_292_235_939,
        ),
        SharedData::new(
            0.894_991_997_878_275_368_851_042_006_782_805,
            0.054_904_695_975_835_191_925_936_891_540_473,
            0.027_475_317_587_851_737_802_948_455_517_811,
        ),
        SharedData::new(
            0.833_442_628_760_834_001_421_021_108_693_570,
            0.068_038_333_812_356_917_207_187_185_656_708,
            0.034_002_130_274_329_337_836_748_795_229_551,
        ),
        SharedData::new(
            0.759_259_263_037_357_630_577_282_865_204_361,
            0.080_140_700_335_001_018_013_234_959_669_111,
            0.040_083_825_504_032_382_074_839_284_467_076,
        ),
        SharedData::new(
            0.673_566_368_473_468_364_485_120_633_247_622,
            0.091_028_261_982_963_649_811_497_220_702_892,
            0.045_502_913_049_921_788_909_870_584_752_660,
        ),
        SharedData::new(
            0.577_662_930_241_222_967_723_689_841_612_654,
            0.100_535_949_067_050_644_202_206_890_392_686,
            0.050_277_679_080_715_671_963_325_259_433_440,
        ),
        SharedData::new(
            0.473_002_731_445_714_960_522_182_115_009_192,
            0.108_519_624_474_263_653_116_093_957_050_117,
            0.054_251_129_888_545_490_144_543_370_459_876,
        ),
        SharedData::new(
            0.361_172_305_809_387_837_735_821_730_127_641,
            0.114_858_259_145_711_648_339_325_545_869_556,
            0.057_437_116_361_567_832_853_582_693_939_506,
        ),
        SharedData::new(
            0.243_866_883_720_988_432_045_190_362_797_452,
            0.119_455_763_535_784_772_228_178_126_512_901,
            0.059_720_340_324_174_059_979_099_291_932_562,
        ),
        SharedData::new(
            0.122_864_692_610_710_396_387_359_818_808_037,
            0.122_242_442_990_310_041_688_959_518_945_852,
            0.061_128_509_717_053_048_305_859_030_416_293,
        ),
    ];

    pub(super) const EXTENDED_DATA: [ExtendedData; 13] = [
        ExtendedData::new(
            0.999_262_104_992_609_834_193_457_486_540_341,
            0.001_987_383_892_330_315_926_507_851_882_843,
        ),
        ExtendedData::new(
            0.988_035_794_534_077_247_637_331_014_577_406,
            0.009_473_973_386_174_151_607_207_710_523_655,
        ),
        ExtendedData::new(
            0.961_614_986_425_842_512_418_130_033_660_167,
            0.016_847_817_709_128_298_231_516_667_536_336,
        ),
        ExtendedData::new(
            0.920_747_115_281_701_561_746_346_084_546_331,
            0.024_009_945_606_953_216_220_092_489_164_881,
        ),
        ExtendedData::new(
            0.865_847_065_293_275_595_448_996_969_588_340,
            0.030_792_300_167_387_488_891_109_020_215_229,
        ),
        ExtendedData::new(
            0.797_873_797_998_500_059_410_410_904_994_307,
            0.037_116_271_483_415_543_560_330_625_367_620,
        ),
        ExtendedData::new(
            0.717_766_406_813_084_388_186_654_079_773_298,
            0.042_872_845_020_170_049_476_895_792_439_495,
        ),
        ExtendedData::new(
            0.626_810_099_010_317_412_788_122_681_624_518,
            0.047_982_537_138_836_713_906_392_255_756_915,
        ),
        ExtendedData::new(
            0.526_325_284_334_719_182_599_623_778_158_010,
            0.052_362_885_806_407_475_864_366_712_137_873,
        ),
        ExtendedData::new(
            0.417_885_382_193_037_748_851_814_394_594_572,
            0.055_950_811_220_412_317_308_240_686_382_747,
        ),
        ExtendedData::new(
            0.303_089_538_931_107_830_167_478_909_980_339,
            0.058_689_680_022_394_207_961_974_175_856_788,
        ),
        ExtendedData::new(
            0.183_718_939_421_048_892_015_969_888_759_528,
            0.060_539_455_376_045_862_945_360_267_517_565,
        ),
        ExtendedData::new(
            0.061_544_483_005_685_078_886_546_392_366_797,
            0.061_471_189_871_425_316_661_544_131_965_264,
        ),
    ];
}

/// Data for the 61-point Gauss-Kronrod integration rule.
mod gk61 {
    use super::{ExtendedData, SharedData};
    pub(super) const KRONROD_CENTRE: f64 = 0.051_494_729_429_451_567_558_340_433_647_099;

    pub(super) const GAUSS_CENTRE: Option<f64> = None;

    pub(super) const EVALUATIONS: usize = 61;

    pub(super) const SHARED_DATA: [SharedData; 15] = [
        SharedData::new(
            0.996_893_484_074_649_540_271_630_050_918_695,
            0.007_968_192_496_166_605_615_465_883_474_674,
            0.003_890_461_127_099_884_051_267_201_844_516,
        ),
        SharedData::new(
            0.983_668_123_279_747_209_970_032_581_605_663,
            0.018_466_468_311_090_959_142_302_131_912_047,
            0.009_273_279_659_517_763_428_441_146_892_024,
        ),
        SharedData::new(
            0.960_021_864_968_307_512_216_871_025_581_798,
            0.028_784_707_883_323_369_349_719_179_611_292,
            0.014_369_729_507_045_804_812_451_432_443_580,
        ),
        SharedData::new(
            0.926_200_047_429_274_325_879_324_277_080_474,
            0.038_799_192_569_627_049_596_801_936_446_348,
            0.019_414_141_193_942_381_173_408_951_050_128,
        ),
        SharedData::new(
            0.882_560_535_792_052_681_543_116_462_530_226,
            0.048_402_672_830_594_052_902_938_140_422_808,
            0.024_191_162_078_080_601_365_686_370_725_232,
        ),
        SharedData::new(
            0.829_565_762_382_768_397_442_898_119_732_502,
            0.057_493_156_217_619_066_481_721_689_402_056,
            0.028_754_048_765_041_292_843_978_785_354_334,
        ),
        SharedData::new(
            0.767_777_432_104_826_194_917_977_340_974_503,
            0.065_974_229_882_180_495_128_128_515_115_962,
            0.032_981_447_057_483_726_031_814_191_016_854,
        ),
        SharedData::new(
            0.697_850_494_793_315_796_932_292_388_026_640,
            0.073_755_974_737_705_206_268_243_850_022_191,
            0.036_882_364_651_821_229_223_911_065_617_136,
        ),
        SharedData::new(
            0.620_526_182_989_242_861_140_477_556_431_189,
            0.080_755_895_229_420_215_354_694_938_460_530,
            0.040_374_538_951_535_959_111_995_279_752_468,
        ),
        SharedData::new(
            0.536_624_148_142_019_899_264_169_793_311_073,
            0.086_899_787_201_082_979_802_387_530_715_126,
            0.043_452_539_701_356_069_316_831_728_117_073,
        ),
        SharedData::new(
            0.447_033_769_538_089_176_780_609_900_322_854,
            0.092_122_522_237_786_128_717_632_707_087_619,
            0.046_059_238_271_006_988_116_271_735_559_374,
        ),
        SharedData::new(
            0.352_704_725_530_878_113_471_037_207_089_374,
            0.096_368_737_174_644_259_639_468_626_351_810,
            0.048_185_861_757_087_129_140_779_492_298_305,
        ),
        SharedData::new(
            0.254_636_926_167_889_846_439_805_129_817_805,
            0.099_593_420_586_795_267_062_780_282_103_569,
            0.049_795_683_427_074_206_357_811_569_379_942,
        ),
        SharedData::new(
            0.153_869_913_608_583_546_963_794_672_743_256,
            0.101_762_389_748_405_504_596_428_952_168_554,
            0.050_881_795_898_749_606_492_297_473_049_805,
        ),
        SharedData::new(
            0.051_471_842_555_317_695_833_025_213_166_723,
            0.102_852_652_893_558_840_341_285_636_705_415,
            0.051_426_128_537_459_025_933_862_879_215_781,
        ),
    ];

    pub(super) const EXTENDED_DATA: [ExtendedData; 15] = [
        ExtendedData::new(
            0.999_484_410_050_490_637_571_325_895_705_811,
            0.001_389_013_698_677_007_624_551_591_226_760,
        ),
        ExtendedData::new(
            0.991_630_996_870_404_594_858_628_366_109_486,
            0.006_630_703_915_931_292_173_319_826_369_750,
        ),
        ExtendedData::new(
            0.973_116_322_501_126_268_374_693_868_423_707,
            0.011_823_015_253_496_341_742_232_898_853_251,
        ),
        ExtendedData::new(
            0.944_374_444_748_559_979_415_831_324_037_439,
            0.016_920_889_189_053_272_627_572_289_420_322,
        ),
        ExtendedData::new(
            0.905_573_307_699_907_798_546_522_558_925_958,
            0.021_828_035_821_609_192_297_167_485_738_339,
        ),
        ExtendedData::new(
            0.857_205_233_546_061_098_958_658_510_658_944,
            0.026_509_954_882_333_101_610_601_709_335_075,
        ),
        ExtendedData::new(
            0.799_727_835_821_839_083_013_668_942_322_683,
            0.030_907_257_562_387_762_472_884_252_943_092,
        ),
        ExtendedData::new(
            0.733_790_062_453_226_804_726_171_131_369_528,
            0.034_979_338_028_060_024_137_499_670_731_468,
        ),
        ExtendedData::new(
            0.660_061_064_126_626_961_370_053_668_149_271,
            0.038_678_945_624_727_592_950_348_651_532_281,
        ),
        ExtendedData::new(
            0.579_345_235_826_361_691_756_024_932_172_540,
            0.041_969_810_215_164_246_147_147_541_285_970,
        ),
        ExtendedData::new(
            0.492_480_467_861_778_574_993_693_061_207_709,
            0.044_814_800_133_162_663_192_355_551_616_723,
        ),
        ExtendedData::new(
            0.400_401_254_830_394_392_535_476_211_542_661,
            0.047_185_546_569_299_153_945_261_478_181_099,
        ),
        ExtendedData::new(
            0.304_073_202_273_625_077_372_677_107_199_257,
            0.049_055_434_555_029_778_887_528_165_367_238,
        ),
        ExtendedData::new(
            0.204_525_116_682_309_891_438_957_671_002_025,
            0.050_405_921_402_782_346_840_893_085_653_585,
        ),
        ExtendedData::new(
            0.102_806_937_966_737_030_147_096_751_318_001,
            0.051_221_547_849_258_772_170_656_282_604_944,
        ),
    ];
}
