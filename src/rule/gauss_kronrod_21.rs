use crate::rule::Rule;

/// The unit struct of the 21-point Gauss-Kronrod quadrature rule.
///
/// The 21-point rule combines a 10-point Gaussian rule with a 21-point Kronrod extension.
#[derive(Clone, Copy, Debug)]
pub struct GaussKronrod21;

impl Rule for GaussKronrod21 {
    type Shared = [f64; 5];
    type Extended = [f64; 5];

    const KRONROD_CENTRE: f64 = 0.149_445_554_002_916_905_664_936_468_389_821;
    const EVALUATIONS: usize = 21;

    fn shared_nodes(&self) -> Self::Shared {
        SHARED_NODES
    }

    fn gauss_weights(&self) -> Self::Shared {
        GAUSS_WEIGHTS
    }

    fn kronrod_weights(&self) -> Self::Shared {
        KRONROD_WEIGHTS
    }

    fn extended_nodes(&self) -> Self::Extended {
        EXTENDED_KRONROD_NODES
    }

    fn extended_weights(&self) -> Self::Extended {
        EXTENDED_KRONROD_WEIGHTS
    }

    fn gauss_centre(&self) -> Option<f64> {
        None
    }
}

static SHARED_NODES: [f64; 5] = [
    0.973_906_528_517_171_720_077_964_012_084_452,
    0.865_063_366_688_984_510_732_096_688_423_493,
    0.679_409_568_299_024_406_234_327_365_114_874,
    0.433_395_394_129_247_190_799_265_943_165_784,
    0.148_874_338_981_631_210_884_826_001_129_720,
];

static GAUSS_WEIGHTS: [f64; 5] = [
    0.066_671_344_308_688_137_593_568_809_893_332,
    0.149_451_349_150_580_593_145_776_339_657_697,
    0.219_086_362_515_982_043_995_534_934_228_163,
    0.269_266_719_309_996_355_091_226_921_569_469,
    0.295_524_224_714_752_870_173_892_994_651_338,
];

static KRONROD_WEIGHTS: [f64; 5] = [
    0.032_558_162_307_964_727_478_818_972_459_390,
    0.075_039_674_810_919_952_767_043_140_916_190,
    0.109_387_158_802_297_641_899_210_590_325_805,
    0.134_709_217_311_473_325_928_054_001_771_707,
    0.147_739_104_901_338_491_374_841_515_972_068,
];

static EXTENDED_KRONROD_NODES: [f64; 5] = [
    0.995_657_163_025_808_080_735_527_280_689_003,
    0.930_157_491_355_708_226_001_207_180_059_508,
    0.780_817_726_586_416_897_063_717_578_345_042,
    0.562_757_134_668_604_683_339_000_099_272_694,
    0.294_392_862_701_460_198_131_126_603_103_866,
];

static EXTENDED_KRONROD_WEIGHTS: [f64; 5] = [
    0.011_694_638_867_371_874_278_064_396_062_192,
    0.054_755_896_574_351_996_031_381_300_244_580,
    0.093_125_454_583_697_605_535_065_465_083_366,
    0.123_491_976_262_065_851_077_958_109_831_074,
    0.142_775_938_577_060_080_797_094_273_138_717,
];

//static const double wg[5] =     /* weights of the 10-point gauss rule */
//{
//};
//
//KRONROD_WEIGHTS
//{
//};
