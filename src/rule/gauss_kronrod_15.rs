use crate::rule::Rule;

/// The unit struct of the 15-point Gauss-Kronrod quadrature rule.
///
/// The 15-point rule combines a 7-point Gaussian rule with a 15-point Kronrod extension.
#[derive(Clone, Copy)]
pub struct GaussKronrod15;

impl Rule for GaussKronrod15 {
    type Shared = [f64; 3];
    type Extended = [f64; 4];

    const KRONROD_CENTRE: f64 = 0.209_482_141_084_727_828_012_999_174_891_714;

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
        Some(0.417_959_183_673_469_387_755_102_040_816_327)
    }
}

static SHARED_NODES: [f64; 3] = [
    0.949_107_912_342_758_524_526_189_684_047_851,
    0.741_531_185_599_394_439_863_864_773_280_788,
    0.405_845_151_377_397_166_906_606_412_076_961,
];

static GAUSS_WEIGHTS: [f64; 3] = [
    0.129_484_966_168_869_693_270_611_432_679_082,
    0.279_705_391_489_276_667_901_467_771_423_780,
    0.381_830_050_505_118_944_950_369_775_488_975,
];

static KRONROD_WEIGHTS: [f64; 3] = [
    0.063_092_092_629_978_553_290_700_663_189_204,
    0.140_653_259_715_525_918_745_189_590_510_238,
    0.190_350_578_064_785_409_913_256_402_421_014,
];

static EXTENDED_KRONROD_NODES: [f64; 4] = [
    0.991_455_371_120_812_639_206_854_697_526_329,
    0.864_864_423_359_769_072_789_712_788_640_926,
    0.586_087_235_467_691_130_294_144_838_258_730,
    0.207_784_955_007_898_467_600_689_403_773_245,
];

static EXTENDED_KRONROD_WEIGHTS: [f64; 4] = [
    0.022_935_322_010_529_224_963_732_008_058_970,
    0.104_790_010_322_250_183_839_876_322_541_518,
    0.169_004_726_639_267_902_826_583_426_598_550,
    0.204_432_940_075_298_892_414_161_999_234_649,
];
