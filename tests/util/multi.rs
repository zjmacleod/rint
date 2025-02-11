// The following test functions can be found in the paper:
// P. van Dooren & L. de Ridder, "An adaptive algorithm for numerical integration over an
// n-dimensional cube", J. Comp. App. Math., Vol. 2, (1976) 207-217

use num_complex::Complex;
use rint::Limits;
use rint::MultiDimensionalIntegrand;

// F1

const DIMF1: usize = 6;
pub const F1_TARGET: f64 = 1.434_761_888_397_263;

pub struct F1 {
    pub limits: [Limits; DIMF1],
}

impl F1 {
    pub fn new() -> Self {
        let x1 = Limits::new(0.0, 2.0);
        let x2 = Limits::new(0.0, 1.0);
        let x3 = Limits::new(0.0, std::f64::consts::PI / 2.0);
        let x4 = Limits::new(-1.0, 1.0);
        let x5 = Limits::new(-1.0, 1.0);
        let x6 = Limits::new(-1.0, 1.0);
        let limits = [x1, x2, x3, x4, x5, x6];
        Self { limits }
    }
}

impl MultiDimensionalIntegrand<DIMF1> for F1 {
    type Scalar = f64;

    fn evaluate(&self, coordinates: &[f64; DIMF1]) -> Self::Scalar {
        let [x1, x2, x3, x4, x5, x6] = coordinates;

        x1 * x2.powi(2) * x3.sin() / (4.0 + x4 + x5 + x6)
    }
}

// F2

const DIMF2: usize = 4;
pub const F2_TARGET: f64 = 5.753_641_449_035_616e-1;

pub struct F2 {
    pub limits: [Limits; DIMF2],
}

impl F2 {
    pub fn new() -> Self {
        let x1 = Limits::new(0.0, 1.0);
        let x2 = Limits::new(0.0, 1.0);
        let x3 = Limits::new(0.0, 1.0);
        let x4 = Limits::new(0.0, 2.0);
        let limits = [x1, x2, x3, x4];
        Self { limits }
    }
}

impl MultiDimensionalIntegrand<DIMF2> for F2 {
    type Scalar = f64;

    fn evaluate(&self, coordinates: &[f64; DIMF2]) -> Self::Scalar {
        let [x1, x2, x3, x4] = coordinates;

        x3.powi(2) * x4 * (x3 * x4).exp() / (x1 + x2 + 1.0).powi(2)
    }
}

// F3

const DIMF3: usize = 3;
pub const F3_TARGET: f64 = 2.152_142_832_595_894;

pub struct F3 {
    pub limits: [Limits; DIMF3],
}

impl F3 {
    pub fn new() -> Self {
        let x1 = Limits::new(0.0, 1.0);
        let x2 = Limits::new(0.0, 1.0);
        let x3 = Limits::new(0.0, 1.0);
        let limits = [x1, x2, x3];
        Self { limits }
    }
}

impl MultiDimensionalIntegrand<DIMF3> for F3 {
    type Scalar = f64;

    fn evaluate(&self, coordinates: &[f64; DIMF3]) -> Self::Scalar {
        let [x1, x2, x3] = coordinates;

        8.0 / (1.0 + 2.0 * (x1 + x2 + x3))
    }
}

// F4

const DIMF4: usize = 5;
pub const F4_TARGET: f64 = 16.0;

pub struct F4 {
    pub limits: [Limits; DIMF4],
}

impl F4 {
    pub fn new() -> Self {
        let x1 = Limits::new(0.0, std::f64::consts::PI);
        let x2 = Limits::new(0.0, std::f64::consts::PI);
        let x3 = Limits::new(0.0, std::f64::consts::PI);
        let x4 = Limits::new(0.0, std::f64::consts::PI);
        let x5 = Limits::new(0.0, std::f64::consts::FRAC_PI_2);
        let limits = [x1, x2, x3, x4, x5];
        Self { limits }
    }
}

impl MultiDimensionalIntegrand<DIMF4> for F4 {
    type Scalar = f64;

    fn evaluate(&self, coordinates: &[f64; DIMF4]) -> Self::Scalar {
        let [x1, x2, x3, x4, x5] = coordinates;

        (x1 + x2 + x3 + x4 + x5).cos()
    }
}

// F5

const DIMF5: usize = 4;
pub const F5_TARGET: f64 = 1.839_071_529_076_452e-1;

pub struct F5 {
    pub limits: [Limits; DIMF5],
}

impl F5 {
    pub fn new() -> Self {
        let x1 = Limits::new(0.0, 1.0);
        let x2 = Limits::new(0.0, 1.0);
        let x3 = Limits::new(0.0, 1.0);
        let x4 = Limits::new(0.0, 1.0);
        let limits = [x1, x2, x3, x4];
        Self { limits }
    }
}

impl MultiDimensionalIntegrand<DIMF5> for F5 {
    type Scalar = f64;

    fn evaluate(&self, coordinates: &[f64; DIMF5]) -> Self::Scalar {
        let [x1, _, _, _] = coordinates;

        (10.0 * x1).sin()
    }
}

// F6

const DIMF6: usize = 2;
pub const F6_TARGET: f64 = -4.0;

pub struct F6 {
    pub limits: [Limits; DIMF6],
}

impl F6 {
    pub fn new() -> Self {
        let x1 = Limits::new(0.0, 3.0 * std::f64::consts::PI);
        let x2 = Limits::new(0.0, 3.0 * std::f64::consts::PI);
        let limits = [x1, x2];
        Self { limits }
    }
}

impl MultiDimensionalIntegrand<DIMF6> for F6 {
    type Scalar = f64;

    fn evaluate(&self, coordinates: &[f64; DIMF6]) -> Self::Scalar {
        let [x1, x2] = coordinates;

        (x1 + x2).cos()
    }
}

// F7

const DIMF7: usize = 3;
pub const F7_TARGET: f64 = 8.630_462_173_553_432e-1;

pub struct F7 {
    pub limits: [Limits; DIMF7],
}

impl F7 {
    pub fn new() -> Self {
        let x1 = Limits::new(0.0, 1.0);
        let x2 = Limits::new(0.0, 1.0);
        let x3 = Limits::new(0.0, 1.0);
        let limits = [x1, x2, x3];
        Self { limits }
    }
}

impl MultiDimensionalIntegrand<DIMF7> for F7 {
    type Scalar = f64;

    fn evaluate(&self, coordinates: &[f64; DIMF7]) -> Self::Scalar {
        let [x1, x2, x3] = coordinates;

        1.0 / (x1 + x2 + x3).powi(2)
    }
}

// F8

const DIMF8: usize = 2;
pub const F8_TARGET: f64 = 1.047_591_113_142_868;

pub struct F8 {
    pub limits: [Limits; DIMF8],
}

impl F8 {
    pub fn new() -> Self {
        let x1 = Limits::new(0.0, 1.0);
        let x2 = Limits::new(0.0, 1.0);
        let limits = [x1, x2];
        Self { limits }
    }
}

impl MultiDimensionalIntegrand<DIMF8> for F8 {
    type Scalar = f64;

    fn evaluate(&self, coordinates: &[f64; DIMF8]) -> Self::Scalar {
        let [x1, x2] = coordinates;

        605.0
            * x2
            * ((1.0 + 120.0 * (1.0 - x2))
                * ((1.0 + 120.0 * (1.0 - x2)).powi(2) + 25.0 * x1.powi(2) * x2.powi(2)))
            .powi(-1)
    }
}

// F9

const DIMF9: usize = 2;
pub const F9_TARGET: f64 = 4.991_249_442_241_215e2;

pub struct F9 {
    pub limits: [Limits; DIMF9],
}

impl F9 {
    pub fn new() -> Self {
        let x1 = Limits::new(0.0, 1.0);
        let x2 = Limits::new(0.0, 1.0);
        let limits = [x1, x2];
        Self { limits }
    }
}

impl MultiDimensionalIntegrand<DIMF9> for F9 {
    type Scalar = f64;

    fn evaluate(&self, coordinates: &[f64; DIMF9]) -> Self::Scalar {
        let [x1, x2] = coordinates;

        ((x1.powi(2) + 0.0001) * ((x2 + 0.25).powi(2) + 0.0001)).powi(-1)
    }
}

// F10

const DIMF10: usize = 2;
pub const F10_TARGET: f64 = 1.436_563_656_918_090;

pub struct F10 {
    pub limits: [Limits; DIMF10],
}

impl F10 {
    pub fn new() -> Self {
        let x1 = Limits::new(0.0, 1.0);
        let x2 = Limits::new(0.0, 1.0);
        let limits = [x1, x2];
        Self { limits }
    }
}

impl MultiDimensionalIntegrand<DIMF10> for F10 {
    type Scalar = f64;

    fn evaluate(&self, coordinates: &[f64; DIMF10]) -> Self::Scalar {
        let [x1, x2] = coordinates;

        (x1 + x2 - 1.0).abs().exp()
    }
}

// C11 = F1 + i F1

const DIMC11: usize = 6;
pub const C11_TARGET_RE: f64 = F1_TARGET;
pub const C11_TARGET_IM: f64 = F1_TARGET;

pub struct C11 {
    pub limits: [Limits; DIMC11],
}

impl C11 {
    pub fn new() -> Self {
        let x1 = Limits::new(0.0, 2.0);
        let x2 = Limits::new(0.0, 1.0);
        let x3 = Limits::new(0.0, std::f64::consts::PI / 2.0);
        let x4 = Limits::new(-1.0, 1.0);
        let x5 = Limits::new(-1.0, 1.0);
        let x6 = Limits::new(-1.0, 1.0);
        let limits = [x1, x2, x3, x4, x5, x6];
        Self { limits }
    }
}

impl MultiDimensionalIntegrand<DIMC11> for C11 {
    type Scalar = Complex<f64>;

    fn evaluate(&self, coordinates: &[f64; DIMC11]) -> Self::Scalar {
        let [x1, x2, x3, x4, x5, x6] = coordinates;

        let re = x1 * x2.powi(2) * x3.sin() / (4.0 + x4 + x5 + x6);
        let im = x1 * x2.powi(2) * x3.sin() / (4.0 + x4 + x5 + x6);

        Complex::new(re, im)
    }
}

// C37 = F3 + i F7

const DIMC37: usize = 3;
pub const C37_TARGET_RE: f64 = 2.152_142_832_595_894;
pub const C37_TARGET_IM: f64 = 8.630_462_173_553_432e-1;

pub struct C37 {
    pub limits: [Limits; DIMC37],
}

impl C37 {
    pub fn new() -> Self {
        let x1 = Limits::new(0.0, 1.0);
        let x2 = Limits::new(0.0, 1.0);
        let x3 = Limits::new(0.0, 1.0);
        let limits = [x1, x2, x3];
        Self { limits }
    }
}

impl MultiDimensionalIntegrand<DIMC37> for C37 {
    type Scalar = Complex<f64>;

    fn evaluate(&self, coordinates: &[f64; DIMC37]) -> Self::Scalar {
        let [x1, x2, x3] = coordinates;

        let re = 8.0 / (1.0 + 2.0 * (x1 + x2 + x3));
        let im = 1.0 / (x1 + x2 + x3).powi(2);

        Complex::new(re, im)
    }
}
