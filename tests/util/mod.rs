#![allow(dead_code)]

#[macro_use]
pub(crate) mod macros;

pub(crate) mod multi;

use num_complex::Complex;
use std::f64::consts::PI;

use rint::{Integrand, Limits};

/* These are the test functions from table 4.1 of the QUADPACK book */

/* f1(x) = x^alpha * log(1/x) */
/* integ(f1,x,0,1) = 1/(alpha + 1)^2 */
pub(crate) struct Function1 {
    pub(crate) alpha: f64,
}

impl Function1 {
    pub(crate) fn new(alpha: f64) -> Self {
        Self { alpha }
    }
}

impl Integrand for Function1 {
    type Point = f64;
    type Scalar = f64;
    fn evaluate(&self, x: Self::Point) -> f64 {
        let alpha = self.alpha;
        x.powf(alpha) * (1.0 / x).ln()
    }
}

/* f2(x) = 4^-alpha / ((x-pi/4)^2 + 16^-alpha) */
/* integ(f2,x,0,1) = arctan((4-pi)4^(alpha-1)) + arctan(pi 4^(alpha-1)) */
pub(crate) struct Function2 {
    pub(crate) alpha: f64,
}

impl Function2 {
    pub(crate) fn new(alpha: f64) -> Self {
        Self { alpha }
    }
}

impl Integrand for Function2 {
    type Point = f64;
    type Scalar = f64;
    fn evaluate(&self, x: Self::Point) -> f64 {
        let alpha = self.alpha;
        4.0_f64.powf(-alpha) / ((x - PI / 4.0).powi(2) + 16.0_f64.powf(-alpha))
    }
}

/* f3(x) = cos(2^alpha * sin(x)) */
/* integ(f3,x,0,pi) = pi J_0(2^alpha) */
pub(crate) struct Function3 {
    pub(crate) alpha: f64,
}

impl Function3 {
    pub(crate) fn new(alpha: f64) -> Self {
        Self { alpha }
    }
}

impl Integrand for Function3 {
    type Point = f64;
    type Scalar = f64;
    fn evaluate(&self, x: Self::Point) -> f64 {
        let alpha = self.alpha;
        (2.0_f64.powf(alpha) * x.sin()).cos()
    }
}

/* f11(x) = log(1/x)^(alpha - 1) */
/* integ(f11,x,0,1) = Gamma(alpha) */
pub(crate) struct Function11 {
    pub(crate) alpha: f64,
}

impl Function11 {
    pub(crate) fn new(alpha: f64) -> Self {
        Self { alpha }
    }
}

impl Integrand for Function11 {
    type Point = f64;
    type Scalar = f64;
    fn evaluate(&self, x: Self::Point) -> f64 {
        let alpha = self.alpha;
        f64::ln(1.0 / x).powf(alpha - 1.0)
    }
}

pub(crate) struct Function15 {
    pub(crate) alpha: i32,
}

impl Function15 {
    pub(crate) fn new(alpha: i32) -> Self {
        Self { alpha }
    }
}

impl Integrand for Function15 {
    type Point = f64;
    type Scalar = f64;
    fn evaluate(&self, x: Self::Point) -> f64 {
        let alpha = self.alpha;
        x.powi(2) * f64::exp(-2.0f64.powi(-alpha) * x)
    }
}

/* f16(x) = x^(alpha - 1) / (1 + 10 x)^2*/
pub(crate) struct Function16 {
    pub(crate) alpha: i32,
}

impl Function16 {
    pub(crate) fn new(alpha: i32) -> Self {
        Self { alpha }
    }
}

impl Integrand for Function16 {
    type Point = f64;
    type Scalar = f64;
    fn evaluate(&self, x: Self::Point) -> f64 {
        let alpha = self.alpha;
        if x.to_bits() == 0f64.to_bits() && alpha == 1 {
            1.0
        } else if x.to_bits() == 0f64.to_bits() && alpha > 1 {
            0.0
        } else {
            x.powi(alpha - 1) / (1.0 + 10.0 * x).powi(2)
        }
    }
}

/* myfn1(x) = exp(-x - x^2) */
/* integ(myfn1,x,-inf,inf) = sqrt(pi) exp(-1/4) */

pub(crate) struct MyFunciton1;

impl MyFunciton1 {
    pub(crate) fn new() -> Self {
        Self
    }
}

impl Integrand for MyFunciton1 {
    type Point = f64;
    type Scalar = f64;
    fn evaluate(&self, x: Self::Point) -> f64 {
        f64::exp(-x - x.powi(2))
    }
}

pub(crate) struct MyFunciton2 {
    pub(crate) alpha: f64,
}

impl MyFunciton2 {
    pub(crate) fn new(alpha: f64) -> Self {
        Self { alpha }
    }
}

impl Integrand for MyFunciton2 {
    type Point = f64;
    type Scalar = f64;
    fn evaluate(&self, x: Self::Point) -> f64 {
        let alpha = self.alpha;
        f64::exp(alpha * x)
    }
}

/* f455(x) = log(x)/(1+100*x^2) */
/* integ(f455,x,0,inf) = -log(10)/20 */
pub(crate) struct Function455;

impl Integrand for Function455 {
    type Point = f64;
    type Scalar = f64;
    fn evaluate(&self, x: Self::Point) -> f64 {
        x.ln() / (1.0 + 100.0 * x.powi(2))
    }
}

impl Function455 {
    pub(crate) fn new() -> Self {
        Self
    }
}

pub(crate) const CATALAN: f64 = 0.915_965_594_177_219_015_054_603_514_932_384_110_774;

fn catalan1(x: f64) -> f64 {
    -x.ln() / (1.0 + x.powi(2))
}

pub(crate) struct Catalan1 {
    lower: f64,
    upper: f64,
}

impl Catalan1 {
    pub(crate) fn new() -> Self {
        Self {
            lower: 0.0,
            upper: 1.0,
        }
    }

    pub(crate) fn limits(&self) -> Limits {
        Limits::new(self.lower, self.upper)
    }
}

impl Integrand for Catalan1 {
    type Point = f64;
    type Scalar = f64;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        catalan1(x)
    }
}

fn catalan2(x: f64) -> f64 {
    let x = x * std::f64::consts::FRAC_PI_2;
    (0.5 * x / x.sin()) * std::f64::consts::FRAC_PI_2
}

pub(crate) struct Catalan2 {
    lower: f64,
    upper: f64,
}

impl Catalan2 {
    pub(crate) fn new() -> Self {
        Self {
            lower: 0.0,
            upper: 1.0,
        }
    }

    pub(crate) fn limits(&self) -> Limits {
        Limits::new(self.lower, self.upper)
    }
}

impl Integrand for Catalan2 {
    type Point = f64;
    type Scalar = f64;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        catalan2(x)
    }
}

fn catalan3(x: f64) -> f64 {
    let x = x * std::f64::consts::FRAC_PI_4;
    (x.cos() / x.sin()).ln() * std::f64::consts::FRAC_PI_4
}

pub(crate) struct Catalan3 {
    lower: f64,
    upper: f64,
}

impl Catalan3 {
    pub(crate) fn new() -> Self {
        Self {
            lower: 0.0,
            upper: 1.0,
        }
    }

    pub(crate) fn limits(&self) -> Limits {
        Limits::new(self.lower, self.upper)
    }
}

impl Integrand for Catalan3 {
    type Point = f64;
    type Scalar = f64;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        catalan3(x)
    }
}

fn catalan4(x: f64) -> f64 {
    x.acos() / (1.0 + x.powi(2)).sqrt()
}

pub(crate) struct Catalan4 {
    lower: f64,
    upper: f64,
}

impl Catalan4 {
    pub(crate) fn new() -> Self {
        Self {
            lower: 0.0,
            upper: 1.0,
        }
    }

    pub(crate) fn limits(&self) -> Limits {
        Limits::new(self.lower, self.upper)
    }
}

impl Integrand for Catalan4 {
    type Point = f64;
    type Scalar = f64;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        catalan4(x)
    }
}

fn catalan5(x: f64) -> f64 {
    x.asinh() / (1.0 - x.powi(2)).sqrt()
}

pub(crate) struct Catalan5 {
    lower: f64,
    upper: f64,
}

impl Catalan5 {
    pub(crate) fn new() -> Self {
        Self {
            lower: 0.0,
            upper: 1.0,
        }
    }

    pub(crate) fn limits(&self) -> Limits {
        Limits::new(self.lower, self.upper)
    }
}

impl Integrand for Catalan5 {
    type Point = f64;
    type Scalar = f64;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        catalan5(x)
    }
}

fn catalan6(x: f64) -> f64 {
    (x.atanh() / (1.0 - x.powf(2.0)).sqrt()) * 0.5
}

pub(crate) struct Catalan6 {
    lower: f64,
    upper: f64,
}

impl Catalan6 {
    pub(crate) fn new() -> Self {
        Self {
            lower: 0.0,
            upper: 1.0,
        }
    }

    pub(crate) fn limits(&self) -> Limits {
        Limits::new(self.lower, self.upper)
    }
}

impl Integrand for Catalan6 {
    type Point = f64;
    type Scalar = f64;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        catalan6(x)
    }
}

fn catalan7(x: f64) -> f64 {
    x.ln() / (1.0 + x.powi(2))
}

pub(crate) struct Catalan7 {
    pub(crate) lower: f64,
}

impl Catalan7 {
    pub(crate) fn new() -> Self {
        Self { lower: 1.0 }
    }
}

impl Integrand for Catalan7 {
    type Point = f64;
    type Scalar = f64;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        catalan7(x)
    }
}

fn catalan8(x: f64) -> f64 {
    0.5 * x.atan() / x / (1.0 + x.powi(2)).sqrt()
}

pub(crate) struct Catalan8 {
    pub(crate) lower: f64,
}

impl Catalan8 {
    pub(crate) fn new() -> Self {
        Self { lower: 0.0 }
    }
}

impl Integrand for Catalan8 {
    type Point = f64;
    type Scalar = f64;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        catalan8(x)
    }
}

fn catalan9(x: f64) -> f64 {
    ((-x).exp()).atan()
}

pub(crate) struct Catalan9 {
    pub(crate) lower: f64,
}

impl Catalan9 {
    pub(crate) fn new() -> Self {
        Self { lower: 0.0 }
    }
}

impl Integrand for Catalan9 {
    type Point = f64;
    type Scalar = f64;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        catalan9(x)
    }
}

fn catalan10(x: f64) -> f64 {
    std::f64::consts::FRAC_PI_2 * (x.powi(4) - 6.0 * x.powi(2) + 1.0) * x.ln().ln()
        / (1.0 + x.powi(2)).powi(3)
}

pub(crate) struct Catalan10 {
    pub(crate) lower: f64,
}

impl Catalan10 {
    pub(crate) fn new() -> Self {
        Self { lower: 1.0 }
    }
}

impl Integrand for Catalan10 {
    type Point = f64;
    type Scalar = f64;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        catalan10(x)
    }
}

fn catalan11(x: f64) -> f64 {
    0.5 * x / x.cosh()
}

pub(crate) struct Catalan11 {
    pub(crate) lower: f64,
}

impl Catalan11 {
    pub(crate) fn new() -> Self {
        Self { lower: 0.0 }
    }
}

impl Integrand for Catalan11 {
    type Point = f64;
    type Scalar = f64;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        catalan11(x)
    }
}

pub(crate) struct ComplexCatalan12 {
    lower: f64,
    upper: f64,
}

impl ComplexCatalan12 {
    pub(crate) fn new() -> Self {
        Self {
            lower: 0.0,
            upper: 1.0,
        }
    }

    pub(crate) fn limits(&self) -> Limits {
        Limits::new(self.lower, self.upper)
    }
}

impl Integrand for ComplexCatalan12 {
    type Point = f64;
    type Scalar = Complex<f64>;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        let re = catalan1(x);
        let im = catalan2(x);
        Complex::new(re, im)
    }
}

pub(crate) struct ComplexCatalan13 {
    lower: f64,
    upper: f64,
}

impl ComplexCatalan13 {
    pub(crate) fn new() -> Self {
        Self {
            lower: 0.0,
            upper: 1.0,
        }
    }

    pub(crate) fn limits(&self) -> Limits {
        Limits::new(self.lower, self.upper)
    }
}

impl Integrand for ComplexCatalan13 {
    type Point = f64;
    type Scalar = Complex<f64>;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        let re = catalan1(x);
        let im = catalan3(x);
        Complex::new(re, im)
    }
}

pub(crate) struct ComplexCatalan14 {
    lower: f64,
    upper: f64,
}

impl ComplexCatalan14 {
    pub(crate) fn new() -> Self {
        Self {
            lower: 0.0,
            upper: 1.0,
        }
    }

    pub(crate) fn limits(&self) -> Limits {
        Limits::new(self.lower, self.upper)
    }
}

impl Integrand for ComplexCatalan14 {
    type Point = f64;
    type Scalar = Complex<f64>;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        let re = catalan1(x);
        let im = catalan4(x);
        Complex::new(re, im)
    }
}

pub(crate) struct ComplexCatalan23 {
    lower: f64,
    upper: f64,
}

impl ComplexCatalan23 {
    pub(crate) fn new() -> Self {
        Self {
            lower: 0.0,
            upper: 1.0,
        }
    }

    pub(crate) fn limits(&self) -> Limits {
        Limits::new(self.lower, self.upper)
    }
}

impl Integrand for ComplexCatalan23 {
    type Point = f64;
    type Scalar = Complex<f64>;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        let re = catalan2(x);
        let im = catalan3(x);
        Complex::new(re, im)
    }
}

pub(crate) struct ComplexCatalan24 {
    lower: f64,
    upper: f64,
}

impl ComplexCatalan24 {
    pub(crate) fn new() -> Self {
        Self {
            lower: 0.0,
            upper: 1.0,
        }
    }

    pub(crate) fn limits(&self) -> Limits {
        Limits::new(self.lower, self.upper)
    }
}

impl Integrand for ComplexCatalan24 {
    type Point = f64;
    type Scalar = Complex<f64>;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        let re = catalan2(x);
        let im = catalan4(x);
        Complex::new(re, im)
    }
}

pub(crate) struct ComplexCatalan34 {
    lower: f64,
    upper: f64,
}

impl ComplexCatalan34 {
    pub(crate) fn new() -> Self {
        Self {
            lower: 0.0,
            upper: 1.0,
        }
    }

    pub(crate) fn limits(&self) -> Limits {
        Limits::new(self.lower, self.upper)
    }
}

impl Integrand for ComplexCatalan34 {
    type Point = f64;
    type Scalar = Complex<f64>;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        let re = catalan3(x);
        let im = catalan4(x);
        Complex::new(re, im)
    }
}

pub(crate) struct ComplexCatalan15 {
    lower: f64,
    upper: f64,
}

impl ComplexCatalan15 {
    pub(crate) fn new() -> Self {
        Self {
            lower: 0.0,
            upper: 1.0,
        }
    }

    pub(crate) fn limits(&self) -> Limits {
        Limits::new(self.lower, self.upper)
    }
}

impl Integrand for ComplexCatalan15 {
    type Point = f64;
    type Scalar = Complex<f64>;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        let re = catalan1(x);
        let im = catalan5(x);
        Complex::new(re, im)
    }
}

pub(crate) struct ComplexCatalan16 {
    lower: f64,
    upper: f64,
}

impl ComplexCatalan16 {
    pub(crate) fn new() -> Self {
        Self {
            lower: 0.0,
            upper: 1.0,
        }
    }

    pub(crate) fn limits(&self) -> Limits {
        Limits::new(self.lower, self.upper)
    }
}

impl Integrand for ComplexCatalan16 {
    type Point = f64;
    type Scalar = Complex<f64>;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        let re = catalan1(x);
        let im = catalan6(x);
        Complex::new(re, im)
    }
}

pub(crate) struct ComplexCatalan45 {
    lower: f64,
    upper: f64,
}

impl ComplexCatalan45 {
    pub(crate) fn new() -> Self {
        Self {
            lower: 0.0,
            upper: 1.0,
        }
    }

    pub(crate) fn limits(&self) -> Limits {
        Limits::new(self.lower, self.upper)
    }
}

impl Integrand for ComplexCatalan45 {
    type Point = f64;
    type Scalar = Complex<f64>;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        let re = catalan4(x);
        let im = catalan5(x);
        Complex::new(re, im)
    }
}

pub(crate) struct ComplexCatalan46 {
    lower: f64,
    upper: f64,
}

impl ComplexCatalan46 {
    pub(crate) fn new() -> Self {
        Self {
            lower: 0.0,
            upper: 1.0,
        }
    }

    pub(crate) fn limits(&self) -> Limits {
        Limits::new(self.lower, self.upper)
    }
}

impl Integrand for ComplexCatalan46 {
    type Point = f64;
    type Scalar = Complex<f64>;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        let re = catalan4(x);
        let im = catalan6(x);
        Complex::new(re, im)
    }
}

pub(crate) struct ComplexCatalan56 {
    lower: f64,
    upper: f64,
}

impl ComplexCatalan56 {
    pub(crate) fn new() -> Self {
        Self {
            lower: 0.0,
            upper: 1.0,
        }
    }

    pub(crate) fn limits(&self) -> Limits {
        Limits::new(self.lower, self.upper)
    }
}

impl Integrand for ComplexCatalan56 {
    type Point = f64;
    type Scalar = Complex<f64>;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        let re = catalan5(x);
        let im = catalan6(x);
        Complex::new(re, im)
    }
}

pub(crate) struct ComplexF1F1 {
    alpha1: f64,
    alpha2: f64,
}

impl ComplexF1F1 {
    pub(crate) fn new(alpha1: f64, alpha2: f64) -> Self {
        Self { alpha1, alpha2 }
    }
}

impl Integrand for ComplexF1F1 {
    type Point = f64;
    type Scalar = Complex<f64>;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        let re = x.powf(self.alpha1) * (1.0 / x).ln();
        let im = x.powf(self.alpha2) * (1.0 / x).ln();
        Complex::new(re, im)
    }
}

pub(crate) struct ComplexF1F3 {
    alpha1: f64,
    alpha2: f64,
}

impl ComplexF1F3 {
    pub(crate) fn new(alpha1: f64, alpha2: f64) -> Self {
        Self { alpha1, alpha2 }
    }
}

impl Integrand for ComplexF1F3 {
    type Point = f64;
    type Scalar = Complex<f64>;

    fn evaluate(&self, x: Self::Point) -> Self::Scalar {
        let re = if (x > 0.0) && (x < 1.0) {
            x.powf(self.alpha1) * (1.0 / x).ln()
        } else {
            0.0
        };

        let im = if (x > 0.3) && (x < 2.71) {
            (2.0_f64.powf(self.alpha2) * x.sin()).cos()
        } else {
            0.0
        };

        Complex::new(re, im)
    }
}

pub(crate) fn test_relative_error(
    calculated: f64,
    target: f64,
    relative_error: f64,
    description: &str,
) -> Result<(), String> {
    if calculated.is_nan() || target.is_nan() {
        if calculated.is_nan() != calculated.is_nan() {
            return Err(format!(
                "Failed test {description}: calculated.is_nan() != target.is_nan()"
            ));
        }
    }

    if calculated.is_infinite() || target.is_infinite() {
        if calculated.is_infinite() != calculated.is_infinite() {
            return Err(format!(
                "Failed test {description}: calculated.is_infinite() != target.is_infinite()"
            ));
        }
    }

    if (target > 0.0 && target < f64::MIN_POSITIVE) || (target < 0.0 && target > -f64::MIN_POSITIVE)
    {
        return Err(format!(
            "Failed test {description}: target value smaller than f64::MIN_POSITIVE"
        ));
    }

    if target != 0.0 {
        if (calculated - target).abs() / target.abs() > relative_error {
            return Err(format!(
                "Failed test {description}: calculated relative error is larger than target"
            ));
        }
    } else {
        if calculated.abs() > relative_error {
            return Err(format!("Failed test {description}: target integral value was zero, but calculated integral value is larger than target relative error"));
        }
    }

    Ok(())
}

//pub(crate) fn test_absolute_error(
//    calculated: f64,
//    target: f64,
//    absolute_error: f64,
//    description: &str,
//) -> Result<(), String> {
//    if calculated.is_nan() || target.is_nan() {
//        if calculated.is_nan() != calculated.is_nan() {
//            return Err(format!(
//                "Failed test {description}: calculated.is_nan() != target.is_nan()"
//            ));
//        }
//    }
//
//    if calculated.is_infinite() || target.is_infinite() {
//        if calculated.is_infinite() != calculated.is_infinite() {
//            return Err(format!(
//                "Failed test {description}: calculated.is_infinite() != target.is_infinite()"
//            ));
//        }
//    }
//
//    if (target > 0.0 && target < f64::MIN_POSITIVE) || (target < 0.0 && target > -f64::MIN_POSITIVE)
//    {
//        return Err(format!(
//            "Failed test {description}: target value smaller than f64::MIN_POSITIVE"
//        ));
//    } else {
//        if (calculated - target).abs() > absolute_error {
//            return Err(format!("Failed test {description}: calculated absolute error is smaller than target relative error"));
//        }
//    }
//    Ok(())
//}

pub(crate) fn test_int(calculated: usize, target: usize, description: &str) -> Result<(), String> {
    if calculated == target {
        Ok(())
    } else {
        Err(format!(
            "Failed test {description}: calculated value: {calculated} target value: {target}"
        ))
    }
}
