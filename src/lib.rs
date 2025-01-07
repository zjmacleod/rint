//! Numerical integration routines written in Rust.
//!
//! This library exposes various numerical integration routines written in pure Rust.
//! The algorithms are heavily inspired by the GNU Scientific Library & QUADPACK numerical integration routines.
//!
#![deny(clippy::pedantic)]
#![allow(clippy::excessive_precision, clippy::module_name_repetitions)]
use num_complex::{Complex, ComplexFloat};
use num_traits::Zero;
use std::fmt::Debug;
use std::ops::{AddAssign, Div, Mul};

mod limits;
//pub mod multi;
pub mod quadrature;

pub use limits::Limits;

/// The integrand of a one-dimensional integral.
///
/// The [`Integrand`] trait is the main entry point of this library. It defines how the integrand,
/// `f(x)`, should be evaluated at the real point `x`. The trait should be implemented on a
/// struct containing all of the constant parameters required by the function.
/// The integrand can be either real valued, with return type [`f64`], or complex valued, with
/// return type [`Complex<f64>`].
///```rust
/// use rint::Integrand;
/// struct GaussianExponential {
///     amplitude: f64,
///     mean: f64,
///     std_dev: f64,
/// }
/// impl Integrand for GaussianExponential {
///     type Scalar = f64;
///     fn evaluate(&self, x: f64) -> Self::Scalar {
///         self.amplitude * f64::exp(-(x - self.mean).powi(2) / (2.0 * self.std_dev.powi(2)))
///     }
/// }
///```
pub trait Integrand {
    /// The output type of the integrand.
    ///
    /// The integrand can be either real valued or complex valued, evaluating to [`f64`] or
    /// [`Complex<f64>`], respectively.
    type Scalar: ScalarF64;

    /// Evaluate the integrand at a real point `x`.
    fn evaluate(&self, x: f64) -> Self::Scalar;
}

pub(crate) mod sealed {
    use super::Complex;
    pub trait Sealed {}

    impl Sealed for f64 {}
    impl Sealed for Complex<f64> {}

    pub trait Max {
        fn max_value() -> Self;
    }

    impl Max for f64 {
        fn max_value() -> Self {
            f64::MAX
        }
    }

    impl Max for Complex<f64> {
        fn max_value() -> Self {
            Complex::new(f64::MAX / 2f64.sqrt(), f64::MAX / 2f64.sqrt())
        }
    }
}

/// A numerical scalar
///
/// The [`Integrand`] trait is implemented for both real and complex valued functions of a real
/// variable. The sealed trait [`ScalarF64`] is implemented for both [`f64`] and
/// [`Complex`].
pub trait ScalarF64:
    PartialEq
    + ComplexFloat<Real = f64>
    + Zero
    + Copy
    + Mul<f64, Output = Self>
    + Div<f64, Output = Self>
    + AddAssign<Self>
    + Debug
    + sealed::Max
    + sealed::Sealed
{
}

impl ScalarF64 for f64 {}
impl ScalarF64 for Complex<f64> {}

impl<I: Integrand> Integrand for &I {
    type Scalar = I::Scalar;

    fn evaluate(&self, x: f64) -> Self::Scalar {
        I::evaluate(self, x)
    }
}

pub trait MultiDimensionalIntegrand<const N: usize> {
    fn evaluate(&self, coordinates: &[f64; N]) -> f64;
}

impl<I: MultiDimensionalIntegrand<N>, const N: usize> MultiDimensionalIntegrand<N> for &I {
    fn evaluate(&self, coordinates: &[f64; N]) -> f64 {
        I::evaluate(self, coordinates)
    }
}
