//! Numerical integration routines written in Rust.
//!
//! This library exposes various numerical integration routines written in pure Rust.
//! The algorithms are heavily inspired by the GNU Scientific Library & QUADPACK numerical integration routines.
//!
#![deny(clippy::pedantic)]
#![allow(
    clippy::excessive_precision,
    clippy::module_name_repetitions,
    clippy::doc_lazy_continuation,
    clippy::cast_possible_truncation,
    clippy::cast_lossless,
    clippy::cast_precision_loss,
    clippy::too_many_lines,
    clippy::struct_field_names,
    clippy::if_not_else
)]
use num_complex::{Complex, ComplexFloat};
use num_traits::Zero;
use std::fmt::Debug;
use std::ops::{Add, AddAssign, Div, Mul, Sub};

mod estimate;
mod limits;
pub mod multi;
pub mod quadrature;

pub use estimate::IntegralEstimate;
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

impl<I: Integrand> Integrand for &I {
    type Scalar = I::Scalar;

    fn evaluate(&self, x: f64) -> Self::Scalar {
        I::evaluate(self, x)
    }
}

pub trait MultiDimensionalIntegrand<const NDIM: usize> {
    type Scalar: ScalarF64;

    fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar;
}

impl<I: MultiDimensionalIntegrand<NDIM>, const NDIM: usize> MultiDimensionalIntegrand<NDIM> for &I {
    type Scalar = I::Scalar;

    fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
        I::evaluate(self, coordinates)
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
    + Add<Self>
    + Sub<Self>
    + Debug
    + sealed::Max
    + sealed::Sealed
{
}

impl ScalarF64 for f64 {}
impl ScalarF64 for Complex<f64> {}

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
