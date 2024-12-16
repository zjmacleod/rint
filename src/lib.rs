//! Numerical integration routines written in Rust.
//!
//! This library exposes various numerical integration routines written in pure Rust.
//! The algorithms are heavily inspired by the GNU Scientific Library & QUADPACK numerical integration routines.
//!
#![deny(clippy::pedantic)]
#![allow(clippy::module_name_repetitions, clippy::excessive_precision)]
mod limits;
pub mod multi;
pub mod quadrature;

pub use limits::Limits;

/// The integrand of a one-dimensional integral.
///
/// The [`Integrand`] trait is the main entry point of this library. It defines how the integrand,
/// i.e. function `f(x)` should be evaluated at a specified variable `x`. The trait should be
/// implemented on a struct containing all of the constant parameters required by the function.
///```rust
/// use rint::Integrand;
/// struct GaussianExponential {
///     amplitude: f64,
///     mean: f64,
///     std_dev: f64,
/// }
/// impl Integrand for GaussianExponential {
///     fn evaluate(&self, x: f64) -> f64 {
///         self.amplitude * f64::exp(-(x - self.mean).powi(2) / (2.0 * self.std_dev.powi(2)))
///     }
/// }
///```
pub trait Integrand {
    fn evaluate(&self, x: f64) -> f64;
}

impl<I: Integrand> Integrand for &I {
    fn evaluate(&self, x: f64) -> f64 {
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
