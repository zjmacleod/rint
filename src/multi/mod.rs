//! Numerical integration routines for multi-dimensional functions.
//!
//!
//! # Overview
//!
//! This module provides numerical integration routines for approximating the $N$-dimensional
//! integrals of the form,
//! $$
//! I = \int_{\Sigma_{N}} f(\mathbf{x}) d\mathbf{x}
//! $$
//! where $\mathbf{x} = (x_{1}, x_{2}, \dots, x_{N})$ and $\Sigma_{N}$ is an $N$-dimensional
//! hypercube. The dimensionality is limited to $2 \le N \le 15$. The functions $f(x)$ can be
//! either real- or complex-valued and implement the [`MultiDimensionalIntegrand`] trait. The
//! routines are based primarily on the [DCUHRE] FORTRAN library (Bernsten, Espelid, Genz), however
//! unlike the original algorithm the routines presented in currently only operate on a single
//! function _not_ a vector of functions. The routines use fully symmetric integration [`Rule`]s,
//! with each rule of a particular order $n$ there is a set of five fully symmetric rules used,
//! where one rule of degree $n=2m+1$ is used to obtain an estimate of the integral, $R\[f\]$,
//! $$
//! I = \int_{\Sigma_{N}} f(\mathbf{x}) d\mathbf{x}
//! \approx R\[f\] = \sum_{i = 1}^{L} w_{i} f(\mathbf{x}\_{i})
//! $$
//! where $L$ is the total number of evaluation points $\mathbf{x}\_{i} = (x_{1},\dots,x_{N})$ and
//! $w_{i}$ are the rule weights. In adition there are four _null rules_ of order $2m-1$, $2m-1$,
//! $2m-3$, and $2m-5$, used to calculate,
//! $$
//! N_{j}\[f\] = \sum_{i = 1}^{L} w_{i}^{j} f(\mathbf{x}\_{i}) ~~~~~~ (j = 1,2,3,4)
//! $$
//! evaluated with the same set of points $\mathbf{x}\_{i}$, however the weights $w_{i}^{j}$ are
//! such that a null rule of degree $d$ will integrate to zero all monomials of degree $\le d$.
//! These are used in the estimation of the error.
//!
//! The algorithm and integration rules are outlined in,
//! - Bernsten, Espelid, & Genz. 1991. Algorithm 698: DCUHRE: an adaptive multidemensional
//! integration routine for a vector of integrals. ACM Trans. Math. Softw. 17, 4 (Dec. 1991),
//! 452–456. <https://doi.org/10.1145/210232.210234>
//! - Bernsten, Espelid, & Genz. 1991. An adaptive algorithm for the approximate calculation of
//! multiple integrals. ACM Trans. Math. Softw. 17, 4 (Dec. 1991), 437–451.
//! <https://doi.org/10.1145/210232.210233>
//!
//! [`MultiDimensionalIntegrand`]: crate::MultiDimensionalIntegrand
//!
//! # Available integrator routines
//!
//! The module provides two classes of routine:
//!
//! - [`Basic`]: A non-adaptive routine which applies a fully-symmetric integration [`Rule`] to a
//! function exactly once. Rules of different order are available, and are generated through the
//! `Rule*::generate` constructors of specific type alias' for each rule.
//!
//! - [`Adaptive`]: A $2 \le N \le 15$ dimensional adaptive routine. On each iteration of the
//! algorithm the axis along which the largest contribution to the error estimate was obtained is
//! used as the bisection axis to bisect the integration region and then calculate new estimates
//! for these newly bisected volumes. This concentrates the integration refinement to the regions
//! with highest error, rapidly reducing the numerical error of the routine. The algorithm uses
//! fully-symmetric integration rules, [`Rule`], of varying order and generality.
//!
//! [DCUHRE]: <https://dl.acm.org/doi/10.1145/210232.210234>
//!
//! # Examples
//!
//! ## [`Basic`] integrator example
//!
//! Here we present a calculation of [Catalan's constant] $G$ using the integral representation:
//! $$
//! G = \int_{0}^{1} \int_{0}^{1} \frac{dxdy}{1 + x^{2} y^{2}},
//! $$
//! which is a smooth integral over the integration region and can be easily integrated with the
//! [`Basic`] routine.
//!
//!```rust
//! use rint::{Limits, MultiDimensionalIntegrand, Tolerance};
//! use rint::multi::{Basic, Rule13};
//!
//! const N: usize = 2;
//! const G: f64 = 0.915_965_594_177_219_015_054_603_514_932_384_110_774;
//!
//! struct Catalan;
//!
//! impl MultiDimensionalIntegrand<N> for Catalan {
//!     type Scalar = f64;
//!     fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
//!         let [x, y] = coordinates;
//!         1.0 / (1.0 + x.powi(2) * y.powi(2))
//!     }
//! }
//!
//! # use std::error::Error;
//! # fn main() -> Result<(), Box<dyn Error>> {
//! let catalan = Catalan;
//! let limits = [Limits::new(0.0,1.0);N];
//! let rule = Rule13::generate();
//! let integral = Basic::new(&catalan, &rule, limits)?.integrate();
//!
//! let result = integral.result();
//! let error = integral.error();
//! let abs_actual_error = (G - result).abs();
//! let iters = integral.iterations();
//! assert_eq!(iters, 1);
//! assert!(abs_actual_error < error);
//! # Ok(())
//! # }
//!```
//!
//!
//! ## [`Adaptive`] integrator example
//!
//! The following example integtates a 4-dimensional function $f(\mathbf{x})$,
//! $$
//! f(\mathbf{x}) = \frac{x_{3}^{2} x_{4} e^{x_{3} x_{4}}}{(1 + x_{1} + x_{2})^{2}}
//! $$
//! where $\mathbf{x} = (x_{1}, x_{2}, x_{3}, x_{4})$ over an $N=4$ dimensional hypercube
//! $((0,1),(0,1),(0,2),(0,1))$ using a fully-symmetric 7-point adaptive algorithm.
//! Adapted from P. van Dooren & L. de Ridder, "An adaptive algorithm for numerical integration over
//! an n-dimensional cube", J. Comp. App. Math., Vol. 2, (1976) 207-217
//!
//!
//!```rust
//! use rint::{Limits, MultiDimensionalIntegrand, Tolerance};
//! use rint::multi::{Adaptive, Rule07};
//!
//! const N: usize = 4;
//!
//! struct F;
//!
//! impl MultiDimensionalIntegrand<N> for F {
//!     type Scalar = f64;
//!     fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
//!         let [x1, x2, x3, x4] = coordinates;
//!         x3.powi(2) * x4 * (x3 * x4).exp() / (x1 + x2 + 1.0).powi(2)
//!     }
//! }
//!
//! # use std::error::Error;
//! # fn main() -> Result<(), Box<dyn Error>> {
//! const TARGET: f64 = 5.753_641_449_035_616e-1;
//! const TOL: f64 = 1e-2;
//!
//! let function = F;
//! let limits = [
//!     Limits::new(0.0, 1.0),
//!     Limits::new(0.0, 1.0),
//!     Limits::new(0.0, 1.0),
//!     Limits::new(0.0, 2.0)
//! ];
//! let rule = Rule07::<N>::generate()?;
//! let tolerance = Tolerance::Relative(TOL);
//!
//! let integral = Adaptive::new(&function, &rule, limits, tolerance, 10000)?.integrate()?;
//!
//! let result = integral.result();
//! let error = integral.error();
//! let actual_error = (result - TARGET).abs();
//! let requested_error = TOL * result.abs();
//!
//! assert!(actual_error < error);
//! assert!(error < requested_error);
//! # Ok(())
//! # }
//!```
mod adaptive;
mod basic;
mod generator;
mod geometry;
mod integrator;
mod region;
pub(crate) mod rule;

pub use adaptive::Adaptive;
pub use basic::Basic;
pub(crate) use integrator::Integrator;
pub(crate) use region::Region;
pub use rule::{Rule, Rule07, Rule09, Rule09N2, Rule11, Rule13};

#[inline]
pub(crate) const fn two_pow_n(n: usize) -> usize {
    let mut exp = n;

    // Never need to check this since N > 2
    //if exp == 0 {
    //    return 1;
    //}
    let mut base = 2;
    let mut acc = 1;

    while exp > 1 {
        if (exp & 1) == 1 {
            acc *= base;
        }
        exp /= 2;
        base *= base;
    }

    acc * base
}

#[inline]
pub(crate) const fn two_pow_n_f64(n: usize) -> f64 {
    two_pow_n(n) as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_two_pow_n() {
        assert_eq!(2, two_pow_n(1));
        assert_eq!(4, two_pow_n(2));
        assert_eq!(8, two_pow_n(3));
        assert_eq!(16, two_pow_n(4));
        assert_eq!(32, two_pow_n(5));
        assert_eq!(64, two_pow_n(6));
        assert_eq!(128, two_pow_n(7));
        assert_eq!(256, two_pow_n(8));
        assert_eq!(512, two_pow_n(9));
        assert_eq!(1024, two_pow_n(10));
        assert_eq!(2048, two_pow_n(11));
        assert_eq!(4096, two_pow_n(12));
        assert_eq!(8192, two_pow_n(13));
        assert_eq!(16384, two_pow_n(14));
        assert_eq!(32768, two_pow_n(15));
    }

    #[test]
    fn test_two_pow_n_f64() {
        assert!((2.0 - two_pow_n_f64(1)).abs() < f64::EPSILON);
        assert!((4.0 - two_pow_n_f64(2)).abs() < f64::EPSILON);
        assert!((8.0 - two_pow_n_f64(3)).abs() < f64::EPSILON);
        assert!((16.0 - two_pow_n_f64(4)).abs() < f64::EPSILON);
        assert!((32.0 - two_pow_n_f64(5)).abs() < f64::EPSILON);
        assert!((64.0 - two_pow_n_f64(6)).abs() < f64::EPSILON);
        assert!((128.0 - two_pow_n_f64(7)).abs() < f64::EPSILON);
        assert!((256.0 - two_pow_n_f64(8)).abs() < f64::EPSILON);
        assert!((512.0 - two_pow_n_f64(9)).abs() < f64::EPSILON);
        assert!((1024.0 - two_pow_n_f64(10)).abs() < f64::EPSILON);
        assert!((2048.0 - two_pow_n_f64(11)).abs() < f64::EPSILON);
        assert!((4096.0 - two_pow_n_f64(12)).abs() < f64::EPSILON);
        assert!((8192.0 - two_pow_n_f64(13)).abs() < f64::EPSILON);
        assert!((16384.0 - two_pow_n_f64(14)).abs() < f64::EPSILON);
        assert!((32768.0 - two_pow_n_f64(15)).abs() < f64::EPSILON);
    }
}
