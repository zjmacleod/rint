//! Numerical integration routines for multi-dimensional functions.
//!
//! This module provides numerical integration routines for integrating functions with
//! dimensionality between $2 \le N \le 15$. The functions can be either real-valued or complex,
//! and are integrated over an $N$-dimensional hypercube. The routines are based primarily on the
//! DCUHRE FORTRAN library (Bernsten, Espelid, Genz) \[1\], however unlike the original algorithm
//! the routines presented in currently only operate on a single function _not_ a vector of
//! functions. The module provides two classes of routine:
//!
//! - [`Adaptive`]: A $2 \le N \le 15$ dimensional adaptive routine with a similar approach to
//! the one-dimensional adaptive routines found in [`crate::quadrature`]. On each iteration of the
//! algorithm the axis along which the largest contribution to the error estimate was obtained is
//! used as the bisection axis to bisect the integration region and then calculate new estimates
//! for these newly bisected volumes. This concentrates the integration refinement to the regions
//! with highest error, rapidly reducing the numerical error of the routine. The algorithm uses
//! fully-symmetric integration rules, [`Rule`], of varying order and generality.
//! These are generated through the `Rule*::generate` constructors of specific type alias' for each
//! rule:
//!
//!     - [`Rule13`]: A 13-point fully symmetric integration rule for functions of $N=2$
//!     dimension.
//!     - [`Rule13`]: An 11-point fully symmetric integration rule for functions of $N=3$
//!     dimension.
//!     - [`Rule09N2`]: A 9-point fully symmetric integration rule for functions of $N=2$
//!     dimension.
//!     - [`Rule09`]: A 9-point fully symmetric integration rule for functions of
//!     $3 \le N \le 15$ dimension.
//!     - [`Rule07`]: A 7-point fully symmetric integration rule for functions of
//!     $2 \le N \le 15$ dimension.
//!
//! - [`Basic`]: A non-adaptive routine which applies a provided Gauss-Kronrod
//! integration [`Rule`] to a function exactly once.
//!
//! \[1\] Jarle Berntsen, Terje O. Espelid, and Alan Genz. 1991. Algorithm 698: DCUHRE: an adaptive
//! multidemensional integration routine for a vector of integrals. ACM Trans. Math. Softw. 17, 4
//! (Dec. 1991), 452–456. <https://doi.org/10.1145/210232.210234>
//!
//! # Example
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
//!```rust
//! use rint::{Limits, MultiDimensionalIntegrand, Tolerance};
//! use rint::multi::{Adaptive, Rule07};
//!
//! const N: usize = 4;
//!
//! struct F {
//!     limits: [Limits; N],
//! }
//!
//! impl F {
//!     fn new() -> Self {
//!         let limits = [
//!             Limits::new(0.0, 1.0),
//!             Limits::new(0.0, 1.0),
//!             Limits::new(0.0, 1.0),
//!             Limits::new(0.0, 2.0)
//!         ];
//!         Self { limits }
//!     }
//! }
//!
//! impl MultiDimensionalIntegrand<N> for F {
//!     type Scalar = f64;
//!     fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
//!         let [x1, x2, x3, x4] = coordinates;
//!         x3.powi(2) * x4 * (x3 * x4).exp() / (x1 + x2 + 1.0).powi(2)
//!     }
//! }
//!
//! const TARGET: f64 = 5.753_641_449_035_616e-1;
//! const TOL: f64 = 1e-2;
//!
//! let function = F::new();
//! let limits = function.limits;
//! let rule = Rule07::<N>::generate().unwrap();
//! let tolerance = Tolerance::Relative(TOL);
//!
//! let integrator = Adaptive::new(
//!     &function,
//!     &rule,
//!     limits,
//!     tolerance,
//!     10000
//! ).unwrap();
//!
//! let integral = integrator.integrate().unwrap();
//! let result = integral.result();
//! let error = integral.error();
//!
//! let actual_error = (result - TARGET).abs();
//! let requested_error = TOL * result.abs();
//!
//! assert!(actual_error < requested_error);
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
#[allow(clippy::float_cmp)]
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
        assert_eq!(2.0, two_pow_n_f64(1));
        assert_eq!(4.0, two_pow_n_f64(2));
        assert_eq!(8.0, two_pow_n_f64(3));
        assert_eq!(16.0, two_pow_n_f64(4));
        assert_eq!(32.0, two_pow_n_f64(5));
        assert_eq!(64.0, two_pow_n_f64(6));
        assert_eq!(128.0, two_pow_n_f64(7));
        assert_eq!(256.0, two_pow_n_f64(8));
        assert_eq!(512.0, two_pow_n_f64(9));
        assert_eq!(1024.0, two_pow_n_f64(10));
        assert_eq!(2048.0, two_pow_n_f64(11));
        assert_eq!(4096.0, two_pow_n_f64(12));
        assert_eq!(8192.0, two_pow_n_f64(13));
        assert_eq!(16384.0, two_pow_n_f64(14));
        assert_eq!(32768.0, two_pow_n_f64(15));
    }
}
