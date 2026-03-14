//! Numerical integration routines for one-dimensional functions.
//!
//!
//! # Overview
//!
//! This module provides a number of numerical quadrature integration routines of a function in one
//! dimension,
//! $$
//! I = \int_{b}^{a} f(x) dx,
//! $$
//! using a Gaussian numerical integration routine. The function $f(x)$ implements the
//! [`Integrand`] trait. A Gaussian numerical integration rule approximates $I$ by performing a
//! weighted sum of the function evaluated at defined points/abscissae,
//! $$
//! I = \int_{b}^{a} f(x) dx \approx I_{n} = \sum_{i = 1}^{n} W_{i} f(X_{i})
//! $$
//! where the $X_{i}$ and $W_{i}$ are the rescaled abscissae and weights,
//! $$
//! X_{i} = \frac{b + a + (a - b) x_{i}}{2} ~~~~~~~~ W_{i} = \frac{(a - b) w_{i}}{2}
//! $$
//! The routines are based primarily on the algorithms presented in the Fortran library [QUADPACK]
//! (Piessens, de Doncker-Kapenga, Ueberhuber and Kahaner) and reimplemented in the [GNU scientific
//! library] (GSL) (Gough, Alken, Gonnet, Holoborodko, Griessinger). The integrators use
//! Gauss-Kronrod integration [`Rule`]s, combining two rules of different order for efficient
//! estimation of the numerical error. The rules for an $n$-point Gauss-Kronrod rule contain $m =
//! (n - 1) / 2$ abscissae _shared_ by the Gaussian and Kronrod rules and an extended set of $n -
//! m$ Kronrod abscissae. The estimate calculated from the order $n$ rule, $I_{n}$, and the
//! estimate calculated from the order $m$ rule, $I_{m}$, are used to determine the error estimate,
//! $$
//! E = |I_{n} - I_{m}|.
//! $$
//! Thus, only $n$ total function evaluations are required to estimate both the integral and the
//! error. To use the integrators, the user implements the trait [`Integrand`] on a type.
//!
//!
//! # Available integrator routines
//!
//! The `rint` library departs from the naming conventions of the [QUADPACK] and [GSL], but
//! provides a selection of comparable implementations:
//!
//! - [`Basic`]: A non-adaptive routine which applies a provided Gauss-Kronrod
//! integration [`Rule`] to a function exactly once.
//!
//! - [`Adaptive`]: A one-dimensional adaptive routine based on the `qag.f` [QUADPACK]
//! and `qag.c` [GSL] implementations. The integration is region is bisected into subregions and an
//! initial estimate is performed. Upon each further iteration of the routine the subregion with
//! the highest numerical error estimate is bisected again and new estimates are calculated for
//! these newly bisected regions. This concentrates the integration refinement to the regions with
//! highest error, rapidly reducing the numerical error of the routine. Gauss-Kronrod integration
//! [`Rule`]s are provided of various order to use with the adaptive algorithm.
//!
//! - [`AdaptiveSingularity`]: A one-dimensional adaptive routine based on the `qags.f`
//! [QUADPACK] and `qags.c` [GSL] implementations. Adaptive routines concentrate new subintervals
//! around the region of highest error. If this region contains an integrable singularity, then the
//! adaptive routine of [`Adaptive`] may fail to obtain a suitable estimate. However,
//! this can be combined with an extrapolation proceedure such as the Wynn epsilon-algorithm to
//! extrapolate the value at these integrable singularities and provide a suitable estimate. As
//! well as handling integrable singularities, [`AdaptiveSingularity`] can be used to
//! calculate integrals with infinite or semi-infinite bounds by using the appropriate constructors.
//!
//! [QUADPACK]: <https://www.netlib.org/quadpack/>
//! [GNU Scientific Library]: <https://www.gnu.org/software/gsl/doc/html/integration.html>
//! [GSL]: <https://www.gnu.org/software/gsl/doc/html/integration.html>
//! [`AdaptiveSingularity`]: crate::quadrature::AdaptiveSingularity
//! [`Adaptive`]: crate::quadrature::Adaptive
//! [`Basic`]: crate::quadrature::Basic
//! [`Integrand`]: crate::Integrand
//!
//! # Examples
//!
//! ## [`Basic`] integrator example
//!
//! Here we present a calculation of the golden ratio $\varphi$ ([`std::f64::consts::GOLDEN_RATIO`])
//! using the integral representation,
//! $$
//! \ln \varphi = \int_{0}^{1/2} \frac{dx}{\sqrt{1 + x^{2}}}
//! $$
//! using the [`Basic`] integrator.
//!```rust
//! use rint::{Integrand, Limits};
//! use rint::quadrature::{Basic, Rule};
//!
//! const PHI: f64 = std::f64::consts::GOLDEN_RATIO;
//!
//! struct GoldenRatio;
//!
//! impl Integrand for GoldenRatio {
//!     type Scalar = f64;
//!
//!     fn evaluate(&self, x: f64) -> Self::Scalar {
//!         1.0 / (1.0 + x.powi(2)).sqrt()
//!     }
//! }
//!
//! let golden_ratio = GoldenRatio;
//! let limits = Limits::new(0.0,0.5);
//! let rule = Rule::gk15();
//! let integral = Basic::new(&golden_ratio, &rule, limits)
//!     .integrate();
//!
//! let result = integral.result();
//! let error = integral.error();
//! let abs_actual_error = (PHI.ln() - result).abs();
//! let iters = integral.iterations();
//! assert_eq!(iters, 1);
//! assert!(abs_actual_error < error);
//!```
//!
//! ## [`Adaptive`] integrator example
//!
//! Here we present a calculation of the integral,
//! $$
//! I = \int_{0}^{1} x^{\alpha} \ln \frac{1}{x} dx = \frac{1}{(1+\alpha)^{2}}
//! $$
//! for different values of $\alpha$ using the [`Adaptive`] integrator.
//!
//!```rust
//! use rint::{Integrand, Limits, Tolerance};
//! use rint::quadrature::{Adaptive, Basic, Rule};
//!
//! use std::f64::consts::*;
//!
//! struct Function1 {
//!     alpha: f64,
//! }
//!
//! impl Integrand for Function1 {
//!     type Scalar = f64;
//!     fn evaluate(&self, x: f64) -> Self::Scalar {
//!         let alpha = self.alpha;
//!         x.powf(alpha) * (1.0 / x).ln()
//!     }
//! }
//!
//! const TOL: f64 = 1.0e-12;
//!
//! let tolerance = Tolerance::Relative(TOL);
//! let limits = Limits::new(0.0, 1.0);
//! let rule = Rule::gk31();
//! let max_iterations = 1000;
//!
//! let alpha_values = [2.6, PI, EULER_GAMMA, LOG10_E, 100.0, PI.powi(3)];
//!
//! for alpha in alpha_values {
//!     let function = Function1 { alpha };
//!
//!     let integral = Adaptive::new(&function, &rule, limits, tolerance, max_iterations)
//!                     .unwrap()
//!                     .integrate()
//!                     .unwrap();
//!
//!     let target = 1.0 / (1.0 + alpha).powi(2);
//!     let result = integral.result();
//!     let error = integral.error();
//!     let abs_actual_error = (result - target).abs();
//!     let tol = TOL * result.abs();
//!
//!     assert!(abs_actual_error < error);
//!     assert!(error < tol);
//! }
//!```
//!
//!
//! ## [`AdaptiveSingularity`] integrator example.
//!
//! Here we present a calculation of [Catalan's constant] $G$ using the integral representation:
//! $$
//! G = \frac{\pi}{2} \int_{1}^{+\infty} \frac{(x^{4} - 6 x^{2} + 1) \ln \ln x}{(1+x^{2})^{3}} d x
//! $$
//! which has a semi-infinite integration region $x \in (1,+\infty)$. We use
//! [`AdaptiveSingularity::semi_infinite_upper`] to approximate $G$.
//!```rust
//! use rint::{Integrand, Limits, Tolerance};
//! use rint::quadrature::AdaptiveSingularity;
//!
//! const G: f64 = 0.915_965_594_177_219_015_054_603_514_932_384_110_774;
//! const TOL: f64 = 1.0e-12;
//!
//! struct Catalan;
//!
//! impl Integrand for Catalan {
//!     type Scalar = f64;
//!
//!     fn evaluate(&self, x: f64) -> Self::Scalar {
//!         let FRAC_PI_2 = std::f64::consts::FRAC_PI_2;
//!         let polynomial = x.powi(4) - 6.0 * x.powi(2) + 1.0;
//!         let denominator = (1.0 + x.powi(2)).powi(3);
//!         let lnlnx = x.ln().ln();
//!         
//!         FRAC_PI_2 * polynomial * lnlnx / denominator
//!     }
//! }
//!
//! let catalan = Catalan;
//! let lower = 1.0;
//! let tolerance = Tolerance::Relative(TOL);
//! let integral = AdaptiveSingularity::semi_infinite_upper(catalan, lower, tolerance, 1000)
//!     .unwrap()
//!     .integrate()
//!     .unwrap();
//!
//! let result = integral.result();
//! let error = integral.error();
//! let abs_actual_error = (G - result).abs();
//! let tol = TOL * result.abs();
//! let iters = integral.iterations();
//! assert_eq!(iters, 32);
//! assert!(abs_actual_error < error);
//! assert!(error < tol);
//!```
//!
//!
//! ## Complex valued [`AdaptiveSingularity`] integrator example.
//!
//! The `rint` library also supports the integration of complex valued integrands. Here we
//! demonstrate the evaluation of the Fourier transform,
//! $$
//! F(\xi) = \int_{-\infty}^{+\infty} f(x) e^{- i \xi x} dx
//! $$
//! where $f(x) = \sin \omega x$ when $|x| < 1$ and $f(x) = 0$ otherwise. This has an exact value,
//! $$
//! F(\xi) = i \frac{\sin(\omega + \xi)}{\omega + \xi} - i \frac{\sin(\omega - \xi)}{\omega - \xi}
//! $$
//! We use [`AdaptiveSingularity::infinite`] to approximate $F(\xi)$ for some choices of
//! $(\omega,\xi)$.
//!```rust
//! use rint::{Integrand, Tolerance};
//! use rint::quadrature::AdaptiveSingularity;
//! use num_complex::{Complex, ComplexFloat};
//! use std::f64::consts::*;
//!
//! const TOL: f64 = 1.0e-12;
//!
//! struct Fourier {
//!     omega: f64,
//!     xi: f64,
//! }
//!
//! impl Fourier {
//!     fn new(omega: f64, xi: f64) -> Self {
//!         Self { omega, xi }
//!     }
//!
//!     fn exact(&self) -> Complex<f64> {
//!         let omega = self.omega;
//!         let xi = self.xi;
//!         let plus = (omega + xi);
//!         let minus = (omega - xi);
//!
//!         Complex::I * ((plus.sin() / plus) - (minus.sin() / minus))
//!     }
//! }
//!
//! impl Integrand for Fourier {
//!     type Scalar = Complex<f64>;
//!
//!     fn evaluate(&self, x: f64) -> Self::Scalar {
//!         if x.abs() >= 1.0 {
//!             Complex::ZERO
//!         } else {
//!             let omega = self.omega;
//!             let xi = self.xi;
//!             
//!             let f = (omega * x).sin();
//!             let e = Complex::cis(-xi * x);
//!
//!             f * e
//!         }
//!     }
//! }
//!
//! let omegavals = [0.0, EULER_GAMMA, LOG10_E, PI.powi(3), GOLDEN_RATIO];
//! let xivals = omegavals;
//! let mut fouriers = Vec::new();
//!
//! for xi in xivals {
//!     for omega in omegavals {
//!         if ((xi - omega).abs() > f64::EPSILON) {
//!             fouriers.push(Fourier::new(omega,xi));
//!         }
//!     }
//! }
//!
//! let tolerance = Tolerance::Relative(TOL);
//!
//! for fourier in fouriers {
//!     let integral = AdaptiveSingularity::infinite(&fourier, tolerance, 1000)
//!         .unwrap()
//!         .integrate()
//!         .unwrap();
//!
//!     let exact = fourier.exact();
//!
//!     let result = integral.result();
//!     let error = integral.error();
//!     let abs_actual_error = (exact - result).abs();
//!     let tol = TOL * result.abs();
//!     let iters = integral.iterations();
//!     assert!(abs_actual_error <= error);
//!     assert!(error <= tol);
//! }
//!```
//!
//! [Catalan's constant]: <https://en.wikipedia.org/wiki/Catalan%27s_constant>
//! [Euler's constant]: <https://en.wikipedia.org/wiki/Euler%27s_constant>
//! [`Basic`]: crate::quadrature::Basic
//! [`Adaptive`]: crate::quadrature::Adaptive
//! [`AdaptiveSingularity`]: crate::quadrature::AdaptiveSingularity
//! [`Tolerance`]: crate::Tolerance
//! [`Error`]: crate::Error
//!
mod adaptive;
mod basic;
mod integrator;
mod region;
pub(crate) mod rule;
mod singularity;

#[cfg(test)]
mod tests;

pub(crate) use integrator::Integrator;
pub(crate) use region::Region;

pub use adaptive::Adaptive;
pub use basic::Basic;
pub use singularity::AdaptiveSingularity;

pub use rule::Rule;

use crate::IntegralEstimate;

pub(crate) fn rescale_error(error: f64, result_abs: f64, result_asc: f64) -> f64 {
    let mut error = error;

    if result_asc != 0.0 && error != 0.0 {
        let scale = (200.0 * error / result_asc).powf(1.5);

        if scale < 1.0 {
            error = result_asc * scale;
        } else {
            error = result_asc;
        }
    }

    if result_abs > f64::MIN_POSITIVE / (50.0 * f64::EPSILON) {
        let min_error = 50.0 * f64::EPSILON * result_abs;

        if min_error > error {
            error = min_error;
        }
    }

    error
}
