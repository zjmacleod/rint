//! Numerical integration routines written in Rust.
//!
//!
//!
//! # Overview
//!
//! This library contains numerical integration routines for both functions of one dimension (see
//! [`quadrature`]) and functions with multiple dimensions up to $N = 15$ (see [`multi`]). The
//! basic principle of the library is to expose all of the routines through the trait system. Each
//! of the one- and multi-dimensional integrators take as a parameter a type implementing the
//! trait [`Integrand`].
//!
//! The integration routines attempt to make approximations to integrals such as the
//! one-dimensional integral,
//! $$
//! I = \int_{b}^{a} f(x) dx
//! $$
//! or the $N$-dimensional
//! $$
//! I = \int_{\Sigma_{N}} f(\mathbf{x}) d\mathbf{x}
//! $$
//! where $\mathbf{x} = (x_{1}, x_{2}, \dots, x_{N})$ and $\Sigma_{N}$ is an $N$-dimensional
//! hypercube. The functions $f(x)$ and $f(\mathbf{x})$ can be real valued, with return type
//! [`f64`] _or_ complex valued with return type [`Complex<f64>`]. The numerical integration
//! routines approximate the integral of a function by performing a weighted sum of the function
//! evaluated at defined points/abscissae. For example, in the one-dimensional case,
//! $$
//! I = \int_{b}^{a} f(x) dx \approx \sum_{i = 1}^{n} W_{i} f(X_{i}) = I_{n}
//! $$
//! where the $X_{i}$ and $W_{i}$ are the rescaled abscissae and weights,
//! $$
//! X_{i} = \frac{b + a + (a - b) x_{i}}{2} ~~~~~~~~ W_{i} = \frac{(a - b) w_{i}}{2}
//! $$
//! Integration rules have a polynomial order. Rules of higher polynomial order use more
//! points/abscissae and weights for evaluating the sum.
//!
//! The library contains both non-adaptive and adaptive integration routines. In the case of an
//! adaptive routine, the user supplies an error [`Tolerance`] which acts as a goal for the
//! numerical error estimate. The numerical integration is considered successful when the numerical
//! error is less than the user supplied tolerance. On success, the output of the numerical
//! integration is an [`IntegralEstimate`].
//!
//!
//! # [`Integrand`] trait
//!
//! The primary entry point for the library is the [`Integrand`] trait. A type implementing the
//! [`Integrand`] trait represents a real or complex valued function which is to be integrated. The
//! trait requires the definition of two associated types, [`Integrand::Point`] and
//! [`Integrand::Scalar`], and an implementation of the method [`Integrand::evaluate`].
//!
//! - [`Integrand::Point`]: This associated type defines the point at which the function is to be
//! evaluated, and determines the types of numerical integrators which are available to the user to
//! integrate the function. Integrators are provided for univariate functions $f(x)$ through the
//! associated type `Point=f64`, while integrators for multivariate functions $f(\mathbf{x})$ are
//! provided through the associated type `Point=[f64;N]` where $N$ is the dimensionality of the
//! point $\mathbf{x}=(x_{1},\dots,x_{N})$ which is limited to $2 \le N \le 15$.
//!
//! - [`Integrand::Scalar`]: This is the output type of the function to be integrated. A _real_
//! valued function should have the output type `Scalar=`[`f64`], while a _complex_ valued function
//! should have output type `Scalar=`[`Complex<f64>`].
//!
//! - [`Integrand::evaluate`]: The trait requires an implementation of an `evaluate` method, which
//! defines how the function takes the input [`Integrand::Point`] and turns this into the output
//! type [`Integrand::Scalar`]. In other words, this method tells the integrators how to evaluate
//! the function at a point, allowing the integration to be done.
//!
//! As an example, consider probability density function of a normal distribution,
//! $$
//! f(x) = \frac{1}{\sqrt{2 \pi \sigma^{2}}} e^{- \frac{(x - \mu)^{2}}{2 \sigma^{2}}}
//! $$
//! which has mean $\mu$ and standard deviation $\sigma$. To integrate this function, one first
//! implements the [`Integrand`] trait,
//!```rust
//! use rint::Integrand;
//!
//! struct NormalDist {
//!     mean: f64,
//!     standard_dev: f64,
//! }
//!
//! impl Integrand for NormalDist {
//!     type Point = f64;
//!     type Scalar = f64;
//!     fn evaluate(&self, x: &Self::Point) -> Self::Scalar {
//!         let mu = self.mean;
//!         let sigma = self.standard_dev;
//!
//!         let prefactor = 1.0 / (2.0 * std::f64::consts::PI * sigma.powi(2));
//!         let exponent = - (x - mu).powi(2) / (2.0 * sigma.powi(2));
//!
//!         prefactor * exponent.exp()
//!     }
//! }
//!```
//! The type `NormalDist` can then be passed as a parameter to a one-dimensional integration
//! routine in the module [`quadrature`] and integrated over a given (possibly infinite) region.
//! Since [`Integrand::evaluate`] does not return a [`Result`] it is the users responsibility to
//! ensure that the evaluation implementation is correct.
//!
//!
//! # Modules
//!
//! ## [`quadrature`]
//!
//! The [`quadrature`] module provides a number of numerical quadrature integration routines of a
//! function in one dimension. The routines are based primarily on the algorithms presented in the
//! Fortran library [QUADPACK] (Piessens, de Doncker-Kapenga, Ueberhuber and Kahaner) and
//! reimplemented in the [GNU scientific library] (GSL) (Gough, Alken, Gonnet, Holoborodko,
//! Griessinger). The integrators use Gauss-Kronrod integration rules, combining two rules of
//! different order for efficient estimation of the numerical error. The rules for an $n$-point
//! Gauss-Kronrod rule contain $m = (n - 1) / 2$ abscissae _shared_ by the Gaussian and Kronrod
//! rules and an extended set of $n - m$ Kronrod abscissae. Thus, only $n$ total function
//! evaluations are required for both the integral and error estimates.
//!
//! The `rint` library departs from the naming conventions of the [QUADPACK] and [GSL], but
//! provides a selection of comparable implementations:
//!
//! - [`quadrature::Basic`]: A non-adaptive routine which applies a provided Gauss-Kronrod
//! integration [`quadrature::Rule`] to a function exactly once.
//!
//! - [`quadrature::Adaptive`]: A one-dimensional adaptive routine based on the `qag.f` [QUADPACK]
//! and `qag.c` [GSL] implementations. The integration is region is bisected into subregions and an
//! initial estimate is performed. Upon each further iteration of the routine the subregion with
//! the highest numerical error estimate is bisected again and new estimates are calculated for
//! these newly bisected regions. This concentrates the integration refinement to the regions with
//! highest error, rapidly reducing the numerical error of the routine. Gauss-Kronrod integration
//! [`quadrature::Rule`]s are provided of various order to use with the adaptive algorithm.
//!
//! - [`quadrature::AdaptiveSingularity`]: A one-dimensional adaptive routine based on the `qags.f`
//! [QUADPACK] and `qags.c` [GSL] implementations. Adaptive routines concentrate new subintervals
//! around the region of highest error. If this region contains an integrable singularity, then the
//! adaptive routine of [`quadrature::Adaptive`] may fail to obtain a suitable estimate. However,
//! this can be combined with an extrapolation proceedure such as the Wynn epsilon-algorithm to
//! extrapolate the value at these integrable singularities and provide a suitable estimate. As
//! well as handling integrable singularities, [`quadrature::AdaptiveSingularity`] can be used to
//! calculate integrals with infinite or semi-infinite bounds by using the appropriate constructors.
//!
//! [QUADPACK]: <https://www.netlib.org/quadpack/>
//! [GNU Scientific Library]: <https://www.gnu.org/software/gsl/doc/html/integration.html>
//! [GSL]: <https://www.gnu.org/software/gsl/doc/html/integration.html>
//!
//! [`quadrature::Rule`]: crate::quadrature::Rule
//!
//! ## [`multi`]
//!
//! The [`multi`] module provides numerical integration routines for integrating functions with
//! dimensionality between $2 \le N \le 15$. The functions can be either real-valued or complex,
//! and are integrated over an $N$-dimensional hypercube. The routines are based primarily on the
//! [DCUHRE] FORTRAN library (Bernsten, Espelid, Genz), which use sets of fully-symmetric
//! integration rules to obtain integral and error estimates. Unlike the original algorithm the
//! routines presented in [`multi`] currently only operate on a single function _not_ a vector of
//! functions. The module provides two classes of routine:
//!
//! - [`multi::Basic`]: A non-adaptive routine which applies a fully-symmetric integration rule
//! [`multi::Rule`] to a multi-dimensional function exactly once.
//!
//! - [`multi::Adaptive`]: A $2 \le N \le 15$ dimensional adaptive routine with a similar approach
//! to the one-dimensional adaptive routines found in [`quadrature`]. On each iteration of the
//! algorithm the axis along which the largest contribution to the error estimate was obtained is
//! used as the bisection axis to bisect the integration region and then calculate new estimates
//! for these newly bisected volumes. This concentrates the integration refinement to the regions
//! with highest error, rapidly reducing the numerical error of the routine. The algorithm uses
//! fully-symmetric integration rules, [`crate::multi::Rule`], of varying order and generality.
//! These are generated through the `Rule*::generate` constructors.
//!
//! [`multi::Rule`]: crate::multi::Rule
//! [DCUHRE]: <https://dl.acm.org/doi/10.1145/210232.210234>
//!
//! # One-dimensional example
//!
//! The following example integrates the function
//! $$
//! f(x) = \frac{\log(x)}{(1 + 100 x^{2})}
//! $$
//! over the
//! semi-infinite interval $0 < x < \infty$ using an adaptive integration routine with singularity
//! detection (see [`quadrature::AdaptiveSingularity`]). Adapted from the [QUADPACK] & [GSL] numerical
//! integration test suites.
//!
//! [`quadrature::AdaptiveSingularity`]: crate::quadrature::AdaptiveSingularity
//!
//!```rust
//! use rint::{Limits, Integrand, Tolerance};
//! use rint::quadrature::AdaptiveSingularity;
//!
//! struct F;
//!
//! impl Integrand for F {
//!     type Point = f64;
//!     type Scalar = f64;
//!     fn evaluate(&self, x: &Self::Point) -> Self::Scalar {
//!         x.ln() / (1.0 + 100.0 * x.powi(2))
//!     }
//! }
//!
//! # use std::error::Error;
//! # fn main() -> Result<(), Box<dyn Error>> {
//! let function = F;
//! let lower = 0.0;
//! let tolerance = Tolerance::Relative(1.0e-3);
//!
//! let integrator = AdaptiveSingularity::semi_infinite_upper(&function, lower, tolerance, 1000)?;
//!
//! let exp_result = -3.616_892_186_127_022_568E-01;
//! let exp_error = 3.016_716_913_328_831_851E-06;
//!
//! let integral = integrator.integrate()?;
//! let result = integral.result();
//! let error = integral.error();
//!
//! let tol = 1.0e-9;
//! assert!((exp_result - result).abs() / exp_result.abs() < tol);
//! assert!((exp_error - error).abs() / exp_error.abs() < tol);
//! # Ok(())
//! # }
//!```
//!
//! # Multi-dimensional example
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
//! use rint::{Limits, Integrand, Tolerance};
//! use rint::multi::{Adaptive, Rule07};
//!
//! const N: usize = 4;
//!
//! struct F;
//!
//! impl Integrand for F {
//!     type Point = [f64; N];
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
//!     Limits::new(0.0, 1.0)?,
//!     Limits::new(0.0, 1.0)?,
//!     Limits::new(0.0, 1.0)?,
//!     Limits::new(0.0, 2.0)?
//! ];
//! let rule = Rule07::<N>::generate()?;
//! let tolerance = Tolerance::Relative(TOL);
//!
//! let integrator = Adaptive::new(&function, &rule, limits, tolerance, 10000)?;
//!
//! let integral = integrator.integrate()?;
//! let result = integral.result();
//! let error = integral.error();
//!
//! let actual_error = (result - TARGET).abs();
//! let requested_error = TOL * result.abs();
//!
//! assert!(actual_error < error);
//! assert!(error < requested_error);
//! # Ok(())
//! # }
//!```
#![deny(clippy::pedantic)]
#![allow(
    clippy::excessive_precision,
    clippy::doc_lazy_continuation,
    clippy::cast_precision_loss,
    clippy::if_not_else
)]
use num_complex::{Complex, ComplexFloat};
use num_traits::Zero;

use std::{error, fmt, ops};

mod estimate;
mod limits;
pub mod multi;
pub mod quadrature;

pub use estimate::IntegralEstimate;
pub use limits::Limits;

/// A real or complex-valued function to be integrated.
///
/// The primary entry point for the library is the [`Integrand`] trait. A type implementing the
/// [`Integrand`] trait represents a real or complex valued function which is to be integrated. The
/// trait requires the definition of two associated types, [`Integrand::Point`] and
/// [`Integrand::Scalar`], and an implementation of the method [`Integrand::evaluate`].
///
/// - [`Integrand::Point`]: This associated type defines the point at which the function is to be
/// evaluated, and determines the types of numerical integrators which are available to the user to
/// integrate the function. Integrators are provided for univariate functions $f(x)$ through the
/// associated type `Point=f64`, while integrators for multivariate functions $f(\mathbf{x})$ are
/// provided through the associated type `Point=[f64;N]` where $N$ is the dimensionality of the
/// point $\mathbf{x}=(x_{1},\dots,x_{N})$ which is limited to $2 \le N \le 15$.
///
/// - [`Integrand::Scalar`]: This is the output type of the function to be integrated. A _real_
/// valued function should have the output type `Scalar=`[`f64`], while a _complex_ valued function
/// should have output type `Scalar=`[`Complex<f64>`].
///
/// - [`Integrand::evaluate`]: The trait requires an implementation of an `evaluate` method, which
/// defines how the function takes the input [`Integrand::Point`] and turns this into the output
/// type [`Integrand::Scalar`]. In other words, this method tells the integrators how to evaluate
/// the function at a point, allowing the integration to be done.
///
/// # Examples
///
/// ## One-dimensional example
///
/// For example, consider probability density function of a normal distribution,
/// $$
/// f(x) = \frac{1}{\sqrt{2 \pi \sigma^{2}}} e^{- \frac{(x - \mu)^{2}}{2 \sigma^{2}}}
/// $$
/// which has mean $\mu$ and standard deviation $\sigma$. To integrate this function, one first
/// implements the [`Integrand`] trait,
///```rust
/// use rint::Integrand;
///
/// struct NormalDist {
///     mean: f64,
///     standard_dev: f64,
/// }
///
/// impl Integrand for NormalDist {
///     type Point = f64;
///     type Scalar = f64;
///     fn evaluate(&self, x: &Self::Point) -> Self::Scalar {
///         let mu = self.mean;
///         let sigma = self.standard_dev;
///
///         let prefactor = 1.0 / (2.0 * std::f64::consts::PI * sigma.powi(2));
///         let exponent = - (x - mu).powi(2) / (2.0 * sigma.powi(2));
///
///         prefactor * exponent.exp()
///     }
/// }
///```
/// The type `NormalDist` can then be used as the `function` parameter in one of the numerical
/// integration routines provided by the [`quadrature`] module.
///
/// ## Multi-dimensional example
///
/// For example, consider a nonlinear 4-dimensional function $f(\mathbf{x})$,
/// $$
/// f(\mathbf{x}) = a \frac{x_{3}^{2} x_{4} e^{x_{3} x_{4}}}{(1 + x_{1} + x_{2})^{2}}
/// $$
/// where $\mathbf{x} = (x_{1}, x_{2}, x_{3}, x_{4})$, which has some amplitude $a$. To integrate
/// this function, one first implements the [`Integrand`] trait,
///```rust
/// use rint::Integrand;
///
/// const N: usize = 4;
///
/// struct F {
///     amplitude: f64,
/// }
///
/// impl Integrand for F {
///     type Point = [f64; N];
///     type Scalar = f64;
///     fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
///         let [x1, x2, x3, x4] = coordinates;
///         let a = self.amplitude;
///         a * x3.powi(2) * x4 * (x3 * x4).exp() / (x1 + x2 + 1.0).powi(2)
///     }
/// }
///```
/// The type `F` can then be used as the `function` parameter in one of the numerical integration
/// routines provided by the [`multi`] module.
pub trait Integrand {
    /// The input point at which the function is evaluated.
    ///
    /// This associated type defines the point at which the function is to be evaluated, and
    /// determines the types of numerical integrators which are available to the user to integrate
    /// the function. Integrators are provided for univariate functions $f(x)$ through the
    /// associated type `Point=f64`, while integrators for multivariate functions $f(\mathbf{x})$
    /// are provided through the associated type `Point=[f64;N]` where $N$ is the dimensionality of
    /// the point $\mathbf{x}=(x_{1},\dots,x_{N})$ which is limited to $2 \le N \le 15$.
    type Point;

    /// The type that the function evaluates to.
    ///
    /// This is the output type of the function to be integrated. A _real_ valued function should
    /// have the output type `Scalar=`[`f64`], while a _complex_ valued function should have output
    /// type `Scalar=`[`Complex<f64>`].
    type Scalar: ScalarF64;

    /// Evaluate the function at the real point `x`
    ///
    /// The trait requires an implementation of an `evaluate` method, which defines how the
    /// function takes the input [`Integrand::Point`] and turns this into the output type
    /// [`Integrand::Scalar`]. In other words, this method tells the integrators how to evaluate
    /// the function at a point, allowing the integration to be done.
    fn evaluate(&self, x: &Self::Point) -> Self::Scalar;
}

impl<I: Integrand> Integrand for &I {
    type Point = I::Point;
    type Scalar = I::Scalar;

    fn evaluate(&self, x: &Self::Point) -> Self::Scalar {
        I::evaluate(self, x)
    }
}

/// A numerical scalar.
///
/// The [`Integrand`] trait is implemented for both real- and complex-valued functions of a _real_
/// variable. The *sealed* trait [`ScalarF64`] is implemented for both [`f64`] and [`Complex<f64>`].
pub trait ScalarF64:
    PartialEq
    + ComplexFloat<Real = f64>
    + Zero
    + Copy
    + ops::Mul<f64, Output = Self>
    + ops::Div<f64, Output = Self>
    + for<'a> ops::Mul<&'a f64, Output = Self>
    + for<'a> ops::Div<&'a f64, Output = Self>
    + ops::AddAssign<Self>
    + ops::SubAssign<Self>
    + ops::Add<Self>
    + ops::Sub<Self>
    + for<'a> ops::AddAssign<&'a Self>
    + for<'a> ops::SubAssign<&'a Self>
    + for<'a> ops::Add<&'a Self, Output = Self>
    + for<'a> ops::Sub<&'a Self, Output = Self>
    + fmt::Display
    + fmt::Debug
    + sealed::Sealed
{
    /// Numerical integration routines require the scalar type to have a sensible definition of a
    /// maximum value. For [`f64`] this is just [`f64::MAX`]. For [`Complex<f64>`] this is chosen
    /// as `Complex(f64::MAX / 2f64.sqrt(), f64::MAX / 2f64.sqrt())`.
    fn max_value() -> Self;
}

impl ScalarF64 for f64 {
    fn max_value() -> Self {
        f64::MAX
    }
}

impl ScalarF64 for Complex<f64> {
    fn max_value() -> Self {
        Complex::new(f64::MAX / 2f64.sqrt(), f64::MAX / 2f64.sqrt())
    }
}

pub(crate) mod sealed {
    use num_complex::Complex;

    pub trait Sealed {}

    impl Sealed for f64 {}
    impl Sealed for Complex<f64> {}
}

/// User selected tolerance type.
///
/// The adaptive routines will return the first approximation, `result`, to the integral which has an
/// absolute `error` smaller than the tolerance set by the choice of [`Tolerance`], where
/// * [`Tolerance::Absolute(abserr)`] specifies an absolute error and returns final
/// [`IntegralEstimate`] when `error <= abserr`,
/// * [`Tolerance::Relative(relerr)`] specifies a relative error and returns final
/// [`IntegralEstimate`] when `error <= relerr * abs(result)`,  
/// * [`Tolerance::Either{ abserr, relerr }`] to return a result as soon as _either_ the relative or
/// absolute error bound has been satisfied.
///
/// [`Tolerance::Absolute(abserr)`]: crate::Tolerance#variant.Absolute
/// [`Tolerance::Relative(relerr)`]: crate::Tolerance#variant.Relative
/// [`Tolerance::Either{ abserr, relerr }`]: crate::Tolerance#variant.Either
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub enum Tolerance {
    /// Specify absolute error and returns final [`IntegralEstimate`] when `error<=abserr`.
    Absolute(f64),
    /// Specify relative error and return final [`IntegralEstimate`] when `error<=relerr*abs(result)`.
    Relative(f64),
    /// Specify _both_ absolute and relative error and return final [`IntegralEstimate`] as soon as
    /// _either_ the relative or absolute error bound has been satisfied.
    Either { absolute: f64, relative: f64 },
}

impl Tolerance {
    #[must_use]
    pub(crate) fn tolerance<T: ScalarF64>(&self, integral_value: &T) -> f64 {
        match *self {
            Tolerance::Absolute(v) => v,
            Tolerance::Relative(v) => v * integral_value.abs(),
            Tolerance::Either { absolute, relative } => {
                f64::max(absolute, relative * integral_value.abs())
            }
        }
    }

    pub(crate) fn check(&self) -> Result<(), InitialisationError> {
        match self {
            Tolerance::Absolute(v) => {
                if *v <= 0.0 {
                    let kind = InitialisationErrorKind::AbsoluteBoundNegativeOrZero(*v);
                    return Err(InitialisationError::new(kind));
                }
            }
            Tolerance::Relative(v) => {
                if *v < 50.0 * f64::EPSILON {
                    let kind = InitialisationErrorKind::RelativeBoundTooSmall(*v);
                    return Err(InitialisationError::new(kind));
                }
            }
            Tolerance::Either { absolute, relative } => {
                if *absolute <= 0.0 && *relative < 50.0 * f64::EPSILON {
                    let kind = InitialisationErrorKind::InvalidTolerance {
                        absolute: *absolute,
                        relative: *relative,
                    };
                    return Err(InitialisationError::new(kind));
                }
            }
        }

        Ok(())
    }
}

/// General error enum for the `rint` crate.
///
/// Errors can occur in two ways in the `rint` crate, either on initialisation/setup of the
/// integrators _before_ numerical integration has been attempted *or* as a result of
/// encountering an error _during_ the running of the numerical integration routine. As a result,
/// the [`Error`] enum provides two variants [`Error::Initialisation`] to notify of errors on
/// initialisation and [`Error::Integration`] to notify of errors during integration. Each variant
/// holds a struct giving more information about the error that occurred---see
/// [`InitialisationError`] and [`IntegrationError`] for more details.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub enum Error<T: ScalarF64> {
    Initialisation(InitialisationError),
    Integration(IntegrationError<T>),
}

impl<T: ScalarF64> From<InitialisationError> for Error<T> {
    fn from(other: InitialisationError) -> Self {
        Error::Initialisation(other)
    }
}

impl<T: ScalarF64> From<IntegrationError<T>> for Error<T> {
    fn from(other: IntegrationError<T>) -> Self {
        Error::Integration(other)
    }
}

impl<T: ScalarF64> error::Error for Error<T> {}

impl<T: ScalarF64> fmt::Display for Error<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Error::Initialisation(err) => err.fmt(f),
            Error::Integration(err) => err.fmt(f),
        }
    }
}

/// Error type for errors that occur during initialisation/setup of an integrator.
///
/// Errors can occur in two ways in the `rint` crate, either on initialisation/setup of the
/// integrators _before_ numerical integration has been attempted *or* as a result of
/// encountering an error _during_ the running of the numerical integration routine.
/// [`InitialisationError`] provides further details about the error that occurred during the
/// initialisation/setup of the integrators. The only reason that an initialisation error can occur
/// is when there is bad user input when generating either the [`Tolerance`], such as a negative
/// absolute tolerance value or a relative tolerance value which is too close to [`f64::EPSILON`],
/// or when an invalid dimensionality $N$ is used to generate a multi-dimensional integrator and/or
/// rule. The kind of error [`InitialisationErrorKind`] which occurred during initialisation is
/// obtained through the [`InitialisationError::kind`] method.
///
/// See also [`Error`] and [`IntegrationError`].
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct InitialisationError {
    kind: InitialisationErrorKind,
}

impl InitialisationError {
    pub(crate) const fn new(kind: InitialisationErrorKind) -> Self {
        Self { kind }
    }

    /// Return the initialisation error kind.
    #[must_use]
    pub const fn kind(&self) -> InitialisationErrorKind {
        self.kind
    }
}

/// The kind of [`InitialisationError`] which occurred.
///
/// The only reason that an initialisation error can occur is when there is bad user input when
/// generating either the [`Tolerance`], such as a negative absolute tolerance value or a relative
/// tolerance value which is too close to [`f64::EPSILON`], or when an invalid dimensionality $N$
/// is used to generate a multi-dimensional integrator and/or rule.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub enum InitialisationErrorKind {
    /// The absolute tolerance bound `tol` requested for the adaptive integration routines must
    /// satisfy `tol > 0.0`.
    AbsoluteBoundNegativeOrZero(f64),

    /// The relative tolerance bound `tol` requested for the adaptive integration routines must
    /// satisfy `tol > 50.0 * f64::EPSILON`.
    RelativeBoundTooSmall(f64),

    /// The user provided values for `upper` and `lower` when generating a set of [`Limits`] using
    /// [`Limits::new`] did not satisfy `(upper - lower).abs() < f64::MAX` and `(upper + lower).abs()
    /// < f64::MAX` . The integration range should be rescaled to satisfy this constraint.
    IntegrationRangeTooLarge { lower: f64, upper: f64 },

    /// An invalid tolerance was requested for the adaptive integration routine.
    /// The absolute tolerance bound `abs_tol` must satisfy `abs_tol > 0.0`.
    /// The relative tolerance bound `rel_tol` must satisfy satisfy `rel_tol > 50.0 * f64::EPSILON`.
    InvalidTolerance { absolute: f64, relative: f64 },

    /// An invalid integration dimensionality `N` was used in a multi-dimensional integration.
    /// This library only provides multi-dimensional integration routines suitable for dimensions
    /// $2 <= N <= 15$.
    InvalidDimension(usize),

    /// An invalid integration dimensionality $N$ was used to construct the 7-point $N$-dimensional
    /// integration rule [`Rule07`].
    /// This rule is only suitable for dimensions $2 \le N \le 15$.
    ///
    /// [`Rule07`]: crate::multi::Rule07
    InvalidDimensionForRule07(usize),

    /// An invalid integration dimensionality $N$ was used to construct the 9-point $N$-dimensional
    /// integration rule [`Rule09`].
    /// This rule is only suitable for dimensions $3 \le N \le 15$.
    ///
    /// [`Rule09`]: crate::multi::Rule09
    InvalidDimensionForRule09(usize),
}

impl fmt::Display for InitialisationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.kind() {
            InitialisationErrorKind::AbsoluteBoundNegativeOrZero(tol) => {
                write!(f, "Invalid tolerance requested:\n`Absolute({tol})`\nThe absolute tolerance bound `tol` requested for the adaptive integration routines must satisfy `tol > 0.0`.")
            }

            InitialisationErrorKind::RelativeBoundTooSmall(tol) => {
                write!(f, "Invalid tolerance requested:\n`Relative({tol})`\nThe relative tolerance bound `tol` requested for the adaptive integration routines must satisfy `tol > 50.0 * f64::EPSILON`.")
            }

            InitialisationErrorKind::InvalidTolerance { absolute, relative } => {
                write!(f, "Invalid tolerance requested:\n`Either '{{' abs_tol: {absolute}, rel_tol: {relative} '}}'`\nThe absolute tolerance bound `abs_tol` requested for the adaptive integration routines must satisfy `abs_tol > 0.0`. The relative tolerance bound `rel_tol` requested for the adaptive integration routines must satisfy `rel_tol > 50.0 * f64::EPSILON`.")
            }

            InitialisationErrorKind::IntegrationRangeTooLarge { lower, upper } => {
                write!(f, "Invalid integration limits requested:\nWhen generating a new set of integration `Limits` using `Limits::new(lower, upper)?` the bounds must satisfy `(upper - lower).abs() < f64::MAX` and  `(upper + lower).abs() < f64::MAX`. Please rescale the integration bounds accordingly. Provided bounds were: `lower = {lower}` and `upper = {upper}`.")
            }

            InitialisationErrorKind::InvalidDimension(ndim) => {
                write!(f, "An invalid integration dimensionality `N = {ndim}` was used in a multi-dimensional integration.  This library only provides multi-dimensional integration routines suitable for dimensions `2 <= N <= 15`.")
            }

            InitialisationErrorKind::InvalidDimensionForRule07(ndim) => {
                write!(f, "An invalid integration dimensionality `N = {ndim}` was used to construct the 7-point multi-dimensional integration rule [`Rule07`]. This rule is only suitable for dimensions `2 <= N <= 15`.  ")
            }

            InitialisationErrorKind::InvalidDimensionForRule09(ndim) => {
                write!(f, "An invalid integration dimensionality `N = {ndim}` was used to construct the 9-point multi-dimensional integration rule [`Rule09`]. This rule is only suitable for dimensions `3 <= N <= 15`.  ")
            }
        }
    }
}

impl error::Error for InitialisationError {}

/// Error type for errors that occur during integration of the user supplied function.
///
/// Errors can occur in two ways in the `rint` crate, either on initialisation/setup of the
/// integrators _before_ numerical integration has been attempted *or* as a result of
/// encountering an error _during_ the running of the numerical integration routine.
/// [`IntegrationError`] provides further details about the error that occurred during the
/// integration of the user supplied function, specifically the reason for the error
/// [`IntegrationErrorKind`] (accessed through the [`IntegrationError::kind`] method) and the
/// [`IntegralEstimate`] which was calculated up to the point that the error occurred (accessed
/// through the [`IntegrationError::estimate`] method). Errors typically occur during integration
/// due to, for example, difficult integration regions, non-integrable singularities, the maximum
/// number of iterations being reached, etc. See [`IntegrationErrorKind`] for more details.
///
/// See also [`Error`] and [`InitialisationError`].
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct IntegrationError<T: ScalarF64> {
    estimate: IntegralEstimate<T>,
    kind: IntegrationErrorKind,
}

impl<T: ScalarF64> IntegrationError<T> {
    pub(crate) const fn new(estimate: IntegralEstimate<T>, kind: IntegrationErrorKind) -> Self {
        Self { estimate, kind }
    }

    /// Return the error [`IntegrationErrorKind`] which was encountered.
    pub const fn kind(&self) -> IntegrationErrorKind {
        self.kind
    }

    /// Return a reference to the best [`IntegralEstimate`] which was calculated before an error
    /// occurred.
    pub fn estimate(&self) -> &IntegralEstimate<T> {
        &self.estimate
    }
}

/// The kind of [`IntegrationError`] which occurred.
///
/// Errors typically occur during integration due to, for example, difficult integration regions,
/// non-integrable singularities, the maximum number of iterations being reached, etc.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub enum IntegrationErrorKind {
    /// The user supplied maximum number of adaptive iterations was reached before the requested
    /// numerical integration tolerance could be achieved.
    MaximumIterationsReached(usize),

    /// Roundoff error detected which prevents the requested tolerance from being achieved.
    /// The approximated numerical error may be under-estimated.
    RoundoffErrorDetected,

    /// Extremely bad integrand behaviour. Possible non-integrable singularity, divergence, or
    /// discontinuity detected between the upper and lower limits.
    BadIntegrandBehaviour(Limits),

    /// The numerical integration routine is not converging. Roundoff error is detected in the
    /// extrapolation table. It is assumed that the requested tolerance cannot be achieved and the
    /// returned result is the best which can be obtained.
    DoesNotConverge,

    /// The integral is probably divergent or slowly convergent. NOTE: divergence can also occur
    /// with any other error kind.
    DivergentOrSlowlyConverging,

    /// The integration Workspace was not properly initialised. This error should not be possible
    /// in user code.
    UninitialisedWorkspace,
}

impl<T: ScalarF64> fmt::Display for IntegrationError<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let estimate = self.estimate();
        let result = estimate.result();
        let error = estimate.error();
        let iterations = estimate.iterations();
        match self.kind() {
            IntegrationErrorKind::MaximumIterationsReached(i) => {
                write!(f, "Maximum number of iterations/subdivisions  ({i}) reached. Try increasing max_iterations. If this yields no improvement it is advised to analyse the integrand to determine integration difficulties. If the position of a local difficulty can be determined, one may gain from splitting the total integration interval and calling the integrator on each sub-interval.\nresult:\t{result}\nerror\t{error:.10e}\niterations:\t{iterations}.")
            }

            IntegrationErrorKind::RoundoffErrorDetected => {
                write!(f, "Roundoff error detected. This prevents the requested tolerance from being achieved and the returned error may be under-estimated.\nresult:\t{result}\nerror\t{error:.10e}\niterations:\t{iterations}.")
            }

            IntegrationErrorKind::BadIntegrandBehaviour(limits) => {
                let lower = limits.lower();
                let upper = limits.upper();
                write!(f, "Extremely bad integrand behaviour. Possible non-integrable singularity, divergence, or discontinuity detected between ({lower},{upper}).\nresult:\t{result}\nerror\t{error:.10e}\niterations:\t{iterations}.\nTry reducing the requested tolerance.")
            }

            IntegrationErrorKind::DoesNotConverge => {
                write!(f, "The algorithm does not converge. Roundoff error is detected in the extrapolation table. It is assumed that the requested tolerance cannot be achieved and the returned result is the best which can be obtained.\nresult:\t{result}\nerror\t{error:.10e}\niterations:\t{iterations}.")
            }

            IntegrationErrorKind::DivergentOrSlowlyConverging => {
                write!(f, "The integral is probably divergent or slowly convergent. NOTE: divergence can also occur with any other error kind.\nresult:\t{result}\nerror\t{error:.10e}\niterations:\t{iterations}.")
            }

            IntegrationErrorKind::UninitialisedWorkspace => {
                write!(f, "The integration Workspace was not properly initialised. This error should not be possible. If this error is returned, contact the crate maintainers.")
            }
        }
    }
}

impl<T: ScalarF64> error::Error for IntegrationError<T> {}
