//! Numerical integration routines written in Rust.
//!
//! This library exposes various numerical integration routines written in pure Rust.
//! The algorithms are heavily inspired by the GNU Scientific Library & QUADPACK numerical integration routines.
//!
#![deny(clippy::pedantic)]
#![allow(
    unused,
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
    + ops::Mul<f64, Output = Self>
    + ops::Div<f64, Output = Self>
    + for<'a> ops::Mul<&'a f64, Output = Self>
    + for<'a> ops::Div<&'a f64, Output = Self>
    + ops::AddAssign<Self>
    + ops::Add<Self>
    + ops::Sub<Self>
    + for<'a> ops::AddAssign<&'a Self>
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
/// * [`Tolerance::Absolute(abserr)`] specifies an absolute error and returns final [`IntegralEstimate`] when `error <= abserr`,
/// * [`Tolerance::Relative(relerr)`] specifies a relative error and returns final [`IntegralEstimate`] when `error <= relerr * abs(result)`,  
/// * [`Tolerance::Either{ abserr, relerr }`] to return a result as soon as _either_ the relative or absolute error bound has been satisfied.
///
/// [`Tolerance::Absolute(abserr)`]: crate::quadrature::Tolerance#variant.Absolute
/// [`Tolerance::Relative(relerr)`]: crate::quadrature::Tolerance#variant.Relative
/// [`Tolerance::Either{ abserr, relerr }`]: crate::quadrature::Tolerance#variant.Either
pub enum Tolerance {
    /// Specify an absolute error and returns final [`IntegralEstimate`] when `error <= abserr`.
    Absolute(f64),
    /// Specify a relative error and return final [`IntegralEstimate`] when `error <= relerr * abs(result)`.
    Relative(f64),
    /// Specify _both_ an absolute and relative error and return final [`IntegralEstimate`] as soon as _either_ the relative or absolute error bound has been satisfied.
    Either { absolute: f64, relative: f64 },
}

impl Tolerance {
    // TODO make const when abs() is const 1.85
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

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub enum Error<T: ScalarF64> {
    InitialisationError(InitialisationError),
    IntegrationError(IntegrationError<T>),
}

impl<T: ScalarF64> From<InitialisationError> for Error<T> {
    fn from(other: InitialisationError) -> Self {
        Error::InitialisationError(other)
    }
}

impl<T: ScalarF64> From<IntegrationError<T>> for Error<T> {
    fn from(other: IntegrationError<T>) -> Self {
        Error::IntegrationError(other)
    }
}

impl<T: ScalarF64> error::Error for Error<T> {}

impl<T: ScalarF64> fmt::Display for Error<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Error::InitialisationError(err) => err.fmt(f),
            Error::IntegrationError(err) => err.fmt(f),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct InitialisationError {
    kind: InitialisationErrorKind,
}

impl InitialisationError {
    pub(crate) const fn new(kind: InitialisationErrorKind) -> Self {
        Self { kind }
    }

    pub(crate) const fn kind(&self) -> InitialisationErrorKind {
        self.kind
    }
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub(crate) enum InitialisationErrorKind {
    /// The absolute tolerance bound `tol` requested for the adaptive integration routines must satisfy `tol > 0.0`.
    AbsoluteBoundNegativeOrZero(f64),

    /// The relative tolerance bound `tol` requested for the adaptive integration routines must satisfy `tol > 50.0 * f64::EPSILON`.
    RelativeBoundTooSmall(f64),

    /// An invalid tolerance was requested for the adaptive integration routine.
    /// The absolute tolerance bound `abs_tol` must satisfy `abs_tol > 0.0`.
    /// The relative tolerance bound `rel_tol` must satisfy satisfy `rel_tol > 50.0 * f64::EPSILON`.
    InvalidTolerance { absolute: f64, relative: f64 },

    /// An invalid integration dimensionality `NDIM` was used in a multi-dimensional integration.
    /// This library only provides multi-dimensional integration routines suitable for dimensions
    /// `2 <= NDIM <= 15`.
    InvalidDimension(usize),

    /// An invalid integration dimensionality `NDIM` was used to construct the 7-point multi-dimensional
    /// integration rule [`Rule07`].
    /// This rule is only suitable for dimensions `2 <= NDIM <= 15`.
    InvalidDimensionForRule07(usize),

    /// An invalid integration dimensionality `NDIM` was used to construct the 9-point multi-dimensional
    /// integration rule [`Rule09`].
    /// This rule is only suitable for dimensions `3 <= NDIM <= 15`.
    InvalidDimensionForRule09(usize),

    /// An invalid integration dimensionality `NDIM` was used to construct the 11-point multi-dimensional
    /// integration rule [`Rule11`].
    /// This rule is only suitable for dimensions `NDIM == 3`.
    InvalidDimensionForRule11(usize),

    /// An invalid integration dimensionality `NDIM` was used to construct the 13-point multi-dimensional
    /// integration rule [`Rule13`].
    /// This rule is only suitable for dimensions `NDIM == 2`.
    InvalidDimensionForRule13(usize),
    ///// The minimum number of function evaluations, `min`, to be used in a multi-dimensional integration
    ///// routine was less than the supplied maximum number of function evaluations, `max`. The
    ///// minimum number of function evaluations must satisfy `min < max`.
    //MinimumPointsLargerThanMaximum { min: usize, max: usize },

    ///// The maximum number of function evaluations, `max`, to be used in a multi-dimensional integration
    ///// routine was less than the required minimum, `minimum_max`, for the given integral dimensionality `NDIM` and integration rule.
    ///// [`Rule07`] requires a `max > minimum_max = 1 + 6 NDIM + 2 NDIM (NDIM - 1) + 2^NDIM`.
    ///// [`Rule09`] requires a `max > minimum_max = 1 + 8 NDIM + 6 NDIM (NDIM - 1) + 4 NDIM (NDIM - 1) (NDIM - 1) / 3 + 2^NDIM`.
    ///// [`Rule11`] requires a `max > minimum_max = 127`.
    ///// [`Rule07`] requires a `max > minimum_max = 65`.
    //MaximumPointsInsufficient {
    //    ndim: usize,
    //    max: usize,
    //    minimum_max: usize,
    //},
}

impl fmt::Display for InitialisationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self.kind() {
                InitialisationErrorKind::AbsoluteBoundNegativeOrZero(tol) => {
                    "Invalid tolerance requested:\n`Absolute({tol})`\nThe absolute tolerance bound `tol` requested for the adaptive integration routines must satisfy `tol > 0.0`."
                }

                InitialisationErrorKind::RelativeBoundTooSmall(tol) => {
                    "Invalid tolerance requested:\n`Relative({tol})`\nThe relative tolerance bound `tol` requested for the adaptive integration routines must satisfy `tol > 50.0 * f64::EPSILON`."
                }

                InitialisationErrorKind::InvalidTolerance { absolute, relative } => {
                    "Invalid tolerance requested:\n`Either '{{' abs_tol: {absolute}, rel_tol: {relative} '}}'`\nThe absolute tolerance bound `abs_tol` requested for the adaptive integration routines must satisfy `abs_tol > 0.0`. The relative tolerance bound `rel_tol` requested for the adaptive integration routines must satisfy `rel_tol > 50.0 * f64::EPSILON`."
                }

                InitialisationErrorKind::InvalidDimension(ndim) => {
                    "An invalid integration dimensionality `NDIM = {ndim}` was used in a multi-dimensional integration.  This library only provides multi-dimensional integration routines suitable for dimensions `2 <= NDIM <= 15`."
                }

                InitialisationErrorKind::InvalidDimensionForRule07(ndim) => {
                    "An invalid integration dimensionality `NDIM = {ndim}` was used to construct the 7-point multi-dimensional integration rule [`Rule07`]. This rule is only suitable for dimensions `2 <= NDIM <= 15`.  "
                }

                InitialisationErrorKind::InvalidDimensionForRule09(ndim) => {
                    "An invalid integration dimensionality `NDIM = {ndim}` was used to construct the 9-point multi-dimensional integration rule [`Rule09`]. This rule is only suitable for dimensions `3 <= NDIM <= 15`.  "
                }

                InitialisationErrorKind::InvalidDimensionForRule11(ndim) => {
                    "An invalid integration dimensionality `NDIM = {ndim}` was used to construct the 11-point multi-dimensional integration rule [`Rule11`]. This rule is only suitable for dimensions `NDIM == 3`."
                }

                InitialisationErrorKind::InvalidDimensionForRule13(ndim) => {
                    "An invalid integration dimensionality `NDIM = {ndim}` was used to construct the 13-point multi-dimensional integration rule [`Rule13`]. This rule is only suitable for dimensions `NDIM == 2`."
                } //    InitialisationErrorKind::MinimumPointsLargerThanMaximum { min, max } => {
                  //            "The minimum number of function evaluations, `min = {min}`, to be used in a multi-dimensional integration routine is larger than the supplied maximum number of function evaluations, `max = {max}`. The minimum number of function evaluations must satisfy `min < max`."
                  //    }

                  //    InitialisationErrorKind::MaximumPointsInsufficient {
                  //        ndim,
                  //        max,
                  //        minimum_max,
                  //    } => {
                  //     "The maximum number of function evaluations, `max = {max}`, to be used in a multi-dimensional integration routine was less than the required minimum, `minimum_max = {minimum_max}`, for the given integral dimensionality `NDIM = {ndim}` and integration rule.\n- [`Rule07`] requires a `max > minimum_max = 1 + 6 NDIM + 2 NDIM (NDIM - 1) + 2^NDIM`.\n- [`Rule09`] requires a `max > minimum_max = 1 + 8 NDIM + 6 NDIM (NDIM - 1) + 4 NDIM (NDIM - 1) (NDIM - 1) / 3 + 2^NDIM`.\n- [`Rule11`] requires a `max > minimum_max = 127`.\n- [`Rule07`] requires a `max > minimum_max = 65`."
                  //    }
            }
        )
    }
}

impl error::Error for InitialisationError {}

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
    pub(crate) const fn kind(&self) -> IntegrationErrorKind {
        self.kind
    }

    /// Return a reference to the best [`IntegralEstimate`] which was calculated before an error
    /// occurred.
    pub fn estimate(&self) -> IntegralEstimate<T> {
        self.estimate
    }
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub(crate) enum IntegrationErrorKind {
    /// The user supplied maximum number of adaptive iterations was reached before the
    /// requested numerical integration tolerance could be achieved.
    MaximumIterationsReached(usize),

    /// Roundoff error detected which prevents the requested tolerance from being achieved.
    /// The approximated numerical error may be under-estimated.
    RoundoffErrorDetected,

    /// Extremely bad integrand behaviour. Possible non-integrable singularity, divergence, or discontinuity detected between the upper and lower limits.
    BadIntegrandBehaviour(Limits),

    /// The numerical integration routine is not converging.
    /// Roundoff error is detected in the extrapolation table.
    /// It is assumed that the requested tolerance cannot be achieved and the returned result is the best which can be obtained.
    DoesNotConverge,

    /// The integral is probably divergent or slowly convergent. NOTE: divergence can also occur with any other error kind.
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
        write!(
            f,
            "{}",
            match self.kind() {
                IntegrationErrorKind::MaximumIterationsReached(i) => {
                    "Maximum number of iterations/subdivisions  ({i}) reached. Try increasing max_iterations. If this yields no improvement it is advised to analyse the integrand to determine integration difficulties. If the position of a local difficulty can be determined, one may gain from splitting the total integration interval and calling the integrator on each sub-interval.\nresult:\t{result}\nerror\t{error:.10e}\niterations:\t{iterations}."
                }

                IntegrationErrorKind::RoundoffErrorDetected => {
                    "Roundoff error detected. This prevents the requested tolerance from being achieved and the returned error may be under-estimated.\nresult:\t{result}\nerror\t{error:.10e}\niterations:\t{iterations}."
                }

                IntegrationErrorKind::BadIntegrandBehaviour(limits) => {
                    let lower = limits.lower();
                    let upper = limits.upper();
                    "Extremely bad integrand behaviour. Possible non-integrable singularity, divergence, or discontinuity detected between ({lower},{upper}).\nresult:\t{result}\nerror\t{error:.10e}\niterations:\t{iterations}.\nTry reducing the requested tolerance."
                }

                IntegrationErrorKind::DoesNotConverge => {
                    "The algorithm does not converge. Roundoff error is detected in the extrapolation table. It is assumed that the requested tolerance cannot be achieved and the returned result is the best which can be obtained.\nresult:\t{result}\nerror\t{error:.10e}\niterations:\t{iterations}."
                }

                IntegrationErrorKind::DivergentOrSlowlyConverging => {
                    "The integral is probably divergent or slowly convergent. NOTE: divergence can also occur with any other error kind.\nresult:\t{result}\nerror\t{error:.10e}\niterations:\t{iterations}."
                }

                IntegrationErrorKind::UninitialisedWorkspace => {
                    "The integration Workspace was not properly initialised. This error should not be possible. If this error is returned, contact the crate maintainers."
                }
            }
        )
    }
}

impl<T: ScalarF64> error::Error for IntegrationError<T> {}
