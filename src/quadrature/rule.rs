//! Non-adaptive Gauss-Kronrod quadrature integration routines.
use std::fmt::Debug;

pub(crate) mod gauss_kronrod_15;
pub(crate) mod gauss_kronrod_21;
pub(crate) mod gauss_kronrod_31;
pub(crate) mod gauss_kronrod_41;
pub(crate) mod gauss_kronrod_51;
pub(crate) mod gauss_kronrod_61;

pub use gauss_kronrod_15::GaussKronrod15;
pub use gauss_kronrod_21::GaussKronrod21;
pub use gauss_kronrod_31::GaussKronrod31;
pub use gauss_kronrod_41::GaussKronrod41;
pub use gauss_kronrod_51::GaussKronrod51;
pub use gauss_kronrod_61::GaussKronrod61;

mod private {
    pub trait Sealed {}

    impl Sealed for super::GaussKronrod15 {}
    impl Sealed for super::GaussKronrod21 {}
    impl Sealed for super::GaussKronrod31 {}
    impl Sealed for super::GaussKronrod41 {}
    impl Sealed for super::GaussKronrod51 {}
    impl Sealed for super::GaussKronrod61 {}
}

/// Gauss-Kronrod integration rule.
///
/// Gaussian quadrature is a method of approximating the definite integral of a function as a
/// weighted sum of function values at specified points within the domain of integration.
/// An `n` point Gaussian quadrature rule consists of a set of points `x_i` and weights `w_i` for
/// `i = 1,...,n`.
/// Gauss-Kronrod quadrature is a variant of Gaussian quadrature which combines an `n`-point
/// Gaussian rule with a `(2n - 1)`-rule such that the resulting combination is of order `3n + 1`.
/// The difference between the lower-order `n`-point rule and the Kronrod extension are used to
/// estimate the approximate error of the numerical integration.
/// The [`Rule`] trait defines a Gauss-Kronrod quadrature rule, returning the sets of points `x_i`
/// and weights `w_i` for the `n`-point Gauss rule and the corresponding `(2n - 1)`-point Kronrod
/// extension.
pub trait Rule: private::Sealed + Clone + Copy + Debug {
    type Shared: IntoIterator<Item = (f64, f64, f64)>;
    type Extended: IntoIterator<Item = (f64, f64)>;

    const KRONROD_CENTRE: f64;
    const EVALUATIONS: usize;

    fn shared_data(&self) -> Self::Shared;
    fn extended_data(&self) -> Self::Extended;
    fn gauss_centre(&self) -> Option<f64>;
    fn kronrod_centre(&self) -> f64 {
        Self::KRONROD_CENTRE
    }
    fn evaluations(&self) -> usize {
        Self::EVALUATIONS
    }
}
