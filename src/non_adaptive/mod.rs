//! Non-adaptive Gauss-Kronrod quadrature integration routines.
pub(crate) mod gauss_kronrod_15;
pub(crate) mod gauss_kronrod_21;
pub(crate) mod gauss_kronrod_31;
pub(crate) mod gauss_kronrod_41;
pub(crate) mod gauss_kronrod_51;
pub(crate) mod gauss_kronrod_61;
pub(crate) mod integration;

pub use gauss_kronrod_15::GaussKronrod15;
pub use gauss_kronrod_21::GaussKronrod21;
pub use gauss_kronrod_31::GaussKronrod31;
pub use gauss_kronrod_41::GaussKronrod41;
pub use gauss_kronrod_51::GaussKronrod51;
pub use gauss_kronrod_61::GaussKronrod61;
pub use integration::GaussKronrod;
pub use integration::GaussKronrodIntegral;

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
pub trait Rule: private::Sealed {
    /// Iterator type for the Gaussian nodes and weights.
    type Gauss: IntoIterator<Item = f64>;
    /// Iterator type for the Kronrod nodes and weights.
    type Kronrod: IntoIterator<Item = f64>;

    /// Return the Gaussian quadrature nodes `x_i`.
    fn gauss_nodes(&self) -> Self::Gauss;
    /// Return the Gaussian quadrature weights `w_i`.
    fn gauss_weights(&self) -> Self::Gauss;
    /// Return the Kronrod extension nodes `x_i`.
    fn kronrod_nodes(&self) -> Self::Kronrod;
    /// Return the Krondod extension weights `w_i`.
    fn kronrod_weights(&self) -> Self::Kronrod;
}
