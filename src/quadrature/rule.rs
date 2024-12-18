//! Non-adaptive Gauss-Kronrod quadrature integration routines.
use crate::quadrature::basic::Region;
use crate::quadrature::rescale_error;
use crate::{Integrand, Limits};
use std::fmt::Debug;

pub(crate) mod gauss_kronrod_15;
pub(crate) mod gauss_kronrod_21;
pub(crate) mod gauss_kronrod_31;
pub(crate) mod gauss_kronrod_41;
pub(crate) mod gauss_kronrod_51;
pub(crate) mod gauss_kronrod_61;

#[derive(Clone, Copy)]
pub(crate) struct SharedData {
    point: f64,
    gauss: f64,
    kronrod: f64,
}

impl SharedData {
    pub(crate) const fn new(point: f64, gauss: f64, kronrod: f64) -> Self {
        Self {
            point,
            gauss,
            kronrod,
        }
    }
    pub(crate) const fn point(&self) -> f64 {
        self.point
    }
    pub(crate) const fn gauss(&self) -> f64 {
        self.gauss
    }
    pub(crate) const fn kronrod(&self) -> f64 {
        self.kronrod
    }
}

#[derive(Clone, Copy)]
pub(crate) struct ExtendedData {
    point: f64,
    kronrod: f64,
}

impl ExtendedData {
    pub(crate) const fn new(point: f64, kronrod: f64) -> Self {
        Self { point, kronrod }
    }
    pub(crate) const fn point(&self) -> f64 {
        self.point
    }
    pub(crate) const fn kronrod(&self) -> f64 {
        self.kronrod
    }
}

#[derive(Clone, Copy, Debug)]
pub enum Rule {
    GaussKronrod15,
    GaussKronrod21,
    GaussKronrod31,
    GaussKronrod41,
    GaussKronrod51,
    GaussKronrod61,
}

impl Rule {
    pub(crate) fn evaluations(self) -> usize {
        match self {
            Rule::GaussKronrod15 => gauss_kronrod_15::GaussKronrod15::evaluations(),
            Rule::GaussKronrod21 => gauss_kronrod_21::GaussKronrod21::evaluations(),
            Rule::GaussKronrod31 => gauss_kronrod_31::GaussKronrod31::evaluations(),
            Rule::GaussKronrod41 => gauss_kronrod_41::GaussKronrod41::evaluations(),
            Rule::GaussKronrod51 => gauss_kronrod_51::GaussKronrod51::evaluations(),
            Rule::GaussKronrod61 => gauss_kronrod_61::GaussKronrod61::evaluations(),
        }
    }

    pub(crate) fn integrate<I: Integrand>(self, limits: &Limits, function: &I) -> Region {
        match self {
            Rule::GaussKronrod15 => gauss_kronrod_15::GaussKronrod15::integrate(limits, function),
            Rule::GaussKronrod21 => gauss_kronrod_21::GaussKronrod21::integrate(limits, function),
            Rule::GaussKronrod31 => gauss_kronrod_31::GaussKronrod31::integrate(limits, function),
            Rule::GaussKronrod41 => gauss_kronrod_41::GaussKronrod41::integrate(limits, function),
            Rule::GaussKronrod51 => gauss_kronrod_51::GaussKronrod51::integrate(limits, function),
            Rule::GaussKronrod61 => gauss_kronrod_61::GaussKronrod61::integrate(limits, function),
        }
    }
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
pub(crate) trait Integrate {
    type Shared: IntoIterator<Item = SharedData>;
    type Extended: IntoIterator<Item = ExtendedData>;

    const KRONROD_CENTRE: f64;
    const EVALUATIONS: usize;

    fn shared_data() -> Self::Shared;
    fn extended_data() -> Self::Extended;
    fn gauss_centre() -> Option<f64>;
    fn kronrod_centre() -> f64 {
        Self::KRONROD_CENTRE
    }
    fn evaluations() -> usize {
        Self::EVALUATIONS
    }

    fn integrate<I: Integrand>(limits: &Limits, function: &I) -> Region {
        let centre = limits.centre();
        let half_length = limits.half_width();
        let abs_half_length = half_length.abs();
        let f_centre = function.evaluate(centre);

        let initial_kronrod = Self::kronrod_centre() * f_centre;
        let initial_gauss = if let Some(v) = Self::gauss_centre() {
            v * f_centre
        } else {
            0.0
        };
        let initial_abs = initial_kronrod.abs();

        // XXX Can we get rid of this additional allocation?
        // Vec<(kronrod_weight, (rate_plus, rate_minus))>
        let mut function_values: Vec<(f64, (f64, f64))> = Vec::with_capacity(61);

        let (gauss_result, kronrod_shared, abs_shared) = Self::shared_data()
            .into_iter()
            .map(|data| {
                let point = data.point();
                let gauss = data.gauss();
                let kronrod = data.kronrod();
                let abscissa = half_length * point;
                let rate_plus = function.evaluate(centre + abscissa);
                let rate_minus = function.evaluate(centre - abscissa);
                let rate = rate_plus + rate_minus;
                let rate_abs = rate_plus.abs() + rate_minus.abs();
                function_values.push((kronrod, (rate_plus, rate_minus)));
                (gauss * rate, kronrod * rate, kronrod * rate_abs)
            })
            .fold((initial_gauss, initial_kronrod, initial_abs), |a, v| {
                (a.0 + v.0, a.1 + v.1, a.2 + v.2)
            });

        let (kronrod_result, abs_result) = Self::extended_data()
            .into_iter()
            .map(|data| {
                let point = data.point();
                let kronrod = data.kronrod();
                let abscissa = half_length * point;
                let rate_plus = function.evaluate(centre + abscissa);
                let rate_minus = function.evaluate(centre - abscissa);
                let rate = rate_plus + rate_minus;
                let rate_abs = rate_plus.abs() + rate_minus.abs();
                function_values.push((kronrod, (rate_plus, rate_minus)));
                (kronrod * rate, kronrod * rate_abs)
            })
            .fold((kronrod_shared, abs_shared), |a, v| (a.0 + v.0, a.1 + v.1));

        let mean = kronrod_result * 0.5;

        let initial_asc = Self::kronrod_centre() * (f_centre - mean).abs();

        let asc_result = function_values
            .iter()
            .map(|(k, (rp, rm))| k * ((rp - mean).abs() + (rm - mean).abs()))
            .fold(initial_asc, |a, v| a + v);

        let error = (kronrod_result - gauss_result) * half_length;

        let result = kronrod_result * half_length;
        let result_abs = abs_result * abs_half_length;
        let result_asc = asc_result * abs_half_length;
        let error = rescale_error(error, result_abs, result_asc);

        Region::unevaluated()
            .with_error(error)
            .with_result(result)
            .with_result_abs(result_abs)
            .with_result_asc(result_asc)
            .with_limits(*limits)
    }
}
