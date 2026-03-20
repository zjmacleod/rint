# rint

Numerical integration routines written in Rust. 

# Overview

This library contains numerical integration routines for both functions of one dimension (see [`quadrature`]) and functions with multiple dimensions up to $N = 15$ (see [`multi`]). The basic principle of the library is to expose all of the routines through the trait system. Each of the one- and multi-dimensional integrators take as a parameter a type implementing the trait [`Integrand`].

The integration routines attempt to make approximations to integrals such as the one-dimensional integral,
```math
I = \int_{b}^{a} f(x) dx
```
or the $N$-dimensional
```math
I = \int_{\Sigma_{N}} f(\mathbf{x}) d\mathbf{x}
```
where $\mathbf{x} = (x_{1}, x_{2}, \dots, x_{N})$ and $\Sigma_{N}$ is an $N$-dimensional hypercube. The functions $f(x)$ and $f(\mathbf{x})$ can be real valued, with return type [`f64`] _or_ complex valued with return type [`Complex<f64>`]. The numerical integration routines approximate the integral of a function by performing a weighted sum of the function evaluated at defined points/abscissae. For example, in the one-dimensional case,
```math
I = \int_{b}^{a} f(x) dx \approx \sum_{i = 1}^{n} W_{i} f(X_{i}) = I_{n}
```
where the $X_{i}$ and $W_{i}$ are the rescaled abscissae and weights,
```math
X_{i} = \frac{b + a + (a - b) x_{i}}{2} ~~~~~~~~ W_{i} = \frac{(a - b) w_{i}}{2}
```
Integration rules have a polynomial order. Rules of higher polynomial order use more points/abscissae and weights for evaluating the sum.

The library contains both non-adaptive and adaptive integration routines. In the case of an adaptive routine, the user supplies an error [`Tolerance`] which acts as a goal for the numerical error estimate. The numerical integration is considered successful when the numerical error is less than the user supplied tolerance. On success, the output of the numerical integration is an [`IntegralEstimate`].


# Traits

The primary entry point for the library is the [`Integrand`] trait. A type implementing the [`Integrand`] trait represents a real or complex valued function which is to be integrated. The trait requires the definition of two associated types, [`Integrand::Point`] and [`Integrand::Scalar`], and an implementation of the method [`Integrand::evaluate`].

- [`Integrand::Point`]: This associated type defines the point at which the function is to be evaluated, and determines the types of numerical integrators which are available to the user to integrate the function. Integrators are provided for univariate functions $f(x)$ through the associated type `Point=f64`, while integrators for multivariate functions $f(\mathbf{x})$ are provided through the associated type `Point=[f64;N]` where $N$ is the dimensionality of the point $\mathbf{x}=(x_{1},\dots,x_{N})$ which is limited to $2 \le N \le 15$.

- [`Integrand::Scalar`]: This is the output type of the function to be integrated. A _real_ valued function should have the output type `Scalar=`[`f64`], while a _complex_ valued function should have output type `Scalar=`[`Complex<f64>`].

- [`Integrand::evaluate`]: The trait requires an implementation of an `evaluate` method, which defines how the function takes the input [`Integrand::Point`] and turns this into the output type [`Integrand::Scalar`]. In other words, this method tells the integrators how to evaluate the function at a point, allowing the integration to be done.


As an example, consider probability density function of a normal distribution,
```math
f(x) = \frac{1}{\sqrt{2 \pi \sigma^{2}}} e^{- \frac{(x - \mu)^{2}}{2 \sigma^{2}}}
```
which has mean $\mu$ and standard deviation $\sigma$. To integrate this function, one first implements the [`Integrand`] trait,
```rust
use rint::Integrand;

struct NormalDist {
    mean: f64,
    standard_dev: f64,
}

impl Integrand for NormalDist {
    type Scalar = f64;
    fn evaluate(&self, x: f64) -> Self::Scalar {
        let mu = self.mean;
        let sigma = self.standard_dev;

        let prefactor = 1.0 / (2.0 * std::f64::consts::PI * sigma.powi(2));
        let exponent = - (x - mu).powi(2) / (2.0 * sigma.powi(2));

        prefactor * exponent.exp()
    }
}
```
The type `NormalDist` can then be passed as a parameter to a one-dimensional integration routine in the module [`quadrature`] and integrated over a given (possibly infinite) region.


# Modules

## [`quadrature`]

The [`quadrature`] module provides a number of numerical quadrature integration routines of a function in one dimension. The routines are based primarily on the algorithms presented in the Fortran library [QUADPACK] (Piessens, de Doncker-Kapenga, Ueberhuber and Kahaner) and reimplemented in the [GNU scientific library] (GSL) (Gough, Alken, Gonnet, Holoborodko, Griessinger). The integrators use Gauss-Kronrod integration rules, combining two rules of different order for efficient estimation of the numerical error. The rules for an $n$-point Gauss-Kronrod rule contain $m = (n - 1) / 2$ abscissae _shared_ by the Gaussian and Kronrod rules and an extended set of $n - m$ Kronrod abscissae. Thus, only $n$ total function evaluations are required for both the integral and error estimates.

The `rint` library departs from the naming conventions of the [QUADPACK] and [GSL], but provides a selection of comparable implementations:

- [`quadrature::Basic`]: A non-adaptive routine which applies a provided Gauss-Kronrod integration [`quadrature::Rule`] to a function exactly once.

- [`quadrature::Adaptive`]: A one-dimensional adaptive routine based on the `qag.f` [QUADPACK] and `qag.c` [GSL] implementations. The integration is region is bisected into subregions and an initial estimate is performed. Upon each further iteration of the routine the subregion with the highest numerical error estimate is bisected again and new estimates are calculated for these newly bisected regions. This concentrates the integration refinement to the regions with highest error, rapidly reducing the numerical error of the routine. Gauss-Kronrod integration [`quadrature::Rule`]s are provided of various order to use with the adaptive algorithm.

- [`quadrature::AdaptiveSingularity`]: A one-dimensional adaptive routine based on the `qags.f` [QUADPACK] and `qags.c` [GSL] implementations. Adaptive routines concentrate new subintervals around the region of highest error. If this region contains an integrable singularity, then the adaptive routine of [`quadrature::Adaptive`] may fail to obtain a suitable estimate. However, this can be combined with an extrapolation proceedure such as the Wynn epsilon-algorithm to extrapolate the value at these integrable singularities and provide a suitable estimate. As well as handling integrable singularities, [`quadrature::AdaptiveSingularity`] can be used to calculate integrals with infinite or semi-infinite bounds by using the appropriate constructors.

## [`multi`]

The [`multi`] module provides numerical integration routines for integrating functions with dimensionality between $2 \le N \le 15$. The functions can be either real-valued or complex, and are integrated over an $N$-dimensional hypercube. The routines are based primarily on the [DCUHRE] FORTRAN library (Bernsten, Espelid, Genz), which use sets of fully-symmetric integration rules to obtain integral and error estimates. Unlike the original algorithm the routines presented in [`multi`] currently only operate on a single function _not_ a vector of functions. The module provides two classes of routine:

- [`multi::Basic`]: A non-adaptive routine which applies a fully-symmetric integration rule [`multi::Rule`] to a multi-dimensional function exactly once.

- [`multi::Adaptive`]: A $2 \le N \le 15$ dimensional adaptive routine with a similar approach to the one-dimensional adaptive routines found in [`quadrature`]. On each iteration of the algorithm the axis along which the largest contribution to the error estimate was obtained is used as the bisection axis to bisect the integration region and then calculate new estimates for these newly bisected volumes. This concentrates the integration refinement to the regions with highest error, rapidly reducing the numerical error of the routine. The algorithm uses fully-symmetric integration rules, [`multi::Rule`], of varying order and generality. These are generated through the `Rule*::generate` constructors.


# One-dimensional example

The following example integrates the function
```math
f(x) = \frac{\log(x)}{(1 + 100 x^{2})}
```
over the semi-infinite interval $0 < x < \infty$ using an adaptive integration routine with singularity detection (see [`quadrature::AdaptiveSingularity`]). Adapted from the [QUADPACK] & [GSL] numerical integration test suites.


```rust
use rint::{Limits, Integrand, Tolerance};
use rint::quadrature::AdaptiveSingularity;
use std::error::Error;

struct F;

impl Integrand for F {
    type Scalar = f64;
    fn evaluate(&self, x: f64) -> Self::Scalar {
        x.ln() / (1.0 + 100.0 * x.powi(2))
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let function = F;
    let lower = 0.0;
    let tolerance = Tolerance::Relative(1.0e-3);
    
    let integrator = AdaptiveSingularity::semi_infinite_upper(&function, lower, tolerance, 1000)?;
    
    let exp_result = -3.616_892_186_127_022_568E-01;
    let exp_error = 3.016_716_913_328_831_851E-06;
    
    let integral = integrator.integrate()?;
    let result = integral.result();
    let error = integral.error();
    
    let tol = 1.0e-9;
    assert!((exp_result - result).abs() / exp_result.abs() < tol);
    assert!((exp_error - error).abs() / exp_error.abs() < tol);
    Ok(())
}
```

# Multi-dimensional example

The following example integtates a 4-dimensional function $f(\mathbf{x})$,
```math
f(\mathbf{x}) = \frac{x_{3}^{2} x_{4} e^{x_{3} x_{4}}}{(1 + x_{1} + x_{2})^{2}}
```
where $\mathbf{x} = (x_{1}, x_{2}, x_{3}, x_{4})$ over an $N=4$ dimensional hypercube $((0,1),(0,1),(0,2),(0,1))$ using a fully-symmetric 7-point adaptive algorithm. Adapted from P. van Dooren & L. de Ridder, "An adaptive algorithm for numerical integration over an n-dimensional cube", J. Comp. App. Math., Vol. 2, (1976) 207-217

```rust
use rint::{Limits, Integrand, Tolerance};
use rint::multi::{Adaptive, Rule07};
use std::error::Error;

const N: usize = 4;

struct F;

impl Integrand for F {
    type Point = [f64; N]; 
    type Scalar = f64;
    fn evaluate(&self, coordinates: &[f64; N]) -> Self::Scalar {
        let [x1, x2, x3, x4] = coordinates;
        x3.powi(2) * x4 * (x3 * x4).exp() / (x1 + x2 + 1.0).powi(2)
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    const TARGET: f64 = 5.753_641_449_035_616e-1;
    const TOL: f64 = 1e-2;
    
    let function = F;
    let limits = [
        Limits::new(0.0, 1.0),
        Limits::new(0.0, 1.0),
        Limits::new(0.0, 1.0),
        Limits::new(0.0, 2.0)
    ];
    let rule = Rule07::<N>::generate()?;
    let tolerance = Tolerance::Relative(TOL);
    
    let integrator = Adaptive::new(&function, &rule, limits, tolerance, 10000)?;
    
    let integral = integrator.integrate()?;
    let result = integral.result();
    let error = integral.error();
    
    let actual_error = (result - TARGET).abs();
    let requested_error = TOL * result.abs();
    
    assert!(actual_error < error);
    assert!(error < requested_error);
    Ok(())
}
```

[GNU Scientific Library]: <https://www.gnu.org/software/gsl/doc/html/integration.html>
[GSL]: <https://www.gnu.org/software/gsl/doc/html/integration.html>
[QUADPACK]: <https://www.netlib.org/quadpack/>
[DCUHRE]: <https://dl.acm.org/doi/10.1145/210232.210234>
[`Integrand`]: <https://docs.rs/rint/latest/rint/trait.Integrand.html>
[`Limits`]: <https://docs.rs/rint/latest/rint/struct.Limits.html>
[`Tolerance`]: <https://docs.rs/rint/latest/rint/enum.Tolerance.html>
[`IntegralEstimate`]: <https://docs.rs/rint/latest/rint/struct.IntegralEstimate.html>
[`quadrature`]: <https://docs.rs/rand/latest/rand/quadrature/>
[`quadrature::Rule`]:  <https://docs.rs/rint/latest/rint/quadrature/struct.Rule.html>
[`quadrature::Basic`]:  <https://docs.rs/rint/latest/rint/quadrature/struct.Basic.html>
[`quadrature::Adaptive`]:  <https://docs.rs/rint/latest/rint/quadrature/struct.Adaptive.html>
[`quadrature::AdaptiveSingularity`]:  <https://docs.rs/rint/latest/rint/quadrature/struct.AdaptiveSingularity.html>
[`multi`]: <https://docs.rs/rand/latest/rand/multi/>
[`multi::Rule`]:  <https://docs.rs/rint/latest/rint/multi/struct.Rule.html>
[`multi::Basic`]:  <https://docs.rs/rint/latest/rint/multi/struct.Basic.html>
[`multi::Adaptive`]:  <https://docs.rs/rint/latest/rint/multi/struct.Adaptive.html>
[`f64`]: <https://doc.rust-lang.org/std/primitive.f64.html>
[`Complex<f64>`]: <https://docs.rs/num/latest/num/complex/struct.Complex.html>
