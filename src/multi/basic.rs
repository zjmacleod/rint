use crate::multi::{Integrator, Rule};
use crate::IntegralEstimate;
use crate::Limits;
use crate::MultiDimensionalIntegrand;
use crate::{InitialisationError, InitialisationErrorKind};

pub struct Basic<I, const NDIM: usize, const FINAL: usize, const TOTAL: usize> {
    function: I,
    rule: Rule<NDIM, FINAL, TOTAL>,
    limits: [Limits; NDIM],
}

impl<I, const NDIM: usize, const FINAL: usize, const TOTAL: usize> Basic<I, NDIM, FINAL, TOTAL>
where
    I: MultiDimensionalIntegrand<NDIM>,
{
    /// Create a new [`Basic`] multi-dimensional integrator.
    ///
    /// The user first defines a `function` which is a `struct` implementing the
    /// [`MultiDimensionalIntegrand`] trait and selects a fully-symmetric multi-dimensional
    /// integration [`Rule`], `rule`, to integrate the function in the hypercube formed by the
    /// [`Limits`], `limits` in each of the `NDIM` integration directions.
    ///
    /// # Errors
    /// Will fail if `NDIM < 2` or `NDIM > 15`. The routines probided in this module are developed
    /// for dimensionalities between `2 <= NDIM <= 15`.
    pub fn new(
        function: I,
        rule: Rule<NDIM, FINAL, TOTAL>,
        limits: [Limits; NDIM],
    ) -> Result<Self, InitialisationError> {
        if NDIM < 2 || NDIM > 15 {
            return Err(InitialisationError::new(
                InitialisationErrorKind::InvalidDimension(NDIM),
            ));
        };
        Ok(Self {
            function,
            rule,
            limits,
        })
    }

    pub fn integrate(&self) -> IntegralEstimate<I::Scalar> {
        let integral = self.integrator().integrate();
        IntegralEstimate::new()
            .with_result(integral.result())
            .with_error(integral.error())
            .with_iterations(1)
            .with_function_evaluations(self.rule.evaluations())
    }

    const fn integrator(&self) -> Integrator<'_, I, NDIM, FINAL, TOTAL> {
        Integrator::new(&self.function, &self.rule, self.limits)
    }
}

#[cfg(test)]
mod testsfoo {
    use super::*;
    use crate::multi::{Rule07, Rule09, Rule11, Rule13};

    // TODO: check the largest error axis selection

    #[test]
    fn compare_basic_7point_with_dcuhre_output_ndim_2() {
        const NDIM: usize = 2;

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let exponent = -x * y;
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.79659970249839818;
            let dcuhre_error = 2.7359932817440052E-005;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-13);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-7);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let exponent = -(x * x + y * y);
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.55774638300971069;
            let dcuhre_error = 1.3049001762496107E-005;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-12);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-7);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let exponent = -((x * x).cos().powi(2) * (y * y).cos().powi(2));
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.51478351071879791;
            let dcuhre_error = 0.37449223525594855;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-10);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-9);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    f64::cos((1.0 + x * x).ln() * (1.0 + y * y).ln())
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.99332124356102158;
            let dcuhre_error = 9.3939393118662768E-005;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-10);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-7);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    (x * x / (2.0 - x.cos())) + (y * y / (2.0 - y.cos()))
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limits = [Limits::new(-2.0, 3.0), Limits::new(1.0, 10.0)];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 887.97458362328393;
            let dcuhre_error = 611.03739938210026;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-7);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-8);
        }
    }

    #[test]
    fn compare_basic_7point_with_dcuhre_output_ndim_3() {
        const NDIM: usize = 3;

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent = -x * y * z;
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.89120896753764178;
            let dcuhre_error = 1.3869851434747299E-005;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-12);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-7);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent = -(x * x + y * y + z * z);
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.41653847897138208;
            let dcuhre_error = 5.8625669977822451E-006;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-12);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 2e-7);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent =
                        -((x * x).cos().powi(2) * (y * y).cos().powi(2) * (z * z).cos().powi(2));
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.57845060677495597;
            let dcuhre_error = 0.47067683202127575;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-11);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-8);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    f64::cos((1.0 + x * x).ln() * (1.0 + y * y).ln() * (1.0 + z * z).ln())
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.99921595669136032;
            let dcuhre_error = 2.5792992219862503E-003;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-12);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-8);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    (x * x / (2.0 - x.cos()))
                        + (y * y / (2.0 - y.cos()))
                        + (z * z / (2.0 - z.cos()))
                }
            }

            let function = Function;
            let rule = Rule07::generate().unwrap();
            let limits = [
                Limits::new(-2.0, 3.0),
                Limits::new(1.0, 10.0),
                Limits::new(0.0, -1.0),
            ];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = -899.77650881070997;
            let dcuhre_error = 543.83282190394414;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-7);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-8);
        }
    }

    #[test]
    fn compare_basic_9point_with_dcuhre_output_ndim_3() {
        const NDIM: usize = 3;

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent = -x * y * z;
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule09::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.89121284475576390;
            let dcuhre_error = 2.8431093456305666E-003;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-13);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent = -(x * x + y * y + z * z);
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule09::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.41653839142919569;
            let dcuhre_error = 7.2117833072453187E-008;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-8);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent =
                        -((x * x).cos().powi(2) * (y * y).cos().powi(2) * (z * z).cos().powi(2));
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule09::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.57854494486523600;
            let dcuhre_error = 4.2597555388569463E-005;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-14);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-10);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    f64::cos((1.0 + x * x).ln() * (1.0 + y * y).ln() * (1.0 + z * z).ln())
                }
            }

            let function = Function;
            let rule = Rule09::generate().unwrap();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.99922153302627292;
            let dcuhre_error = 7.5347754069994284E-004;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-12);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    (x * x / (2.0 - x.cos()))
                        + (y * y / (2.0 - y.cos()))
                        + (z * z / (2.0 - z.cos()))
                }
            }

            let function = Function;
            let rule = Rule09::generate().unwrap();
            let limits = [
                Limits::new(-2.0, 3.0),
                Limits::new(1.0, 10.0),
                Limits::new(0.0, -1.0),
            ];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = -905.91915061432690;
            let dcuhre_error = 144.00666640485426;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-14);
        }
    }

    #[test]
    fn compare_basic_11point_with_dcuhre_output_ndim_3() {
        const NDIM: usize = 3;

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent = -x * y * z;
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule11::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.89121279793446351;
            let dcuhre_error = 1.5894785143767116E-008;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-20);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-20);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent = -(x * x + y * y + z * z);
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule11::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.41653838596655557;
            let dcuhre_error = 3.2348830851556911E-006;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-20);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-20);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    let exponent =
                        -((x * x).cos().powi(2) * (y * y).cos().powi(2) * (z * z).cos().powi(2));
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule11::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.57853848358343507;
            let dcuhre_error = 2.3874271216029219E-003;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-20);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-13);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    f64::cos((1.0 + x * x).ln() * (1.0 + y * y).ln() * (1.0 + z * z).ln())
                }
            }

            let function = Function;
            let rule = Rule11::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.99922120779864210;
            let dcuhre_error = 1.3732744589107246E-007;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-20);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-9);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let z = coordinates[2];
                    (x * x / (2.0 - x.cos()))
                        + (y * y / (2.0 - y.cos()))
                        + (z * z / (2.0 - z.cos()))
                }
            }

            let function = Function;
            let rule = Rule11::generate();
            let limits = [
                Limits::new(-2.0, 3.0),
                Limits::new(1.0, 10.0),
                Limits::new(0.0, -1.0),
            ];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = -904.73349759612938;
            let dcuhre_error = 36.531147257907918;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-20);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-14);
        }
    }

    #[test]
    fn compare_basic_13point_with_dcuhre_output_ndim_2() {
        const NDIM: usize = 2;

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let exponent = -x * y;
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule13::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.79659959929708202;
            let dcuhre_error = 3.2673674642535072E-013;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-20);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-20);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let exponent = -(x * x + y * y);
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule13::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.55774628535131576;
            let dcuhre_error = 9.7456124932917485E-012;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-20);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-20);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    let exponent = -((x * x).cos().powi(2) * (y * y).cos().powi(2));
                    f64::exp(exponent)
                }
            }

            let function = Function;
            let rule = Rule13::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.51478983417297963;
            let dcuhre_error = 4.4473459718302440E-004;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-20);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-12);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    f64::cos((1.0 + x * x).ln() * (1.0 + y * y).ln())
                }
            }

            let function = Function;
            let rule = Rule13::generate();
            let limit = Limits::new(0.0, 1.0);
            let limits = [limit; NDIM];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 0.99331725120688252;
            let dcuhre_error = 8.0707600531214359E-010;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-20);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-20);
        }

        {
            struct Function;

            impl MultiDimensionalIntegrand<NDIM> for Function {
                type Scalar = f64;
                fn evaluate(&self, coordinates: &[f64; NDIM]) -> Self::Scalar {
                    let x = coordinates[0];
                    let y = coordinates[1];
                    (x * x / (2.0 - x.cos())) + (y * y / (2.0 - y.cos()))
                }
            }

            let function = Function;
            let rule = Rule13::generate();
            let limits = [Limits::new(-2.0, 3.0), Limits::new(1.0, 10.0)];

            let integral = Basic::new(function, rule, limits).unwrap();

            let integral_result = integral.integrate();

            let result = integral_result.result();
            let error = integral_result.error();
            let dcuhre_result = 911.85973569354837;
            let dcuhre_error = 129.79778763945998;

            assert!((result - dcuhre_result).abs() / dcuhre_result.abs() < 1e-15);
            assert!((error - dcuhre_error).abs() / dcuhre_error.abs() < 1e-20);
        }
    }
}
