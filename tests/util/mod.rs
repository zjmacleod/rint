use std::f64::consts::PI;

use gauss_kronrod_integration::Integrand;

/* These are the test functions from table 4.1 of the QUADPACK book */

/* f1(x) = x^alpha * log(1/x) */
/* integ(f1,x,0,1) = 1/(alpha + 1)^2 */
pub(crate) struct Function1 {
    pub(crate) alpha: f64,
}

impl Integrand for Function1 {
    fn evaluate(&self, x: &f64) -> f64 {
        let alpha = self.alpha;
        x.powf(alpha) * (1.0 / x).ln()
    }
}

/* f2(x) = 4^-alpha / ((x-pi/4)^2 + 16^-alpha) */
/* integ(f2,x,0,1) = arctan((4-pi)4^(alpha-1)) + arctan(pi 4^(alpha-1)) */
pub(crate) struct Function2 {
    pub(crate) alpha: f64,
}

impl Integrand for Function2 {
    fn evaluate(&self, x: &f64) -> f64 {
        let alpha = self.alpha;
        4.0_f64.powf(-alpha) / ((x - PI / 4.0).powi(2) + 16.0_f64.powf(-alpha))
    }
}

/* f3(x) = cos(2^alpha * sin(x)) */
/* integ(f3,x,0,pi) = pi J_0(2^alpha) */
pub(crate) struct Function3 {
    pub(crate) alpha: f64,
}

impl Integrand for Function3 {
    fn evaluate(&self, x: &f64) -> f64 {
        let alpha = self.alpha;
        (2.0_f64.powf(alpha) * x.sin()).cos()
    }
}

/* f16(x) = x^(alpha - 1) / (1 + 10 x)^2*/
pub(crate) struct Function16 {
    pub(crate) alpha: i32,
}

impl Integrand for Function16 {
    fn evaluate(&self, x: &f64) -> f64 {
        let alpha = self.alpha;
        if x.to_bits() == 0f64.to_bits() && alpha == 1 {
            1.0
        } else if x.to_bits() == 0f64.to_bits() && alpha > 1 {
            0.0
        } else {
            x.powi(alpha - 1) / (1.0 + 10.0 * x).powi(2)
        }
    }
}

pub(crate) fn test_relative_error(
    calculated: f64,
    target: f64,
    relative_error: f64,
    description: &str,
) -> Result<(), String> {
    if calculated.is_nan() || target.is_nan() {
        if calculated.is_nan() != calculated.is_nan() {
            return Err(format!(
                "Failed test {description}: calculated.is_nan() != target.is_nan()"
            ));
        }
    }

    if calculated.is_infinite() || target.is_infinite() {
        if calculated.is_infinite() != calculated.is_infinite() {
            return Err(format!(
                "Failed test {description}: calculated.is_infinite() != target.is_infinite()"
            ));
        }
    }

    if (target > 0.0 && target < f64::MIN_POSITIVE) || (target < 0.0 && target > -f64::MIN_POSITIVE)
    {
        return Err(format!(
            "Failed test {description}: target value smaller than f64::MIN_POSITIVE"
        ));
    }

    if target != 0.0 {
        if (calculated - target).abs() / target.abs() > relative_error {
            return Err(format!(
                "Failed test {description}: calculated relative error is larger than target"
            ));
        }
    } else {
        if calculated.abs() > relative_error {
            return Err(format!("Failed test {description}: calculated integral value is smaller than calculated relative error"));
        }
    }

    Ok(())
}

pub(crate) fn test_absolute_error(
    calculated: f64,
    target: f64,
    absolute_error: f64,
    description: &str,
) -> Result<(), String> {
    if calculated.is_nan() || target.is_nan() {
        if calculated.is_nan() != calculated.is_nan() {
            return Err(format!(
                "Failed test {description}: calculated.is_nan() != target.is_nan()"
            ));
        }
    }

    if calculated.is_infinite() || target.is_infinite() {
        if calculated.is_infinite() != calculated.is_infinite() {
            return Err(format!(
                "Failed test {description}: calculated.is_infinite() != target.is_infinite()"
            ));
        }
    }

    if (target > 0.0 && target < f64::MIN_POSITIVE) || (target < 0.0 && target > -f64::MIN_POSITIVE)
    {
        return Err(format!(
            "Failed test {description}: target value smaller than f64::MIN_POSITIVE"
        ));
    } else {
        if (calculated - target).abs() > absolute_error {
            return Err(format!("Failed test {description}: calculated absolute error is smaller than target relative error"));
        }
    }
    Ok(())
}

pub(crate) fn test_int(calculated: usize, target: usize, description: &str) -> Result<(), String> {
    if calculated == target {
        Ok(())
    } else {
        Err(format!(
            "Failed test {description}: calculated value: {calculated} target value: {target}"
        ))
    }
}
