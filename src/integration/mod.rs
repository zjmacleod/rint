pub mod adaptive;
pub mod basic;
pub mod singularity;

pub use adaptive::Adaptive;
pub use adaptive::GaussKronrodAdaptive;
pub use basic::Basic;
pub use basic::GaussKronrodBasic;

/// User selected error bound type.
pub enum ErrorBound {
    Absolute(f64),
    Relative(f64),
    Either { absolute: f64, relative: f64 },
}

impl ErrorBound {
    #[must_use]
    pub fn tolerance(&self, integral_value: &f64) -> f64 {
        match *self {
            ErrorBound::Absolute(v) => v,
            ErrorBound::Relative(v) => v * integral_value.abs(),
            ErrorBound::Either { absolute, relative } => {
                f64::max(absolute, relative * integral_value.abs())
            }
        }
    }
}

fn rescale_error(error: f64, result_abs: f64, result_asc: f64) -> f64 {
    let mut error = error.abs();

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
