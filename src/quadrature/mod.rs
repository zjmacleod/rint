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
