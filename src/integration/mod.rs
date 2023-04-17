pub mod basic;
pub mod non_adaptive;

pub use basic::GaussKronrod;
pub use basic::GaussKronrodBasic;

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
