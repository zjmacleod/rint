mod adaptive;
mod basic;
mod generator;
mod geometry;
mod integrator;
mod region;
pub(crate) mod rule;

pub use adaptive::Adaptive;
pub use basic::Basic;
pub(crate) use integrator::Integrator;
pub(crate) use region::Region;
pub use rule::{Rule, Rule07, Rule09, Rule09N2, Rule11, Rule13};

#[inline]
pub(crate) const fn two_pow_n(n: usize) -> usize {
    let mut exp = n;

    // Never need to check this since NDIM > 2
    //if exp == 0 {
    //    return 1;
    //}
    let mut base = 2;
    let mut acc = 1;

    while exp > 1 {
        if (exp & 1) == 1 {
            acc *= base;
        }
        exp /= 2;
        base *= base;
    }

    acc * base
}

#[inline]
pub(crate) const fn two_pow_n_f64(n: usize) -> f64 {
    two_pow_n(n) as f64
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    #[test]
    fn test_two_pow_n() {
        assert_eq!(2, two_pow_n(1));
        assert_eq!(4, two_pow_n(2));
        assert_eq!(8, two_pow_n(3));
        assert_eq!(16, two_pow_n(4));
        assert_eq!(32, two_pow_n(5));
        assert_eq!(64, two_pow_n(6));
        assert_eq!(128, two_pow_n(7));
        assert_eq!(256, two_pow_n(8));
        assert_eq!(512, two_pow_n(9));
        assert_eq!(1024, two_pow_n(10));
        assert_eq!(2048, two_pow_n(11));
        assert_eq!(4096, two_pow_n(12));
        assert_eq!(8192, two_pow_n(13));
        assert_eq!(16384, two_pow_n(14));
        assert_eq!(32768, two_pow_n(15));
    }

    #[test]
    fn test_two_pow_n_f64() {
        assert_eq!(2.0, two_pow_n_f64(1));
        assert_eq!(4.0, two_pow_n_f64(2));
        assert_eq!(8.0, two_pow_n_f64(3));
        assert_eq!(16.0, two_pow_n_f64(4));
        assert_eq!(32.0, two_pow_n_f64(5));
        assert_eq!(64.0, two_pow_n_f64(6));
        assert_eq!(128.0, two_pow_n_f64(7));
        assert_eq!(256.0, two_pow_n_f64(8));
        assert_eq!(512.0, two_pow_n_f64(9));
        assert_eq!(1024.0, two_pow_n_f64(10));
        assert_eq!(2048.0, two_pow_n_f64(11));
        assert_eq!(4096.0, two_pow_n_f64(12));
        assert_eq!(8192.0, two_pow_n_f64(13));
        assert_eq!(16384.0, two_pow_n_f64(14));
        assert_eq!(32768.0, two_pow_n_f64(15));
    }
}
