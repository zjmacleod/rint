#[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
pub struct Limits {
    lower: f64,
    upper: f64,
}

impl Limits {
    #[must_use]
    pub const fn new(lower: f64, upper: f64) -> Self {
        Self { lower, upper }
    }

    #[must_use]
    pub const fn lower(&self) -> f64 {
        self.lower
    }

    #[must_use]
    pub const fn upper(&self) -> f64 {
        self.upper
    }

    pub(crate) const fn centre(&self) -> f64 {
        (self.upper + self.lower) * 0.5
    }

    pub(crate) const fn width(&self) -> f64 {
        self.upper - self.lower
    }

    pub(crate) const fn half_width(&self) -> f64 {
        self.width() * 0.5
    }

    pub(crate) const fn bisect(&self) -> [Self; 2] {
        let upper = self.upper();
        let lower = self.lower();
        let midpoint = (upper + lower) * 0.5;

        [Self::new(lower, midpoint), Self::new(midpoint, upper)]
    }

    #[inline]
    pub(crate) fn subinterval_too_small(&self) -> bool {
        let lower = self.lower();
        let upper = self.upper();
        let midpoint = self.centre();

        let eps = f64::EPSILON;
        let min = f64::MIN_POSITIVE;

        let tmp = (1.0 + 100.0 * eps) * (midpoint.abs() + 1000.0 * min);

        lower.abs() <= tmp && upper.abs() <= tmp
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bisection() {
        let limit = Limits::new(0.0, 1.0);
        let [lower, upper] = limit.bisect();

        assert_eq!(lower, Limits::new(0.0, 0.5));
        assert_eq!(upper, Limits::new(0.5, 1.0));
    }
}
