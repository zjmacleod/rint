/// Integration limits for an integration axis.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Limits {
    lower: f64,
    upper: f64,
}

impl Limits {
    /// Generate a new set of integration [`Limits`].
    #[must_use]
    pub const fn new(lower: f64, upper: f64) -> Self {
        Self { lower, upper }
    }

    /// Return the lower integration limit.
    #[must_use]
    pub const fn lower(&self) -> f64 {
        self.lower
    }

    /// Return the upper integration limit.
    #[must_use]
    pub const fn upper(&self) -> f64 {
        self.upper
    }

    /// Return the centre of the integration region.
    pub(crate) const fn centre(&self) -> f64 {
        // Overflow guard
        (self.upper).midpoint(self.lower)
    }

    /// Return the width of the integration region.
    pub(crate) const fn width(&self) -> f64 {
        self.upper - self.lower
    }

    /// Return the half width of the integration region.
    pub(crate) const fn half_width(&self) -> f64 {
        self.width() * 0.5
    }

    /// Bisect the integration region into two new regions.
    #[must_use]
    pub const fn bisect(&self) -> [Self; 2] {
        let upper = self.upper();
        let lower = self.lower();
        let midpoint = (upper + lower) * 0.5;

        [Self::new(lower, midpoint), Self::new(midpoint, upper)]
    }

    /// Determine if a subinterval is too small.
    ///
    /// If an integral has a singularity or other area of difficulty, it is possible for a
    /// subregion to become too small upon further bisection. This function does a guard check of
    /// this condition and returns true iff the subregion can be further subdivided.
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

    /// Scale the limits by a constant factor.
    #[must_use]
    pub fn scale(self, scale: f64) -> Self {
        let lower = self.lower() * scale;
        let upper = self.upper() * scale;

        Limits::new(lower, upper)
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
