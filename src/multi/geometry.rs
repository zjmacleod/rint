use crate::Limits;

pub(crate) struct Geometry<const NDIM: usize> {
    pub(crate) centre: [f64; NDIM],
    pub(crate) half_widths: [f64; NDIM],
    pub(crate) volume: f64,
    pub(crate) largest_axis: usize,
}

impl<const NDIM: usize> Geometry<NDIM> {
    pub(crate) const fn new(limits: &[Limits; NDIM]) -> Self {
        let mut centre = [0.0; NDIM];
        let mut half_widths = [0.0; NDIM];

        let mut volume = 1.0;

        let mut j = 0;
        let mut largest_axis = 0;
        while j < NDIM {
            centre[j] = limits[j].centre();
            half_widths[j] = limits[j].half_width();
            volume *= half_widths[j];
            if half_widths[j] > half_widths[largest_axis] {
                largest_axis = j;
            };
            j += 1;
        }

        Self {
            centre,
            half_widths,
            volume,
            largest_axis,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn geometry_tests() {
        let limits = [
            Limits::new(0.0, 10.0),
            Limits::new(-1.0, 2.0),
            Limits::new(-100.0, 0.0),
        ];

        let centre_should_be = [5.0, 0.5, -50.0];
        let half_widths_should_be = [5.0, 1.5, 50.0];
        let volume_should_be = 375.0;
        let largest_axis_should_be = 2;

        let Geometry {
            centre,
            half_widths,
            volume,
            largest_axis,
        } = Geometry::new(&limits);

        for (a, b) in centre.iter().zip(centre_should_be.iter()) {
            assert_eq!(a, b);
        }

        for (a, b) in half_widths.iter().zip(half_widths_should_be.iter()) {
            assert_eq!(a, b);
        }

        assert_eq!(volume, volume_should_be);
        assert_eq!(largest_axis, largest_axis_should_be);
    }
}
