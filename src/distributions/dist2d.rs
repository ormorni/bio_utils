use itertools::{iproduct, Itertools};
use crate::distributions::convolve_2d;

/// A struct holding the probabilities of a random variable over the non-negative reals.
#[derive(Clone)]
pub struct Dist2D {
    /// The probability distribution.
    pub data: Vec<Vec<f64>>,
    /// The maximal value in a cell before it is dropped.
    drop: f64,
    /// The resolution in which the distribution is kept.
    scale: (f64, f64),
    /// The shift to the data in the ProbArray.
    shift: (isize, isize),
}

impl Dist2D {
    /// Initializes a new empty Dist2D.
    pub fn new(drop: f64, scale: (f64, f64)) -> Dist2D {
        Dist2D {
            data: Vec::new(),
            drop,
            scale,
            shift: (0, 0),
        }
    }

    /// Maps a float to an index of the array.
    fn map_index(&self, idx: (f64, f64)) -> (isize, isize) {
        (
            (idx.0 * self.scale.0).round() as isize + self.shift.0,
            (idx.1 * self.scale.1).round() as isize + self.shift.1,
        )
    }

    /// Initializes a ProbArray from an f64.
    pub fn from_f64(val: (f64, f64), drop: f64, scale: (f64, f64)) -> Dist2D {
        let mut res = Dist2D::new(drop, scale);

        let target_ind = (val.0 * res.scale.0, val.1 * res.scale.1);

        // The index is mapped to the average of two values.
        let higher_frac = (
            target_ind.0 - target_ind.0.floor(),
            target_ind.1 - target_ind.1.floor(),
        );
        let lower_frac = (1. - higher_frac.0, 1. - higher_frac.1);

        res.data = match (higher_frac.0 == 0., higher_frac.1 == 0.) {
            (true, true) => vec![vec![1.]],
            (false, true) => vec![vec![lower_frac.0], vec![higher_frac.0]],
            (true, false) => vec![vec![lower_frac.1, higher_frac.1]],
            (false, false) => vec![
                vec![lower_frac.0 * lower_frac.1, lower_frac.0 * higher_frac.1],
                vec![higher_frac.0 * lower_frac.1, higher_frac.0 * higher_frac.1],
            ],
        };

        res.shift = (
            -target_ind.0.floor() as isize,
            -target_ind.1.floor() as isize,
        );

        res
    }

    /// Trims empty cells from the ProbArray, while keeping the ProbArray in a rectangular shape.
    pub fn trim(&mut self) {
        let mut start = (0, 0);
        let mut end = self.shape();

        while (start.0 < end.0) && (start.1 < end.1) {
            // Testing removing the first row.
            while self.data[start.0][start.1..end.1].iter().sum::<f64>() < self.drop {
                start.0 += 1;
                continue;
            }
            // Testing removing the last row.
            if self.data[end.0 - 1].iter().sum::<f64>() < self.drop {
                end.0 -= 1;
                continue;
            }
            // Testing removing the first column.
            if self.data[start.0..end.0]
                .iter()
                .map(|row| row[start.1])
                .sum::<f64>()
                < self.drop
            {
                start.1 += 1;
                continue;
            }
            // Testing removing the last column.
            if self.data[start.0..end.0]
                .iter()
                .map(|row| row[end.1 - 1])
                .sum::<f64>()
                < self.drop
            {
                end.1 -= 1;
                continue;
            }
            break;
        }

        if (start.0 == end.0) || (start.1 == end.1) {
            self.data = vec![];
            self.shift = (0, 0);
        } else {
            self.data = self.data[start.0..end.0]
                .iter()
                .map(|v| v[start.1..end.1].to_vec())
                .collect_vec();

            self.shift.0 -= start.0 as isize;
            self.shift.1 -= start.1 as isize;
        }
    }

    /// Returns the distribution corresponding to the sum of two independent distributions.
    pub fn convolve(&self, other: &Dist2D) -> Dist2D {
        for (i, j) in iproduct!(0..self.data.len(), 0..self.data[0].len()) {
            assert!(self.data[i][j].is_finite());
        }
        for (i, j) in iproduct!(0..other.data.len(), 0..other.data[0].len()) {
            assert!(other.data[i][j].is_finite());
        }

        let mut res = Dist2D::new(self.drop, self.scale);
        res.data = convolve_2d(self.data.clone(), other.data.clone());
        res.shift = (self.shift.0 + other.shift.0, self.shift.1 + other.shift.1);
        for (i, j) in iproduct!(0..res.data.len(), 0..res.data[0].len()) {
            assert!(res.data[i][j].is_finite());
        }

        res.trim();
        res
    }

    /// The length of the internal array.
    pub fn shape(&self) -> (usize, usize) {
        if self.data.len() == 0 {
            return (0, 0);
        }
        (self.data.len(), self.data[0].len())
    }

    /// Returns the distribution matching the weighted average of the two distributions.
    pub fn weighted_average(dist_1: &Dist2D, p1: f64, dist_2: &Dist2D, p2: f64) -> Dist2D {
        assert!(p1.is_finite());
        assert!(p2.is_finite());

        for (i, j) in iproduct!(0..dist_1.data.len(), 0..dist_1.data[0].len()) {
            assert!(dist_1.data[i][j].is_finite());
        }
        for (i, j) in iproduct!(0..dist_2.data.len(), 0..dist_2.data[0].len()) {
            assert!(dist_2.data[i][j].is_finite());
        }

        let mut res = Dist2D::new(dist_1.drop, dist_1.scale);
        res.shift = (
            dist_1.shift.0.max(dist_2.shift.0),
            dist_1.shift.1.max(dist_2.shift.1),
        );

        let delta_1 = (res.shift.0 - dist_1.shift.0, res.shift.1 - dist_1.shift.1);
        let delta_2 = (res.shift.0 - dist_2.shift.0, res.shift.1 - dist_2.shift.1);

        let shape = (
            (dist_1.shape().0 as isize + delta_1.0).max(dist_2.shape().0 as isize + delta_2.0)
                as usize,
            (dist_1.shape().1 as isize + delta_1.1).max(dist_2.shape().1 as isize + delta_2.1)
                as usize,
        );

        res.data = vec![vec![0.; shape.1]; shape.0];

        for (i, j) in iproduct!(0..dist_1.shape().0, 0..dist_1.shape().1) {
            res.data[i + delta_1.0 as usize][j + delta_1.1 as usize] += p1 * dist_1.data[i][j];
        }
        for (i, j) in iproduct!(0..dist_2.shape().0, 0..dist_2.shape().1) {
            res.data[i + delta_2.0 as usize][j + delta_2.1 as usize] += p2 * dist_2.data[i][j];
        }

        for (i, j) in iproduct!(0..res.data.len(), 0..res.data[0].len()) {
            assert!(res.data[i][j].is_finite());
        }

        res.trim();



        res
    }

    /// Returns an integral over the function using distribution as a measure.
    pub fn integrate(&self, function: impl Fn((f64, f64)) -> f64) -> f64 {
        let mut res = 0.;
        for idx in iproduct!(0..self.shape().0, 0..self.shape().1) {
            let unmapped = (
                (idx.0 as isize - self.shift.0) as f64 / self.scale.0,
                (idx.1 as isize - self.shift.1) as f64 / self.scale.1,
            );
            res += function(unmapped) * self.data[idx.0][idx.1];
        }
        res
    }
}

#[cfg(test)]
mod tests_2d {
    use std::cmp::Ordering;
    use std::f64::consts::PI;
    use crate::distributions::dist2d::Dist2D;
    use crate::test_utils::assert_close;

    #[test]
    fn test_conv() {
        let f1 = Dist2D::from_f64((0.5, 0.5), 0., (1., 1.));
        let f2 = Dist2D::from_f64((-0.5, -0.5), 0., (1., 1.));
        let conv = f1.convolve(&f2);

        assert_eq!(conv.shift, (1, 1));
        assert_eq!(
            conv.data,
            vec![
                vec![0.0625, 0.125, 0.0625],
                vec![0.1250, 0.250, 0.1250],
                vec![0.0625, 0.125, 0.0625]
            ]
        );
    }

    #[test]
    fn test_av() {
        let f1 = Dist2D::from_f64((0.5, 0.5), 0., (1., 1.));
        let f2 = Dist2D::from_f64((-0.5, -0.5), 0., (1., 1.));
        let av = Dist2D::weighted_average(&f1, 0.5, &f2, 0.5);

        assert_eq!(av.shift, (1, 1));
        assert_eq!(
            av.data,
            vec![
                vec![0.125, 0.125, 0.000],
                vec![0.125, 0.250, 0.125],
                vec![0.000, 0.125, 0.125]
            ]
        );
    }

    #[test]
    fn test_integrate() {
        const CONV_ROUNDS: i32 = 12;
        const EXPECTATION: i32 = 1 << (CONV_ROUNDS - 1);

        let mut f = Dist2D::from_f64((0.5, 0.5), 1e-10, (1., 1.));
        for _ in 0..CONV_ROUNDS {
            f = f.convolve(&f);
        }

        assert_close(f.integrate(|_| 1.), 1.);
        let one_sixth_normal = f.integrate(|(x, y)| {
            match ((y - EXPECTATION as f64) / (x - EXPECTATION as f64))
                .atan()
                .total_cmp(&(PI / 3.))
            {
                Ordering::Less => 0.,
                Ordering::Equal => 0.5,
                Ordering::Greater => 1.,
            }
        });
        assert!((one_sixth_normal - 1. / 6.).abs() < 1e-3);
    }
}
