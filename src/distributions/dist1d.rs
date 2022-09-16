use itertools::{iproduct, Itertools};
use std::fmt::{Debug, Formatter};
use crate::distributions::{convolve, convolve_2d};

/// A struct holding the probabilities of a random variable over the non-negative reals.
#[derive(Clone)]
pub struct Dist1D {
    /// The probability distribution.
    pub(crate) data: Vec<f64>,
    /// The maximal value in a cell before it is dropped.
    drop: f64,
    /// The resolution in which the distribution is kept.
    scale: f64,
    /// A shift, used to include negative numbers.
    shift: isize,
}

impl Dist1D {
    pub fn new(drop: f64, scale: f64) -> Dist1D {
        Dist1D {
            data: Vec::new(),
            drop,
            scale,
            shift: 0,
        }
    }

    /// Maps a float to an index of the array.
    pub(crate) fn map_index(&self, idx: f64) -> isize {
        // Infinite indices are passed as-is.
        if idx.is_finite() {
            return (idx * self.scale).round() as isize + self.shift;
        }
        if idx.is_infinite() && idx > 0. {
            return isize::MAX;
        }
        if idx.is_infinite() && idx < 0. {
            return isize::MIN;
        }
        panic!("Mapped illegal index: {}", idx)
    }

    /// Initializes a ProbArray from an f64.
    pub fn from_f64(val: f64, drop: f64, scale: f64) -> Dist1D {
        let mut res = Dist1D::new(drop, scale);

        let target_ind = val * res.scale;

        // The index is mapped to the average of two values.
        let higher_frac = target_ind - target_ind.floor();
        let lower_frac = 1. - higher_frac;

        if higher_frac != 0. {
            res.data = vec![lower_frac, higher_frac];
        } else {
            res.data = vec![lower_frac];
        }
        res.shift = -(target_ind.floor() as isize);

        res
    }

    /// Trims empty cells from the ProbArray.
    pub fn trim(&mut self) {
        let mut new_start = 0;
        while new_start < self.data.len() && self.data[new_start].abs() < self.drop {
            new_start += 1;
        }
        let mut new_end = self.data.len();
        while new_end > 0 && self.data[new_end - 1].abs() < self.drop {
            new_end -= 1;
        }
        if new_start < new_end {
            self.data = self.data[new_start..new_end].to_vec();
            self.shift -= new_start as isize;
        } else {
            self.data = vec![0.];
            self.shift = 0;
        }
    }

    /// Returns the distribution corresponding to the sum of two independent distributions.
    pub fn convolve(&self, other: &Dist1D) -> Dist1D {
        let mut res = Dist1D::new(self.drop, self.scale);
        res.data = convolve(&self.data, &other.data);
        for i in 0..res.len() {
            res.data[i] = res.data[i].max(0.);
        }
        res.shift = self.shift + other.shift;
        res.trim();
        res
    }

    /// The length of the internal array.
    pub(crate) fn len(&self) -> usize {
        self.data.len()
    }

    /// Returns the distribution matching the weighted average of the two distributions.
    pub fn weighted_average(dist_1: &Dist1D, p1: f64, dist_2: &Dist1D, p2: f64) -> Dist1D {
        assert!(p1.is_finite());
        assert!(p2.is_finite());

        let mut res = Dist1D::new(dist_1.drop, dist_1.scale);
        res.shift = dist_1.shift.max(dist_2.shift);
        let delta_1 = res.shift - dist_1.shift;
        let delta_2 = res.shift - dist_2.shift;
        res.data = vec![
            0.;
            (dist_1.len() as isize + delta_1).max(dist_2.len() as isize + delta_2)
                as usize
        ];
        for i in 0..dist_1.len() {
            res.data[(i as isize + delta_1) as usize] += p1 * dist_1.data[i];
        }
        for i in 0..dist_2.len() {
            res.data[(i as isize + delta_2) as usize] += p2 * dist_2.data[i];
        }

        res.trim();
        res
    }

    /// Returns an integral over the function using distribution as a measure.
    pub fn integrate(&self, function: impl Fn(f64) -> f64) -> f64 {
        let mut res = 0.;
        for (idx, prob) in self.data.iter().enumerate() {
            let unmapped = (idx as isize - self.shift) as f64 / self.scale;
            res += function(unmapped) * prob;
        }
        res
    }
}

impl Debug for Dist1D {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str(
            format!(
                "Dist1D(scale={}, shift={}, len={})",
                self.scale,
                self.shift,
                self.len()
            )
            .as_str(),
        )
    }
}


#[cfg(test)]
mod tests_1d {
    use crate::distributions::dist1d::Dist1D;
    use itertools::Itertools;
    use rand::distributions::Uniform;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};
    use std::cmp::Ordering;
    use crate::test_utils::{assert_close, SEED};

    #[test]
    fn test_conv() {
        let f1 = Dist1D::from_f64(1.5, 0., 6.);
        let f2 = Dist1D::from_f64(1.5, 0., 6.);
        let conv = f1.convolve(&f2);
        assert_eq!(conv.map_index(3.), 0);
        assert_eq!(conv.data, vec![1.]);
    }

    #[allow(unused)]
    fn test_conv_large() {
        let mut f1 = Dist1D::new(0., 10.);
        let mut f2 = Dist1D::new(0., 10.);

        f1.shift = 10;
        f1.data = vec![1., 1., 1.];
        f2.data = vec![1., 1., 1.];

        let conv = f1.convolve(&f2);
        println!("{:?}", conv.data);
        println!("{:?}", conv.shift);
        println!(
            "{:?}",
            (-10..60)
                .map(|i| conv.map_index(i as f64 * 0.1))
                .collect_vec()
        );
    }

    #[test]
    fn test_av() {
        let f1 = Dist1D::from_f64(2., 0.01, 1.);
        assert_eq!(f1.shift, -2);
        assert_eq!(f1.data, vec![1.]);
        let f2 = Dist1D::from_f64(0.25, 0.01, 1.);
        assert_eq!(f2.shift, 0);
        assert_eq!(f2.data, vec![0.75, 0.25]);

        let av = Dist1D::weighted_average(&f1, 0.25, &f2, 0.75);
        assert_eq!(av.data, vec![0.5625, 0.1875, 0.25]);
        assert_eq!(av.shift, 0);
    }

    #[test]
    /// Tests a single example of trimming.
    fn trim_example() {
        let mut dist = Dist1D::new(0.01, 1.);
        dist.data = vec![0.0, 0.0, 1., 0.0, 0.0];
        dist.trim();
        // assert_close(dist.sum_range(1.5, 2.5), 1.);
        assert_close(dist.integrate(|f|if (1.5..2.5).contains(&f) {1.} else {0.}), 1.);
        assert_eq!(dist.data.len(), 1);
        assert_eq!(dist.shift, -2);
        assert_eq!(dist.map_index(1.9), 0);
    }

    /// Creates a random distribution.
    fn random_dist(rng: &mut StdRng) -> Dist1D {
        let mut res = Dist1D::new(0., rng.gen::<f64>());
        res.data = rng
            .sample_iter(Uniform::new(-1., 1.))
            .take(10)
            .collect_vec();
        res
    }

    #[test]
    /// Tests multiple randomized trims.
    fn trim_test() {
        let mut rng = StdRng::from_seed(SEED);
        for _ in 0..100 {
            let mut p1 = random_dist(&mut rng);
            let mut p2 = random_dist(&mut rng);
            for i in 0..3 {
                p1.data[i] *= 0.01;
                p1.data[9 - i] *= 0.01;
                p2.data[i] *= 0.01;
                p2.data[9 - i] *= 0.01;
            }

            p1.drop = 0.;
            p2.scale = p1.scale;
            p2.drop = 1e-3;

            let f1 = p1.convolve(&p2);
            let f2 = p2.convolve(&p1);

            // println!("{:?} {:?}", f1, f2);
            // println!("{:?} {:?}", f1.data, f2.data);

            for i in 0..f1.len() {
                let tar_ind = (i as isize + f2.shift - f1.shift) as usize;
                if (0..f2.len()).contains(&tar_ind) {
                    assert_close(f1.data[i], f2.data[tar_ind]);
                }
            }
        }
    }

    #[test]
    fn test_integrate() {
        const CONV_ROUNDS: i32 = 13;
        let expect = 0.5 * 2f64.powi(CONV_ROUNDS);
        let std = (0.25 * 2f64.powi(CONV_ROUNDS)).sqrt();

        let mut f = Dist1D::from_f64(0.5, 0., 3.);
        for _ in 0..CONV_ROUNDS {
            f = f.convolve(&f);
        }

        assert_close(f.integrate(|_| 1.), 1.);
        let norm_cdf_1 = f.integrate(|x| match ((x - expect as f64) / std).total_cmp(&(1.)) {
            Ordering::Less => 1.,
            Ordering::Equal => 0.5,
            Ordering::Greater => 0.,
        });
        assert_close(norm_cdf_1, 0.8413447460685429);
    }
}

