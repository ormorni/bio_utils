use itertools::Itertools;
use rustfft::num_traits::MulAddAssign;

use crate::distributions::convolve;
use std::{
    fmt::{Debug, Formatter},
    ops::{Add, AddAssign, Mul, MulAssign},
};

/// A struct holding the probabilities of a random variable over the non-negative reals.
#[derive(Clone)]
pub struct Dist1D {
    /// The probability distribution.
    pub data: Vec<f64>,
    /// The maximal value in a cell before it is dropped.
    pub drop: f64,
    /// The resolution in which the distribution is kept.
    pub scale: f64,
    /// A shift, used to include negative numbers.
    pub shift: isize,
    /// Tracks the meaningful length of the data in the distribution.
    length: usize,
}

impl Dist1D {
    pub fn new(drop: f64, scale: f64) -> Dist1D {
        Dist1D {
            data: Vec::new(),
            drop,
            scale,
            shift: 0,
            length: 0,
        }
    }

    /// Maps a float to an index of the array.
    pub fn map_index(&self, idx: f64) -> isize {
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
        let higher_part = target_ind - target_ind.floor();
        let lower_part = 1. - higher_part;

        if higher_part != 0. {
            res.data = vec![lower_part, higher_part];
            res.length = 2;
        } else {
            res.data = vec![lower_part];
            res.length = 1;
        }
        res.shift = -(target_ind.floor() as isize);
        assert!(res.length <= res.data.len());
        res
    }

    pub fn from_vec(data: &Vec<f64>, shift: isize, drop: f64, scale: f64) -> Dist1D {
        let mut res_data = Vec::with_capacity(data.len().next_power_of_two());
        res_data.extend(data);
        res_data.resize(res_data.len().next_power_of_two(), 0.);

        Dist1D {
            data: res_data,
            drop,
            scale,
            shift,
            length: data.len(),
        }
    }

    /// Trims empty cells from the distribution.
    /// The distribution is only trimmed to lengths that are powers of 2.
    pub fn trim(&mut self) {
        let mut new_start = 0;
        while new_start < self.data.len() && self.data[new_start].abs() < self.drop {
            new_start += 1;
        }
        let mut new_end = self.length;
        while new_end > 0 && self.data[new_end - 1].abs() < self.drop {
            new_end -= 1;
        }
        if new_start >= new_end {
            self.data = vec![];
            self.shift = 0;
            self.length = 0;
            return;
        }
        let new_len = new_end - new_start;
        if new_len < self.length {
            for i in 0..new_len {
                self.data[i] = self.data[i + new_start];
            }
            self.shift -= new_start as isize;
        }
        self.length = new_end - new_start;
        self.data.resize(self.length.next_power_of_two(), 0.);
        assert!(self.length <= self.data.len());
    }

    /// Returns the distribution corresponding to the sum of two independent distributions.
    pub fn convolve(&self, other: &Dist1D) -> Dist1D {
        let mut res = Dist1D::new(self.drop, self.scale);
        res.data = convolve(&self.data, &other.data);
        for i in 0..res.len() {
            res.data[i] = res.data[i].max(0.);
        }
        res.shift = self.shift + other.shift;
        res.length = res.data.len();
        res.trim();
        res
    }

    /// The length of the internal array.
    pub fn len(&self) -> usize {
        self.length
    }

    /// Returns the distribution matching the weighted average of the two distributions.
    pub fn weighted_average(dist_1: &Dist1D, p1: f64, dist_2: &Dist1D, p2: f64) -> Dist1D {
        &(dist_1 * p1) + &(dist_2 * p2)
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

    pub fn add_mul_assign(&mut self, rhs: &Dist1D, factor: f64) {
        let low = (-self.shift).min(-rhs.shift);
        let high = (self.len() as isize - self.shift).max(rhs.len() as isize - rhs.shift);
        if high - low > self.data.len() as isize {
            // The current allocated array will not suffice.
            self.data
                .resize(((high - low) as usize).next_power_of_two(), 0.);
        }
        // Shifting the internal data.
        if rhs.shift > self.shift {
            let diff = (rhs.shift - self.shift) as usize;
            for i in (diff..self.data.len()).rev() {
                self.data[i] = self.data[i - diff];
            }
            self.data[..diff].fill(0.);
            self.shift = rhs.shift;
        }
        let diff = (self.shift - rhs.shift) as usize;
        for i in 0..rhs.len() {
            self.data[i + diff] += factor * rhs.data[i];
        }
        self.length = (high - low) as usize;
        assert!(self.length <= self.data.len());
    }
}

impl Mul<f64> for &Dist1D {
    type Output = Dist1D;

    fn mul(self, rhs: f64) -> Self::Output {
        let data = self.data.iter().map(|i| i * rhs).collect_vec();
        Dist1D {
            data,
            drop: self.drop,
            scale: self.scale,
            shift: self.shift,
            length: self.length,
        }
    }
}

impl MulAssign<f64> for Dist1D {
    fn mul_assign(&mut self, rhs: f64) {
        for i in self.data.iter_mut() {
            *i *= rhs;
        }
    }
}

impl Add for &Dist1D {
    type Output = Dist1D;

    fn add(self, rhs: Self) -> Self::Output {
        let mut res = Dist1D::new(self.drop, self.scale);
        res.shift = self.shift.max(rhs.shift);
        let res_len = ((self.len() as isize - self.shift).max(rhs.len() as isize - rhs.shift)
            + res.shift) as usize;
        let alloc_len = res_len.next_power_of_two();

        let delta_1 = res.shift - self.shift;
        let delta_2 = res.shift - rhs.shift;
        res.data = vec![0.; alloc_len];
        for i in 0..self.len() {
            res.data[(i as isize + delta_1) as usize] += self.data[i];
        }
        for i in 0..rhs.len() {
            res.data[(i as isize + delta_2) as usize] += rhs.data[i];
        }
        res.length = res_len;
        assert!(res.length <= res.data.len());
        res
    }
}

impl AddAssign<&Dist1D> for Dist1D {
    fn add_assign(&mut self, rhs: &Self) {
        let low = (-self.shift).min(-rhs.shift);
        let high = (self.len() as isize - self.shift).max(rhs.len() as isize - rhs.shift);
        if high - low > self.data.len() as isize {
            // The current allocated array will not suffice.
            self.data
                .resize(((high - low) as usize).next_power_of_two(), 0.);
        }
        // Shifting the internal data.
        if rhs.shift > self.shift {
            let diff = (rhs.shift - self.shift) as usize;
            for i in (diff..self.data.len()).rev() {
                self.data[i] = self.data[i - diff];
            }
            self.data[..diff].fill(0.);
            self.shift = rhs.shift;
        }
        let diff = (self.shift - rhs.shift) as usize;
        for i in 0..rhs.len() {
            self.data[i + diff] += rhs.data[i];
        }
        self.length = (high - low) as usize;
        assert!(self.length <= self.data.len());
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
    use crate::utils::tests::{assert_close, SEED};
    use itertools::Itertools;
    use rand::distributions::Uniform;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};
    use std::cmp::Ordering;

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
        assert_eq!(f1.data[..f1.len()], vec![1.]);
        let f2 = Dist1D::from_f64(0.25, 0.01, 1.);
        assert_eq!(f2.shift, 0);
        assert_eq!(f2.data[..f2.len()], vec![0.75, 0.25]);

        let av = Dist1D::weighted_average(&f1, 0.25, &f2, 0.75);
        assert_eq!(av.data[..av.len()], vec![0.5625, 0.1875, 0.25]);
        assert_eq!(av.shift, 0);
    }

    #[test]
    /// Tests a single example of trimming.
    fn trim_example() {
        let mut dist = Dist1D::new(0.01, 1.);
        dist.data = vec![0.0, 0.0, 1., 0.0, 0.0];
        dist.length = 5;
        dist.trim();
        // assert_close(dist.sum_range(1.5, 2.5), 1.);
        assert_close(
            dist.integrate(|f| if (1.5..2.5).contains(&f) { 1. } else { 0. }),
            1.,
        );
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

    #[test]
    fn test_ops() {
        let mut d1 = Dist1D {
            data: vec![1., 1., 2., 3.],
            drop: 1e-6,
            scale: 1.,
            shift: 0,
            length: 4,
        };
        assert_eq!((&d1 * 2.).data, vec![2., 2., 4., 6.]);

        let d2 = Dist1D {
            data: vec![1., 1., 2., 3.],
            drop: 1e-6,
            scale: 1.,
            shift: 1,
            length: 4,
        };

        let mut d3 = Dist1D {
            data: vec![2., 1., 2., 3.],
            drop: 1e-6,
            scale: 1.,
            shift: -1,
            length: 3,
        };

        assert_eq!(
            (&d1 + &d2).data,
            vec![1.0, 2.0, 3.0, 5.0, 3.0, 0.0, 0.0, 0.0]
        );

        assert_eq!((&d1 + &d3).data, vec![1.0, 3.0, 3.0, 5.0]);

        d1 += &d3;
        assert_eq!(d1.data, vec![1.0, 3.0, 3.0, 5.0]);

        d3 += &d1;
        assert_eq!(d3.data, vec![1.0, 5.0, 4.0, 7.0]);
    }

    #[test]
    fn test_add() {
        let mut rng = StdRng::seed_from_u64(0x12345678);

        for _ in 0..100 {
            let v1 = (&mut rng)
                .sample_iter(rand_distr::Normal::new(0., 1.).unwrap())
                .take(100)
                .collect::<Vec<f64>>();
            let v2 = (&mut rng)
                .sample_iter(rand_distr::Normal::new(0., 1.).unwrap())
                .take(100)
                .collect::<Vec<f64>>();

            let shift_1 = rng.gen_range(-10..10);
            let shift_2 = rng.gen_range(-10..10);
            let drop = 1e-6;
            let scale = rng.gen::<f64>() + 1.;
            let d1 = Dist1D::from_vec(&v1, shift_1, drop, scale);
            let d2 = Dist1D::from_vec(&v2, shift_2, drop, scale);

            let a = rng.gen::<f64>();
            let b = rng.gen::<f64>();
            let c = rng.gen::<f64>();
            let d = rng.gen::<f64>();

            let d3 = &d1 + &d2;

            assert_close(
                d1.integrate(|f| a * f.powi(3) + b * f.powi(2) + c * f + d)
                    + d2.integrate(|f| a * f.powi(3) + b * f.powi(2) + c * f + d),
                d3.integrate(|f| a * f.powi(3) + b * f.powi(2) + c * f + d),
            )
        }
    }

    #[test]
    fn test_add_assign() {
        let mut rng = StdRng::seed_from_u64(0x12345678);
        let v1 = (&mut rng)
            .sample_iter(rand_distr::Normal::new(0., 1.).unwrap())
            .take(100)
            .collect::<Vec<f64>>();
        let v2 = (&mut rng)
            .sample_iter(rand_distr::Normal::new(0., 1.).unwrap())
            .take(100)
            .collect::<Vec<f64>>();

        let shift_1 = rng.gen_range(-10..10);
        let shift_2 = rng.gen_range(-10..10);
        let drop = 1e-6;
        let scale = rng.gen::<f64>() + 1.;
        let mut d1 = Dist1D::from_vec(&v1, shift_1, drop, scale);
        let d2 = Dist1D::from_vec(&v2, shift_2, drop, scale);

        let a = rng.gen::<f64>();
        let b = rng.gen::<f64>();
        let c = rng.gen::<f64>();
        let d = rng.gen::<f64>();

        let s1 = d1.integrate(|f| a * f.powi(3) + b * f.powi(2) + c * f + d)
            + d2.integrate(|f| a * f.powi(3) + b * f.powi(2) + c * f + d);

        d1 += &d2;

        assert_close(
            s1,
            d1.integrate(|f| a * f.powi(3) + b * f.powi(2) + c * f + d),
        )
    }
}
