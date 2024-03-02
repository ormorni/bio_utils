use itertools::Itertools;
use rustfft::FftNum;

use crate::distributions::convolve;
use std::{
    fmt::{Debug, Display, Formatter},
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign},
};

pub trait Float: num_traits::Float + AddAssign + Display + FftNum {}

impl<T: num_traits::Float + AddAssign + Display + FftNum> Float for T {}
/// A struct holding the probabilities of a random variable over the non-negative reals.
#[derive(Clone)]
pub struct Dist1D<F: Float> {
    /// The probability distribution.
    pub data: Vec<F>,
    /// The maximal value in a cell before it is dropped.
    pub drop: F,
    /// The resolution in which the distribution is kept.
    pub scale: F,
    /// A shift, used to include negative numbers.
    pub shift: isize,
    /// Tracks the meaningful length of the data in the distribution.
    length: usize,
}

impl<F: Float> Dist1D<F> {
    pub fn new(drop: F, scale: F) -> Dist1D<F> {
        Dist1D {
            data: Vec::new(),
            drop,
            scale,
            shift: 0,
            length: 0,
        }
    }

    /// Maps a float to an index of the array.
    pub fn map_index(&self, idx: F) -> isize {
        // Infinite indices are passed as-is.
        if idx.is_finite() {
            return (idx * self.scale).round().to_isize().unwrap() + self.shift;
        }
        if idx.is_infinite() && idx.is_sign_positive() {
            return isize::MAX;
        }
        if idx.is_infinite() && idx.is_sign_negative() {
            return isize::MIN;
        }
        panic!("Mapped illegal index: {}", idx)
    }

    /// Initializes a ProbArray from an F.
    pub fn from_value(val: F, drop: F, scale: F) -> Dist1D<F> {
        let mut res = Dist1D::new(drop, scale);

        let target_ind = val * res.scale;

        // The index is mapped to the average of two values.
        let higher_part = target_ind - target_ind.floor();
        let lower_part = F::one() - higher_part;

        if higher_part.is_zero() {
            res.data = vec![lower_part];
            res.length = 1;
        } else {
            res.data = vec![lower_part, higher_part];
            res.length = 2;
        }
        res.shift = -(target_ind.floor().to_isize().unwrap());
        assert!(res.length <= res.data.len());
        res
    }

    pub fn from_vec(data: &Vec<F>, shift: isize, drop: F, scale: F) -> Dist1D<F> {
        let mut res_data = Vec::with_capacity(data.len().next_power_of_two());
        res_data.extend(data);
        res_data.resize(res_data.len().next_power_of_two(), F::zero());

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
        self.data.resize(self.length.next_power_of_two(), F::zero());
        assert!(self.length <= self.data.len());
    }

    /// Returns the distribution corresponding to the sum of two independent distributions.
    pub fn convolve(&self, other: &Dist1D<F>) -> Dist1D<F> {
        let mut res = Dist1D::new(self.drop, self.scale);
        res.data = convolve(&self.data[..self.len()], &other.data[..other.len()]);
        for i in 0..res.len() {
            res.data[i] = res.data[i].max(F::zero());
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
    pub fn weighted_average(dist_1: &Dist1D<F>, p1: F, dist_2: &Dist1D<F>, p2: F) -> Dist1D<F> {
        &(dist_1 * p1) + &(dist_2 * p2)
    }

    pub fn unmap_index(&self, idx: usize) -> F {
        F::from_isize(idx as isize - self.shift).unwrap() / self.scale
    }

    /// Returns an integral over the function using distribution as a measure.
    pub fn integrate(&self, function: impl Fn(F) -> F) -> F {
        let mut res = F::zero();
        for (idx, prob) in self.data.iter().enumerate() {
            let unmapped = self.unmap_index(idx);
            res += function(unmapped) * *prob;
        }
        res
    }

    pub fn add_mul_assign(&mut self, rhs: &Dist1D<F>, factor: F) {
        let low = (-self.shift).min(-rhs.shift);
        let high = (self.len() as isize - self.shift).max(rhs.len() as isize - rhs.shift);
        if high - low > self.data.len() as isize {
            // The current allocated array will not suffice.
            self.data
                .resize(((high - low) as usize).next_power_of_two(), F::zero());
        }
        // Shifting the internal data.
        if rhs.shift > self.shift {
            let diff = (rhs.shift - self.shift) as usize;
            for i in (diff..self.data.len()).rev() {
                self.data[i] = self.data[i - diff];
            }
            self.data[..diff].fill(F::zero());
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

impl<F: Float> Mul<F> for &Dist1D<F> {
    type Output = Dist1D<F>;

    fn mul(self, rhs: F) -> Self::Output {
        let data = self.data.iter().map(|i| *i * rhs).collect_vec();
        Dist1D::from_vec(&data, self.shift, self.drop, self.scale)
    }
}

impl<F: Float + MulAssign> MulAssign<F> for Dist1D<F> {
    fn mul_assign(&mut self, rhs: F) {
        for i in self.data.iter_mut() {
            *i *= rhs;
        }
    }
}

impl<F: Float> Div<F> for &Dist1D<F> {
    type Output = Dist1D<F>;

    fn div(self, rhs: F) -> Self::Output {
        let data = self.data.iter().map(|i| *i / rhs).collect_vec();
        Dist1D::from_vec(&data, self.shift, self.drop, self.scale)
    }
}

impl<F: Float + DivAssign> DivAssign<F> for Dist1D<F> {
    fn div_assign(&mut self, rhs: F) {
        for i in self.data.iter_mut() {
            *i /= rhs;
        }
    }
}

impl<F: Float + AddAssign> Add for &Dist1D<F> {
    type Output = Dist1D<F>;

    fn add(self, rhs: Self) -> Self::Output {
        let mut res = Dist1D::new(self.drop, self.scale);
        res.shift = self.shift.max(rhs.shift);
        let res_len = ((self.len() as isize - self.shift).max(rhs.len() as isize - rhs.shift)
            + res.shift) as usize;
        let alloc_len = res_len.next_power_of_two();

        let delta_1 = res.shift - self.shift;
        let delta_2 = res.shift - rhs.shift;
        res.data = vec![F::zero(); alloc_len];
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

impl<F: Float + AddAssign> AddAssign<&Dist1D<F>> for Dist1D<F> {
    fn add_assign(&mut self, rhs: &Self) {
        let low = (-self.shift).min(-rhs.shift);
        let high = (self.len() as isize - self.shift).max(rhs.len() as isize - rhs.shift);
        if high - low > self.data.len() as isize {
            // The current allocated array will not suffice.
            self.data
                .resize(((high - low) as usize).next_power_of_two(), F::zero());
        }
        // Shifting the internal data.
        if rhs.shift > self.shift {
            let diff = (rhs.shift - self.shift) as usize;
            for i in (diff..self.data.len()).rev() {
                self.data[i] = self.data[i - diff];
            }
            self.data[..diff].fill(F::zero());
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

impl<F: Float + Display> Debug for Dist1D<F> {
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
        let f1 = Dist1D::from_value(1.5, 0., 6.);
        let f2 = Dist1D::from_value(1.5, 0., 6.);
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
        let f1 = Dist1D::from_value(2., 0.01, 1.);
        assert_eq!(f1.shift, -2);
        assert_eq!(f1.data[..f1.len()], vec![1.]);
        let f2 = Dist1D::from_value(0.25, 0.01, 1.);
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
    fn random_dist(rng: &mut StdRng) -> Dist1D<f64> {
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

        let mut f = Dist1D::from_value(0.5, 0., 3.);
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
