use itertools::{iproduct, Itertools};
use realfft::num_complex::Complex;
use realfft::num_traits::{Float, FromPrimitive, Zero};
use realfft::{FftNum, RealFftPlanner};
use std::mem::size_of;

/// Computes the smallest power of 2^n greater or equal to the length.
fn fft_array_size(length: usize) -> usize {
    if length == 0 {
        0
    } else {
        1 << (size_of::<usize>() * 8 - (length - 1).leading_zeros() as usize)
    }
}

/// Transposes a 2D array.
fn transposed<T: Clone>(arr: &Vec<Vec<T>>) -> Vec<Vec<T>> {
    assert!(arr.len() > 0);
    assert!(arr[0].len() > 0);

    let mut res = vec![vec![arr[0][0].clone(); arr.len()]; arr[0].len()];
    for (i, j) in iproduct!(0..arr.len(), 0..arr[0].len()) {
        res[j][i] = arr[i][j].clone();
    }
    res
}

/// Convolves two 1D arrays naively.
#[allow(unused)]
pub fn convolve_naive<T: Float>(v1: &[T], v2: &[T]) -> Vec<T> {
    let mut res = vec![T::zero(); v1.len() + v2.len() - 1];

    for (i, j) in iproduct!(0..v1.len(), 0..v2.len()) {
        res[i + j] = res[i + j] + v1[i] * v2[j];
    }

    res
}

pub fn convolve<T: FftNum>(v1: &[T], v2: &[T]) -> Vec<T> {
    let res_length = v1.len() + v2.len() - 1;
    let alloc_length = fft_array_size(res_length);

    let mut planner = realfft::RealFftPlanner::new();
    let forward = planner.plan_fft_forward(alloc_length);
    let backward = planner.plan_fft_inverse(alloc_length);

    let mut forward_inp = Vec::from(v1);
    forward_inp.resize(alloc_length, T::zero());
    let mut backward_inp = backward.make_input_vec();

    forward
        .process(&mut forward_inp, &mut backward_inp)
        .unwrap();

    let mut forward_inp = Vec::from(v2);
    forward_inp.resize(alloc_length, T::zero());

    let mut forward_out = forward.make_output_vec();
    forward.process(&mut forward_inp, &mut forward_out).unwrap();

    for i in 0..forward_out.len() {
        backward_inp[i] = backward_inp[i] * forward_out[i] / T::from_usize(alloc_length).unwrap();
    }
    let length = backward_inp.len();
    backward_inp[0].im = T::zero();
    backward_inp[length - 1].im = T::zero();
    let mut backward_out = backward.make_output_vec();
    backward
        .process(&mut backward_inp, &mut backward_out)
        .expect(&*format!(
            "Backward FFT failed on array of length {} -> {}",
            backward_inp.len(),
            backward_out.len()
        ));

    backward_out.resize(res_length, T::zero());
    backward_out
}

/// Convolves two 2D arrays naively.
#[allow(unused)]
pub fn convolve_2d_naive(v1: &Vec<Vec<f64>>, v2: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let x1 = v1.len();
    let y1 = v1[0].len();

    let x2 = v2.len();
    let y2 = v2[0].len();

    let mut res = vec![vec![0.; y1 + y2 - 1]; x1 + x2 - 1];

    for (i1, j1, i2, j2) in iproduct!(0..x1, 0..y1, 0..x2, 0..y2) {
        res[i1 + i2][j1 + j2] += v1[i1][j1] * v2[i2][j2];
    }

    res
}

/// Convolves two 2D arrays naively.
pub fn convolve_2d<T: FftNum>(mut v1: Vec<Vec<T>>, mut v2: Vec<Vec<T>>) -> Vec<Vec<T>> {
    let x1 = v1.len();
    let y1 = v1[0].len();

    let x2 = v2.len();
    let y2 = v2[0].len();

    let xr = x1 + x2 - 1;
    let yr = y1 + y2 - 1;

    let xa = fft_array_size(xr);
    let ya = fft_array_size(yr);

    let mut r2c_planner = RealFftPlanner::new();
    let mut c2c_planner = rustfft::FftPlanner::new();

    // Performing FFT on the rows.
    let forward = r2c_planner.plan_fft_forward(ya);

    let v1 = v1
        .iter_mut()
        .map(|v| {
            v.resize(ya, T::zero());
            let mut res = forward.make_output_vec();
            forward.process(v, &mut res).unwrap();
            res
        })
        .collect_vec();

    let v2 = v2
        .iter_mut()
        .map(|v| {
            v.resize(ya, T::zero());
            let mut res = forward.make_output_vec();
            forward.process(v, &mut res).unwrap();
            res
        })
        .collect_vec();

    // Transposing.
    let mut v1 = transposed(&v1);
    let mut v2 = transposed(&v2);

    // Performing FFT on the columns.
    let forward = c2c_planner.plan_fft_forward(xa);

    for v in v1.iter_mut() {
        v.resize(xa, Complex::zero());
        forward.process(v);
    }

    for v in v2.iter_mut() {
        v.resize(xa, Complex::zero());
        forward.process(v);
    }

    // Multiplying pointwise.
    for (i, j) in iproduct!(0..v1.len(), 0..v1[0].len()) {
        v1[i][j] = v1[i][j] * v2[i][j] / Complex::from_usize(xa * ya).unwrap();
    }

    // Performing iFFT on the columns.
    let backward = c2c_planner.plan_fft_inverse(xa);
    for v in v1.iter_mut() {
        backward.process(v);
        v.resize(xr, Complex::zero());
    }

    // Transposing and performing FFT on the rows.
    let backward = r2c_planner.plan_fft_inverse(ya);
    transposed(&v1)
        .iter_mut()
        .map(|v| {
            let mut u = backward.make_output_vec();
            // Manually setting positions to 0 to satisfy RealFFT.
            let last = v.len() - 1;
            v[0].im = T::zero();
            v[last].im = T::zero();

            backward.process(v, &mut u).unwrap();
            u.resize(yr, T::zero());
            u
        })
        .collect_vec()
}

#[cfg(test)]
mod tests {
    use crate::distributions::{convolve, convolve_2d, convolve_2d_naive, convolve_naive};
    use itertools::{izip, Itertools};
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};
    use crate::test_utils::{assert_close, SEED};

    /// Tests that the FFT-based convolution and the regular convolution return the same results.
    #[test]
    fn test_conv() {
        let mut rng = StdRng::from_seed(SEED);
        for _ in 0..100 {
            let x1 = rng.gen_range(1..100);
            let x2 = rng.gen_range(1..100);

            let v1 = (0..x1).map(|_| rng.gen::<f64>()).collect_vec();
            let v2 = (0..x2).map(|_| rng.gen::<f64>()).collect_vec();

            let conv1 = convolve_naive(&v1, &v2);
            let conv2 = convolve(&v1, &v2);

            assert_eq!(conv1.len(), conv2.len());
            for (i, j) in izip!(conv1, conv2) {
                assert_close(i, j);
            }
        }
    }

    /// Tests that the FFT-based 2D convolution and the regular 2D convolution return the same results.
    #[test]
    fn test_conv_2d() {
        let mut rng = StdRng::from_seed(SEED);
        for _ in 0..100 {
            let x1: usize = rng.gen_range(1..33);
            let y1: usize = rng.gen_range(1..33);
            let x2: usize = rng.gen_range(1..33);
            let y2: usize = rng.gen_range(1..33);

            let v1 = (0..y1)
                .map(|_| (0..x1).map(|_| rng.gen::<f64>()).collect_vec())
                .collect_vec();
            let v2 = (0..y2)
                .map(|_| (0..x2).map(|_| rng.gen::<f64>()).collect_vec())
                .collect_vec();

            let conv1 = convolve_2d_naive(&v1, &v2);
            let conv2 = convolve_2d(v1, v2);

            assert_eq!(conv1.len(), conv2.len());
            assert_eq!(conv1[0].len(), conv2[0].len());
            for (i, j) in izip!(conv1, conv2) {
                for (i, j) in izip!(i, j) {
                    assert_close(i, j);
                }
            }
        }
    }
}
