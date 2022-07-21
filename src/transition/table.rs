extern crate nalgebra as na;

use super::defs::{Mat21, MAT_SIZE};

use crate::transition::defs::zero_mat;
use std::ops::Mul;

// Defining constants related to the exponentiation process.
// The length of the array of the taylor series factors. Used to reduce the number of matrix multiplications.
const EXP_ARRAY_LENGTH: usize = 10;
// The maximum exponent calculated using the taylor series instead of repeated squaring.
const MAX_EXPONENT: f64 = 0.25;
// The precision used when taking log.
const PREC: f64 = 1e-9;

#[derive(Clone, Debug)]
pub struct TransitionTable {
    /// The table rates of the given table.
    pub rates: Mat21,
    /// An array of matrices used to quickly calculate the exponent of the rates times a constant.
    exp_array: [Mat21; EXP_ARRAY_LENGTH],
    /// A normalization factor used in matrix exponentiation.
    exp_factor: i32,
}

impl TransitionTable {
    pub fn from_rates(rates: &Mat21) -> TransitionTable {
        let mut exp_array = [zero_mat(); EXP_ARRAY_LENGTH];
        exp_array[0].fill_with_identity();

        let mut base_table = rates.clone();
        let mut max_val = base_table
            .iter()
            .map(|x| x.abs())
            .fold(1. as f64, |a, b| a.max(b));
        let mut exp_factor = 0;

        while max_val > 0.1 {
            base_table /= 2.;
            max_val /= 2.;
            exp_factor += 1;
        }

        let mut taylor_factor = 1.;
        for i in 1..EXP_ARRAY_LENGTH {
            taylor_factor *= i as f64;
            exp_array[i] = exp_array[i - 1] * base_table / taylor_factor;
        }

        TransitionTable {
            rates: rates.clone(),
            exp_array,
            exp_factor,
        }
    }

    /// Initializes a table table from the table frequencies.
    pub fn from_freqs(freqs: &Mat21) -> TransitionTable {
        let mut base_freqs = freqs.clone();

        // Subtracting one, to use the taylor approximation of log(x+1).
        for i in 0..MAT_SIZE {
            base_freqs[(i, i)] -= 1.;
        }
        let mut freqs = base_freqs.clone();

        let exp = freqs
            .iter()
            .map(|x| x.abs())
            .fold(0. as f64, |a, b| a.max(b));
        let mut max_val = exp;
        let mut exp_index = 1;

        if max_val > 0.99 {
            panic!("Maximum value of frequency matrix is too large! Log taylor approximation won't converge!")
        }

        let mut rates = zero_mat();

        while max_val > PREC {
            rates += freqs / exp_index as f64;
            freqs *= -base_freqs;
            max_val *= exp;
            exp_index += 1;
        }

        TransitionTable::from_rates(&rates)
    }

    /// Calculates a table rate table from a matrix containing table counts.
    /// `counts[(a, b)]` is the number of observed transitions from `a` to `b`.
    pub fn from_counts(counts: &Mat21) -> TransitionTable {
        let mut freqs = zero_mat();
        for i in 0..MAT_SIZE {
            let tot = counts.row(i).sum();
            for j in 0..MAT_SIZE {
                freqs[(i, j)] = counts[(i, j)] / tot;
            }
        }

        TransitionTable::from_freqs(&freqs)
    }

    /// Gets the table of table probabilities for an edge of the given length.
    pub fn get_table(&self, mut length: f64) -> Mat21 {
        let mut res = zero_mat();
        let mut squaring_count = self.exp_factor;
        while length > MAX_EXPONENT {
            squaring_count += 1;
            length /= 2.;
        }
        let mut factor = 1.;
        for i in 0..EXP_ARRAY_LENGTH {
            res += self.exp_array[i] * factor;
            factor *= length;
        }

        for _ in 0..squaring_count {
            res *= res;
        }

        res
        // (length * self.rates).exp()
    }
}

impl Mul<f64> for TransitionTable {
    type Output = TransitionTable;

    fn mul(self, rhs: f64) -> Self::Output {
        TransitionTable::from_rates(&(self.rates * rhs))
    }
}

#[cfg(test)]
mod tests {
    use crate::transition::defs::zero_mat;
    use crate::transition::table::{TransitionTable, MAT_SIZE};

    /// Testing that taking log works.
    #[test]
    fn test_matrix_log() {
        let mut test = zero_mat();
        for i in 0..MAT_SIZE {
            test[(i, i)] = (i as f64 + 0.5) / (MAT_SIZE as f64);
        }
        let trans = TransitionTable::from_freqs(&test);

        for i in 0..MAT_SIZE {
            assert!((test[(i, i)].ln() - trans.rates[(i, i)]).abs() < 1e-5);
        }
    }

    /// Testing that calculating the exponent of a matrix works.
    #[test]
    fn test_matrix_exp() {
        let mut test = zero_mat();
        for i in 0..MAT_SIZE {
            test[(i, i)] = (i as f64) / (MAT_SIZE as f64) - 0.5;
        }
        let trans = TransitionTable::from_rates(&test);
        let trans_table = trans.get_table(1.);
        for i in 0..MAT_SIZE {
            assert!(
                (test[(i, i)].exp() - trans_table[(i, i)]).abs() < 1e-5,
                "Error too large: exp({}) = {} vs {}",
                test[(i, i)],
                test[(i, i)].exp(),
                trans_table[(i, i)]
            );
        }
    }
}

const _PAM1_DATA: [[f64; 21]; 21] = [
    [
        9867., 2., 9., 10., 3., 8., 17., 21., 2., 6., 4., 2., 6., 2., 22., 35., 32., 0., 2., 18.,
        0.1,
    ],
    [
        1., 9913., 1., 0., 1., 10., 0., 0., 10., 3., 1., 19., 4., 1., 4., 6., 1., 8., 0., 1., 0.1,
    ],
    [
        4., 1., 9822., 36., 0., 4., 6., 6., 21., 3., 1., 13., 0., 1., 2., 20., 9., 1., 4., 1., 0.1,
    ],
    [
        6., 0., 42., 9859., 0., 6., 53., 6., 4., 1., 0., 3., 0., 0., 1., 5., 3., 0., 0., 1., 0.1,
    ],
    [
        1., 1., 0., 0., 9973., 0., 0., 0., 1., 1., 0., 0., 0., 0., 1., 5., 1., 0., 3., 2., 0.1,
    ],
    [
        3., 9., 4., 5., 0., 9876., 27., 1., 23., 1., 3., 6., 4., 0., 6., 2., 2., 0., 0., 1., 0.1,
    ],
    [
        10., 0., 7., 56., 0., 35., 9865., 4., 2., 3., 1., 4., 1., 0., 3., 4., 2., 0., 1., 2., 0.1,
    ],
    [
        21., 1., 12., 11., 1., 3., 7., 9935., 1., 0., 1., 2., 1., 1., 3., 21., 3., 0., 0., 5., 0.1,
    ],
    [
        1., 8., 18., 3., 1., 20., 1., 0., 9912., 0., 1., 1., 0., 2., 3., 1., 1., 1., 4., 1., 0.1,
    ],
    [
        2., 2., 3., 1., 2., 1., 2., 0., 0., 9872., 9., 2., 12., 7., 0., 1., 7., 0., 1., 33., 0.1,
    ],
    [
        3., 1., 3., 0., 0., 6., 1., 1., 4., 22., 9947., 2., 45., 13., 3., 1., 3., 4., 2., 15., 0.1,
    ],
    [
        2., 37., 25., 6., 0., 12., 7., 2., 2., 4., 1., 9926., 20., 0., 3., 8., 11., 0., 1., 1., 0.1,
    ],
    [
        1., 1., 0., 0., 0., 2., 0., 0., 0., 5., 8., 4., 9874., 1., 0., 1., 2., 0., 0., 4., 0.1,
    ],
    [
        1., 1., 1., 0., 0., 0., 0., 1., 2., 8., 6., 0., 4., 9946., 0., 2., 1., 3., 28., 0., 0.1,
    ],
    [
        13., 5., 2., 1., 1., 8., 3., 2., 5., 1., 2., 2., 1., 1., 9926., 12., 4., 0., 0., 2., 0.1,
    ],
    [
        28., 11., 34., 7., 11., 4., 6., 16., 2., 2., 1., 7., 4., 3., 17., 9840., 38., 5., 2., 2.,
        0.1,
    ],
    [
        22., 2., 13., 4., 1., 3., 2., 2., 1., 11., 2., 8., 6., 1., 5., 32., 9871., 0., 2., 9., 0.1,
    ],
    [
        0., 2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 1., 0., 9976., 1., 0., 0.1,
    ],
    [
        1., 0., 3., 0., 3., 0., 1., 0., 4., 1., 1., 0., 0., 21., 0., 1., 1., 2., 9945., 1., 0.1,
    ],
    [
        13., 2., 1., 1., 3., 2., 2., 3., 3., 57., 11., 1., 17., 1., 3., 2., 10., 0., 2., 9901., 0.1,
    ],
    [
        0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 10000.,
    ],
];

/// Returns a table table built according to the pam1 substitution data.
pub fn pam1() -> TransitionTable {
    TransitionTable::from_counts(&Mat21::from(_PAM1_DATA))
}

/// Returns a table table built according to the pam1 substitution data.
pub fn pam100() -> TransitionTable {
    TransitionTable::from_counts(&Mat21::from(_PAM1_DATA)) * 100.
}
