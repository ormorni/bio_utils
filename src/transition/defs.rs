extern crate nalgebra as na;

/// A module containing useful definitions of matrices and vectors.

// Defining the properties of the table matrices for easy customization.
pub const MAT_SIZE: usize = 21;
pub type MatSize = na::U21;
pub type Float = f64;

// Defining the table matrix types.
pub type Mat21Storage = na::ArrayStorage<Float, MAT_SIZE, MAT_SIZE>;
pub type Mat21 = na::Matrix<Float, MatSize, MatSize, Mat21Storage>;

// Defining the vector types.
pub type Vec21Storage = na::ArrayStorage<Float, MAT_SIZE, 1>;
pub type Vec21 = na::Matrix<Float, MatSize, na::U1, Vec21Storage>;

/// Returns the zero matrix.
pub fn zero_mat() -> Mat21 {
    Mat21::default()
}

/// Returns the zero vector.
pub fn zero_vec() -> Vec21 {
    Vec21::default()
}

/// Returns the zero vector.
pub fn one_vec() -> Vec21 {
    let mut res = Vec21::default();
    for i in 0..21 {
        res[i] = 1.;
    }
    res
}

/// Returns the zero vector.
pub fn uniform() -> Vec21 {
    let mut res = Vec21::default();
    for i in 0..MAT_SIZE {
        res[i] = (1 as Float) / (MAT_SIZE as Float);
    }
    res
}

/// Returns the identity matrix.
pub fn identity() -> Mat21 {
    let mut res = Mat21::default();
    res.fill_with_identity();
    res
}

#[cfg(test)]
mod tests {
    use crate::transition::defs::{zero_mat, MAT_SIZE};

    /// Testing the zero function, as it relies on some undefined default behavior.
    #[test]
    fn test_zero_mat() {
        let mat = zero_mat();
        for i in 0..MAT_SIZE {
            for j in 0..MAT_SIZE {
                assert_eq!(mat[(i, j)], 0.);
            }
        }
    }
}
