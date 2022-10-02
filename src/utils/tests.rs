pub const ABS_TOL: f64 = 1e-10;

/// Asserts that two values are close to each other.
pub fn assert_close(a: f64, b: f64) {
    if a.abs() > ABS_TOL && b.abs() > ABS_TOL {
        // If both values are larger than rounding errors, perform a check.
        assert!(
            (a - b).abs() / (a.abs() + b.abs()) < 0.1,
            "Values are not close: {} {}",
            a,
            b
        );
    } else {
        // Assert that both are zero.
        assert!(a.abs() <= ABS_TOL && b.abs() <= ABS_TOL, "Values are not close: {}, {}", a, b);
    }
}

pub const SEED: [u8; 32] = [
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
    27, 28, 29, 30, 31, 32,
];
