use argmin::prelude::*;
use argmin_testfunctions::{rosenbrock_2d, rosenbrock_2d_derivative, rosenbrock_2d_hessian};

/// First, create a struct for your problem
pub struct Rosenbrock {
    pub a: f64,
    pub b: f64,
}

/// Implement `ArgminOp` for `Rosenbrock`
impl ArgminOp for Rosenbrock {
    /// Type of the parameter vector
    type Param = Vec<f64>;
    /// Type of the return value computed by the cost function
    type Output = f64;
    /// Type of the Hessian. Can be `()` if not needed.
    type Hessian = Vec<Vec<f64>>;
    /// Type of the Jacobian. Can be `()` if not needed.
    type Jacobian = ();
    /// Floating point precision
    type Float = f64;

    /// Apply the cost function to a parameter `p`
    fn apply(&self, p: &Self::Param) -> Result<Self::Output, Error> {
        Ok(rosenbrock_2d(p, self.a, self.b))
    }

    /// Compute the gradient at parameter `p`.
    fn gradient(&self, p: &Self::Param) -> Result<Self::Param, Error> {
        Ok(rosenbrock_2d_derivative(p, self.a, self.b))
    }

    /// Compute the Hessian at parameter `p`.
    fn hessian(&self, p: &Self::Param) -> Result<Self::Hessian, Error> {
        let t = rosenbrock_2d_hessian(p, self.a, self.b);
        Ok(vec![vec![t[0], t[1]], vec![t[2], t[3]]])
    }
}
