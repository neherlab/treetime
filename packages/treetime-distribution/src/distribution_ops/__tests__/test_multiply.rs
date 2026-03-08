#[cfg(test)]
mod tests {
  use crate::DistributionPlain as Distribution;
  use crate::distribution_core::formula::DistributionFormula;
  use crate::distribution_ops::multiply::distribution_multiplication;
  use approx::assert_ulps_eq;
  use ndarray::array;

  /// Formula * Function returns a Function with correct pointwise products.
  ///
  /// Regression test for the fix in commit b3d97b34: `multiply_formula_function`
  /// must evaluate the Formula on the Function's grid and return a concrete
  /// Function, not a lazy Formula that panics on grid-based operations.
  #[test]
  fn test_multiply_formula_function_returns_function() {
    // Formula: f(t) = 2*t over [0, 10]
    let formula = DistributionFormula::new(|t| Ok(2.0 * t), 0.0, 10.0);
    let formula_dist = Distribution::Formula(formula);

    // Function on a 5-point grid [1, 3, 5, 7, 9] with values [1, 2, 3, 4, 5]
    let t = array![1.0, 3.0, 5.0, 7.0, 9.0];
    let y = array![1.0, 2.0, 3.0, 4.0, 5.0];
    let function_dist = Distribution::function(t, y).unwrap();

    let result = distribution_multiplication(&formula_dist, &function_dist).unwrap();

    // Result must be a Function, not a Formula
    let Distribution::Function(result_fn) = result else {
      panic!("Expected Function variant, got {result:?}");
    };

    // The grid is resampled to 5 points over [1, 9].
    // At each grid point t_i, the product is formula(t_i) * function(t_i) = 2*t_i * interp(t_i).
    // Grid: t_i = 1 + (9-1)*i/4 = 1, 3, 5, 7, 9
    // formula(1)=2, formula(3)=6, formula(5)=10, formula(7)=14, formula(9)=18
    // function values: 1, 2, 3, 4, 5
    // products: 2, 12, 30, 56, 90
    let expected = array![2.0, 12.0, 30.0, 56.0, 90.0];
    assert_eq!(result_fn.y().len(), expected.len());
    for (&actual, &exp) in result_fn.y().iter().zip(expected.iter()) {
      assert_ulps_eq!(exp, actual, max_ulps = 10);
    }
  }

  /// Commutativity: Function * Formula gives the same result as Formula * Function.
  #[test]
  fn test_multiply_function_formula_commutative() {
    let formula = DistributionFormula::new(|t| Ok(t * t), 0.0, 5.0);
    let formula_dist = Distribution::Formula(formula);

    let t = array![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    let y = array![1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
    let function_dist = Distribution::function(t, y).unwrap();

    let result_ff = distribution_multiplication(&formula_dist, &function_dist).unwrap();
    let result_fxf = distribution_multiplication(&function_dist, &formula_dist).unwrap();

    let (Distribution::Function(ff_fn), Distribution::Function(fxf_fn)) = (&result_ff, &result_fxf) else {
      panic!("Both results should be Function variants");
    };

    assert_eq!(ff_fn.y().len(), fxf_fn.y().len());
    for (&a, &b) in ff_fn.y().iter().zip(fxf_fn.y().iter()) {
      assert_ulps_eq!(a, b, max_ulps = 10);
    }
  }
}
