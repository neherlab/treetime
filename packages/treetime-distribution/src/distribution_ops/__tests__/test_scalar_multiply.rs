#[cfg(test)]
mod tests {
  use crate::DistributionPlain as Distribution;
  use crate::distribution_core::formula::DistributionFormula;
  use crate::distribution_core::function::DistributionFunction;
  use crate::distribution_ops::scalar_multiply::distribution_scalar_multiplication;
  use approx::assert_ulps_eq;
  use ndarray::array;
  use treetime_grid::{BoundaryBehavior, SoftTailLaw};
  use treetime_utils::assert_error;

  #[test]
  fn test_scalar_multiply_formula_returns_error() {
    // Oracle: kb/issues/H-distribution-result-api-panics-on-formula.md.
    let formula = Distribution::Formula(DistributionFormula::new(|_| Ok(1.0), 0.0, 1.0));

    assert_error!(
      distribution_scalar_multiplication(&formula, 2.0),
      "Cannot multiply Formula by scalar: operation not implemented"
    );
  }

  #[test]
  fn test_distribution_scalar_multiplication_function_positive_scalar() {
    let t = array![0.0, 1.0, 2.0];
    let y = array![1.0, 2.0, 3.0];
    let dist = Distribution::function(t.clone(), y.clone()).unwrap();

    let result = distribution_scalar_multiplication(&dist, 2.5).unwrap();

    if let Distribution::Function(f) = result {
      assert_eq!(f.t(), &t);
      assert_ulps_eq!(f.y(), &(y * 2.5), max_ulps = 4);
    } else {
      panic!("Expected Function distribution");
    }
  }

  #[test]
  fn test_distribution_scalar_multiplication_function_zero_scalar() {
    let t = array![0.0, 1.0, 2.0];
    let y = array![1.0, 2.0, 3.0];
    let dist = Distribution::function(t.clone(), y).unwrap();

    let result = distribution_scalar_multiplication(&dist, 0.0).unwrap();

    if let Distribution::Function(f) = result {
      assert_eq!(f.t(), &t);
      assert_ulps_eq!(f.y(), &array![0.0, 0.0, 0.0], max_ulps = 4);
    } else {
      panic!("Expected Function distribution");
    }
  }

  #[test]
  fn test_distribution_scalar_multiplication_function_negative_scalar() {
    let t = array![0.0, 1.0, 2.0];
    let y = array![1.0, 2.0, 3.0];
    let dist = Distribution::function(t.clone(), y.clone()).unwrap();

    let result = distribution_scalar_multiplication(&dist, -1.5).unwrap();

    if let Distribution::Function(f) = result {
      assert_eq!(f.t(), &t);
      assert_ulps_eq!(f.y(), &(y * -1.5), max_ulps = 4);
    } else {
      panic!("Expected Function distribution");
    }
  }

  #[test]
  fn test_distribution_scalar_multiplication_point_positive_scalar() {
    let dist = Distribution::point(5.0, 3.0);

    let result = distribution_scalar_multiplication(&dist, 2.0).unwrap();

    if let Distribution::Point(p) = result {
      assert_ulps_eq!(p.t(), 5.0, max_ulps = 4);
      assert_ulps_eq!(p.amplitude(), 6.0, max_ulps = 4);
    } else {
      panic!("Expected Point distribution");
    }
  }

  #[test]
  fn test_distribution_scalar_multiplication_point_zero_scalar() {
    let dist = Distribution::point(5.0, 3.0);

    let result = distribution_scalar_multiplication(&dist, 0.0).unwrap();

    if let Distribution::Point(p) = result {
      assert_ulps_eq!(p.t(), 5.0, max_ulps = 4);
      assert_ulps_eq!(p.amplitude(), 0.0, max_ulps = 4);
    } else {
      panic!("Expected Point distribution");
    }
  }

  #[test]
  fn test_distribution_scalar_multiplication_point_negative_scalar() {
    let dist = Distribution::point(5.0, 3.0);

    let result = distribution_scalar_multiplication(&dist, -2.0).unwrap();

    if let Distribution::Point(p) = result {
      assert_ulps_eq!(p.t(), 5.0, max_ulps = 4);
      assert_ulps_eq!(p.amplitude(), -6.0, max_ulps = 4);
    } else {
      panic!("Expected Point distribution");
    }
  }

  #[test]
  fn test_distribution_scalar_multiplication_range_positive_scalar() {
    let dist = Distribution::range((1.0, 3.0), 2.0);

    let result = distribution_scalar_multiplication(&dist, 1.5).unwrap();

    if let Distribution::Range(r) = result {
      assert_ulps_eq!(r.start(), 1.0, max_ulps = 4);
      assert_ulps_eq!(r.end(), 3.0, max_ulps = 4);
      assert_ulps_eq!(r.amplitude(), 3.0, max_ulps = 4);
    } else {
      panic!("Expected Range distribution");
    }
  }

  #[test]
  fn test_distribution_scalar_multiplication_range_zero_scalar() {
    let dist = Distribution::range((1.0, 3.0), 2.0);

    let result = distribution_scalar_multiplication(&dist, 0.0).unwrap();

    if let Distribution::Range(r) = result {
      assert_ulps_eq!(r.start(), 1.0, max_ulps = 4);
      assert_ulps_eq!(r.end(), 3.0, max_ulps = 4);
      assert_ulps_eq!(r.amplitude(), 0.0, max_ulps = 4);
    } else {
      panic!("Expected Range distribution");
    }
  }

  #[test]
  fn test_distribution_scalar_multiplication_range_negative_scalar() {
    let dist = Distribution::range((1.0, 3.0), 2.0);

    let result = distribution_scalar_multiplication(&dist, -0.5).unwrap();

    if let Distribution::Range(r) = result {
      assert_ulps_eq!(r.start(), 1.0, max_ulps = 4);
      assert_ulps_eq!(r.end(), 3.0, max_ulps = 4);
      assert_ulps_eq!(r.amplitude(), -1.0, max_ulps = 4);
    } else {
      panic!("Expected Range distribution");
    }
  }

  #[test]
  fn test_distribution_scalar_multiplication_empty() {
    let dist = Distribution::empty();

    let result = distribution_scalar_multiplication(&dist, 5.0).unwrap();

    assert!(matches!(result, Distribution::Empty));
  }

  #[test]
  fn test_distribution_scalar_multiplication_function_fractional_scalar() {
    let t = array![0.0, 1.0, 2.0];
    let y = array![10.0, 20.0, 30.0];
    let dist = Distribution::function(t.clone(), y).unwrap();

    let result = distribution_scalar_multiplication(&dist, 0.1).unwrap();

    if let Distribution::Function(f) = result {
      assert_eq!(f.t(), &t);
      assert_ulps_eq!(f.y(), &array![1.0, 2.0, 3.0], max_ulps = 4);
    } else {
      panic!("Expected Function distribution");
    }
  }

  /// Scaling the probability by a constant leaves both fitted tail laws unchanged.
  ///
  /// A constant probability scale is a constant shift of the negative-log ordinate, which both the
  /// soft slope and the hard edge are invariant to. The input carries a soft `Linear` left law
  /// (slope `-0.7`) and a `Hard` right edge; after scaling by `2.5` both are byte-identical. Before
  /// the fix the rebuild dropped both sides to `Error`.
  #[test]
  fn test_distribution_scalar_multiplication_preserves_tail_laws() {
    let left = BoundaryBehavior::Linear(SoftTailLaw { slope: -0.7 });
    let input = DistributionFunction::from_start_dx_values(0.0, 1.0, array![1.0, 2.0, 3.0, 4.0])
      .unwrap()
      .with_left_extrap(left)
      .unwrap()
      .with_right_extrap(BoundaryBehavior::Hard)
      .unwrap();

    let result = distribution_scalar_multiplication(&Distribution::Function(input), 2.5).unwrap();

    let Distribution::Function(f) = result else {
      panic!("Expected Function distribution");
    };
    assert_eq!(left, f.left_extrap());
    assert_eq!(BoundaryBehavior::Hard, f.right_extrap());
    assert_ulps_eq!(f.y(), &array![2.5, 5.0, 7.5, 10.0], max_ulps = 4);
  }
}
