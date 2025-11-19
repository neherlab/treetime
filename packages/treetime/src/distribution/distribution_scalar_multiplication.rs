use crate::distribution::distribution::Distribution;
use eyre::Report;

/// Multiplies a distribution by a scalar value.
///
/// For all distribution types, multiplies the amplitude values by the scalar.
/// This scales the probability density or likelihood represented by the distribution.
pub fn distribution_scalar_multiplication(dist: &Distribution, scalar: f64) -> Result<Distribution, Report> {
  match dist {
    Distribution::Function(f) => Distribution::function(f.t(), f.y() * scalar),
    Distribution::Point(p) => Ok(Distribution::point(p.t(), p.amplitude() * scalar)),
    Distribution::Range(r) => Ok(Distribution::range((r.start(), r.end()), r.amplitude() * scalar)),
    Distribution::Empty => Ok(Distribution::empty()),
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_abs_diff_eq;
  use ndarray::array;

  #[test]
  fn test_distribution_scalar_multiplication_function_positive_scalar() {
    let t = array![0.0, 1.0, 2.0];
    let y = array![1.0, 2.0, 3.0];
    let dist = Distribution::function(t.clone(), y.clone()).unwrap();

    let result = distribution_scalar_multiplication(&dist, 2.5).unwrap();

    if let Distribution::Function(f) = result {
      assert_eq!(f.t(), &t);
      assert_abs_diff_eq!(f.y(), &(y * 2.5), epsilon = 1e-10);
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
      assert_abs_diff_eq!(f.y(), &array![0.0, 0.0, 0.0], epsilon = 1e-10);
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
      assert_abs_diff_eq!(f.y(), &(y * -1.5), epsilon = 1e-10);
    } else {
      panic!("Expected Function distribution");
    }
  }

  #[test]
  fn test_distribution_scalar_multiplication_point_positive_scalar() {
    let dist = Distribution::point(5.0, 3.0);

    let result = distribution_scalar_multiplication(&dist, 2.0).unwrap();

    if let Distribution::Point(p) = result {
      assert_abs_diff_eq!(p.t(), 5.0, epsilon = 1e-10);
      assert_abs_diff_eq!(p.amplitude(), 6.0, epsilon = 1e-10);
    } else {
      panic!("Expected Point distribution");
    }
  }

  #[test]
  fn test_distribution_scalar_multiplication_point_zero_scalar() {
    let dist = Distribution::point(5.0, 3.0);

    let result = distribution_scalar_multiplication(&dist, 0.0).unwrap();

    if let Distribution::Point(p) = result {
      assert_abs_diff_eq!(p.t(), 5.0, epsilon = 1e-10);
      assert_abs_diff_eq!(p.amplitude(), 0.0, epsilon = 1e-10);
    } else {
      panic!("Expected Point distribution");
    }
  }

  #[test]
  fn test_distribution_scalar_multiplication_point_negative_scalar() {
    let dist = Distribution::point(5.0, 3.0);

    let result = distribution_scalar_multiplication(&dist, -2.0).unwrap();

    if let Distribution::Point(p) = result {
      assert_abs_diff_eq!(p.t(), 5.0, epsilon = 1e-10);
      assert_abs_diff_eq!(p.amplitude(), -6.0, epsilon = 1e-10);
    } else {
      panic!("Expected Point distribution");
    }
  }

  #[test]
  fn test_distribution_scalar_multiplication_range_positive_scalar() {
    let dist = Distribution::range((1.0, 3.0), 2.0);

    let result = distribution_scalar_multiplication(&dist, 1.5).unwrap();

    if let Distribution::Range(r) = result {
      assert_abs_diff_eq!(r.start(), 1.0, epsilon = 1e-10);
      assert_abs_diff_eq!(r.end(), 3.0, epsilon = 1e-10);
      assert_abs_diff_eq!(r.amplitude(), 3.0, epsilon = 1e-10);
    } else {
      panic!("Expected Range distribution");
    }
  }

  #[test]
  fn test_distribution_scalar_multiplication_range_zero_scalar() {
    let dist = Distribution::range((1.0, 3.0), 2.0);

    let result = distribution_scalar_multiplication(&dist, 0.0).unwrap();

    if let Distribution::Range(r) = result {
      assert_abs_diff_eq!(r.start(), 1.0, epsilon = 1e-10);
      assert_abs_diff_eq!(r.end(), 3.0, epsilon = 1e-10);
      assert_abs_diff_eq!(r.amplitude(), 0.0, epsilon = 1e-10);
    } else {
      panic!("Expected Range distribution");
    }
  }

  #[test]
  fn test_distribution_scalar_multiplication_range_negative_scalar() {
    let dist = Distribution::range((1.0, 3.0), 2.0);

    let result = distribution_scalar_multiplication(&dist, -0.5).unwrap();

    if let Distribution::Range(r) = result {
      assert_abs_diff_eq!(r.start(), 1.0, epsilon = 1e-10);
      assert_abs_diff_eq!(r.end(), 3.0, epsilon = 1e-10);
      assert_abs_diff_eq!(r.amplitude(), -1.0, epsilon = 1e-10);
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
      assert_abs_diff_eq!(f.y(), &array![1.0, 2.0, 3.0], epsilon = 1e-10);
    } else {
      panic!("Expected Function distribution");
    }
  }
}
