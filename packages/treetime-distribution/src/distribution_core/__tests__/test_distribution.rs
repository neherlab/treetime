#[cfg(test)]
mod tests {
  use crate::policy::{NegLog, Plain};
  use crate::{Distribution, DistributionFormula};
  use approx::assert_relative_eq;
  use ndarray::{Array1, array};
  use treetime_utils::pretty_assert_ulps_eq;

  #[test]
  fn test_distribution_max_value() {
    let point = Distribution::<Plain>::point(1.0, 5.0);
    assert_relative_eq!(point.max_value(), 5.0);

    let range = Distribution::<Plain>::range((0.0, 2.0), 3.0);
    assert_relative_eq!(range.max_value(), 3.0);

    let func = Distribution::<Plain>::function(array![0.0, 1.0, 2.0], array![1.0, 4.0, 2.0]).unwrap();
    assert_relative_eq!(func.max_value(), 4.0);

    let empty = Distribution::<Plain>::Empty;
    assert_relative_eq!(empty.max_value(), 0.0);
  }

  #[test]
  fn test_distribution_scale_by() {
    let point = Distribution::<Plain>::point(1.0, 5.0);
    let scaled = point.scale_by(2.0);
    if let Distribution::Point(p) = scaled {
      assert_relative_eq!(p.amplitude(), 10.0);
    } else {
      panic!("Expected Point");
    }

    let range = Distribution::<Plain>::range((0.0, 2.0), 3.0);
    let scaled = range.scale_by(0.5);
    if let Distribution::Range(r) = scaled {
      assert_relative_eq!(r.amplitude(), 1.5);
    } else {
      panic!("Expected Range");
    }
  }

  #[test]
  fn test_distribution_normalize() {
    let func = Distribution::<Plain>::function(array![0.0, 1.0, 2.0], array![2.0, 8.0, 4.0]).unwrap();
    let normalized = func.normalize();

    if let Distribution::Function(f) = normalized {
      assert_relative_eq!(f.y()[0], 0.25);
      assert_relative_eq!(f.y()[1], 1.0);
      assert_relative_eq!(f.y()[2], 0.5);
    } else {
      panic!("Expected Function");
    }
  }

  #[test]
  fn test_distribution_normalize_empty_on_zero() {
    let func = Distribution::<Plain>::function(array![0.0, 1.0, 2.0], array![0.0, 0.0, 0.0]).unwrap();
    let normalized = func.normalize();
    assert!(matches!(normalized, Distribution::Empty));
  }

  #[test]
  fn test_distribution_to_plain_normalized_empty() {
    let actual = Distribution::<NegLog>::Empty.to_plain_normalized();
    assert_eq!(Distribution::<Plain>::Empty, actual);
  }

  #[test]
  fn test_distribution_to_plain_normalized_point() {
    let actual = Distribution::<NegLog>::point(2.0, 1000.0).to_plain_normalized();
    let expected = Distribution::<Plain>::point(2.0, 1.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_distribution_to_plain_normalized_range() {
    let actual = Distribution::<NegLog>::range((1.0, 3.0), 1000.0).to_plain_normalized();
    let expected = Distribution::<Plain>::range((1.0, 3.0), 1.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_distribution_to_plain_normalized_function_preserves_likelihood_ratios() {
    let distribution = Distribution::<NegLog>::function(array![0.0, 1.0, 2.0], array![1004.0, 1000.0, 1003.0]).unwrap();
    let actual = distribution.to_plain_normalized();
    let expected = array![(-4.0_f64).exp(), 1.0, (-3.0_f64).exp()];
    pretty_assert_ulps_eq!(expected, actual.y(), max_ulps = 4);
  }

  #[test]
  fn test_distribution_to_plain_normalized_formula_normalizes_constant() {
    let formula = DistributionFormula::new(|_| Ok(1000.0), 0.0, 2.0);
    let actual = Distribution::<NegLog>::Formula(formula).to_plain_normalized();
    let expected = Array1::ones(200);
    pretty_assert_ulps_eq!(expected, actual.y(), max_ulps = 4);
  }

  #[test]
  fn test_distribution_to_plain_normalized_rejects_nonfinite_minimum() {
    let distribution = Distribution::<NegLog>::function(
      array![0.0, 1.0, 2.0],
      array![f64::INFINITY, f64::INFINITY, f64::INFINITY],
    )
    .unwrap();
    let actual = distribution.to_plain_normalized();
    assert_eq!(Distribution::<Plain>::Empty, actual);
  }
}
