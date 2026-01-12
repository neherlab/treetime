use crate::distribution::distribution_multiplication::distribution_multiplication;
use crate::distribution::scaled_distribution::ScaledDistribution;
use eyre::Report;

/// Multiply two scaled distributions.
///
/// Delegates to distribution_multiplication for the inner (normalized) distributions,
/// then combines log_scales.
pub fn scaled_distribution_multiplication(
  a: &ScaledDistribution,
  b: &ScaledDistribution,
) -> Result<ScaledDistribution, Report> {
  if a.is_empty() || b.is_empty() {
    return Ok(ScaledDistribution::default());
  }

  let product_inner = distribution_multiplication(a.inner(), b.inner())?;

  let max_val = product_inner.max_value();
  if max_val <= 0.0 || !max_val.is_finite() {
    return Ok(ScaledDistribution::default());
  }

  let combined_log_scale = a.log_scale() + b.log_scale() + max_val.ln();
  let normalized_product = product_inner.normalize();

  Ok(ScaledDistribution::from_parts(combined_log_scale, normalized_product))
}

/// Multiply many scaled distributions.
///
/// More numerically stable than repeated pairwise multiplication
/// because scales accumulate in log space.
pub fn scaled_distribution_multiply_many(distributions: &[&ScaledDistribution]) -> Result<ScaledDistribution, Report> {
  match distributions {
    [] => Ok(ScaledDistribution::default()),
    [single] => Ok((*single).clone()),
    [first, rest @ ..] => {
      let mut result = (*first).clone();
      for dist in rest {
        result = scaled_distribution_multiplication(&result, dist)?;
        if result.is_empty() {
          return Ok(ScaledDistribution::default());
        }
      }
      Ok(result)
    },
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::distribution::distribution::Distribution;
  use crate::distribution::y_axis_policy::Plain;
  use approx::assert_relative_eq;
  use ndarray::{Array1, array};

  fn make_point(t: f64, amplitude: f64) -> ScaledDistribution {
    let dist = Distribution::<Plain>::point(t, amplitude);
    ScaledDistribution::from_plain(&dist)
  }

  fn make_function(x: Array1<f64>, y: Array1<f64>) -> ScaledDistribution {
    let dist = Distribution::<Plain>::function(x, y).unwrap();
    ScaledDistribution::from_plain(&dist)
  }

  #[test]
  fn test_scaled_distribution_multiply_points_same_location() {
    let a = make_point(2.0, 3.0);
    let b = make_point(2.0, 4.0);
    let result = scaled_distribution_multiplication(&a, &b).unwrap();

    assert_relative_eq!(result.peak_value(), 12.0, epsilon = 1e-10);
  }

  #[test]
  fn test_scaled_distribution_multiply_points_different_location() {
    let a = make_point(2.0, 3.0);
    let b = make_point(5.0, 4.0);
    let result = scaled_distribution_multiplication(&a, &b).unwrap();

    assert!(result.is_empty());
  }

  #[test]
  fn test_scaled_distribution_multiply_functions() {
    let a = make_function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 1.0]);
    let b = make_function(array![0.0, 1.0, 2.0], array![2.0, 1.0, 2.0]);
    let result = scaled_distribution_multiplication(&a, &b).unwrap();

    assert!(!result.is_empty());
    assert!(result.log_scale().is_finite());
  }

  #[test]
  fn test_scaled_distribution_multiply_many_no_underflow() {
    let dists: Vec<ScaledDistribution> = std::iter::repeat_with(|| make_point(0.0, 0.0001)).take(100).collect();
    let refs: Vec<&ScaledDistribution> = dists.iter().collect();

    let result = scaled_distribution_multiply_many(&refs).unwrap();

    assert!(!result.is_empty());
    assert!(result.log_scale().is_finite());

    let expected_log_scale = 100.0 * 0.0001_f64.ln();
    assert_relative_eq!(result.log_scale(), expected_log_scale, epsilon = 1e-10);
  }

  #[test]
  fn test_scaled_distribution_multiply_empty_propagates() {
    let a = ScaledDistribution::default();
    let b = make_point(1.0, 2.0);
    let result = scaled_distribution_multiplication(&a, &b).unwrap();
    assert!(result.is_empty());
  }

  #[test]
  fn test_scaled_distribution_multiply_preserves_normalization() {
    let a = make_function(array![0.0, 1.0, 2.0], array![10.0, 40.0, 20.0]);
    let b = make_function(array![0.0, 1.0, 2.0], array![5.0, 10.0, 5.0]);

    let result = scaled_distribution_multiplication(&a, &b).unwrap();

    assert_relative_eq!(result.inner().max_value(), 1.0, epsilon = 1e-10);
  }

  #[test]
  fn test_scaled_distribution_multiply_many_empty() {
    let result = scaled_distribution_multiply_many(&[]).unwrap();
    assert!(result.is_empty());
  }

  #[test]
  fn test_scaled_distribution_multiply_many_single() {
    let a = make_point(1.0, 5.0);
    let result = scaled_distribution_multiply_many(&[&a]).unwrap();

    assert_relative_eq!(result.log_scale(), a.log_scale(), epsilon = 1e-10);
  }
}
