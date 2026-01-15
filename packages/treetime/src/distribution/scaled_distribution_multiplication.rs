use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_function::DistributionFunction;
use crate::distribution::distribution_multiplication::distribution_multiplication;
use crate::distribution::scaled_distribution::ScaledDistribution;
use crate::distribution::y_axis_policy::Plain;
use approx::ulps_eq;
use eyre::Report;
use itertools::Itertools;
use ndarray::Array1;
use treetime_ops::multiply_many_lazy_normalize;

struct AlignedFunctionArrays<'a> {
  arrays: Vec<&'a Array1<f64>>,
  x_min: f64,
  dx: f64,
}

fn try_extract_aligned_function_arrays<'a>(
  distributions: &'a [&'a ScaledDistribution],
) -> Option<AlignedFunctionArrays<'a>> {
  let first = distributions.first()?;
  let Distribution::Function(first_func) = first.inner() else {
    return None;
  };

  let x_min = first_func.x_min();
  let dx = first_func.dx();
  let len = first_func.len();

  let mut arrays = Vec::with_capacity(distributions.len());

  for dist in distributions {
    match dist.inner() {
      Distribution::Function(f) => {
        if !ulps_eq!(f.x_min(), x_min, max_ulps = 10)
          || !ulps_eq!(f.dx(), dx, max_ulps = 10)
          || f.len() != len
        {
          return None;
        }
        arrays.push(f.y());
      },
      _ => return None,
    }
  }

  Some(AlignedFunctionArrays { arrays, x_min, dx })
}

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
/// When all inputs are Function distributions with matching grids, delegates
/// to `treetime_ops::multiply_many_lazy_normalize` for validated numerical
/// stability. Otherwise falls back to pairwise multiplication.
pub fn scaled_distribution_multiply_many(distributions: &[&ScaledDistribution]) -> Result<ScaledDistribution, Report> {
  match distributions {
    [] => Ok(ScaledDistribution::default()),
    [single] => Ok((*single).clone()),
    _ => {
      if distributions.iter().any(|d| d.is_empty()) {
        return Ok(ScaledDistribution::default());
      }

      let total_input_log_scale: f64 = distributions.iter().map(|d| d.log_scale()).sum();

      if let Some(aligned) = try_extract_aligned_function_arrays(distributions) {
        let array_refs = aligned.arrays.iter().copied().collect_vec();
        let result = multiply_many_lazy_normalize(&array_refs);

        if !result.log_scale.is_finite() {
          return Ok(ScaledDistribution::default());
        }

        let result_func =
          DistributionFunction::<f64, Plain>::from_start_dx_values(aligned.x_min, aligned.dx, result.normalized)?;
        let total_log_scale = total_input_log_scale + result.log_scale;

        return Ok(ScaledDistribution::from_parts(total_log_scale, Distribution::Function(result_func)));
      }

      let mut product = distributions[0].inner().clone();
      for dist in distributions.iter().skip(1) {
        product = distribution_multiplication(&product, dist.inner())?;
      }

      let max_val = product.max_value();
      if max_val <= 0.0 || !max_val.is_finite() {
        return Ok(ScaledDistribution::default());
      }

      Ok(ScaledDistribution::from_parts(
        total_input_log_scale + max_val.ln(),
        product.normalize(),
      ))
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

  #[test]
  fn test_scaled_distribution_multiply_many_functions_same_grid() {
    let a = make_function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 1.0]);
    let b = make_function(array![0.0, 1.0, 2.0], array![0.5, 1.0, 0.5]);
    let c = make_function(array![0.0, 1.0, 2.0], array![0.5, 1.0, 0.5]);

    let result = scaled_distribution_multiply_many(&[&a, &b, &c]).unwrap();

    assert!(!result.is_empty());
    assert!(result.log_scale().is_finite());
    assert_relative_eq!(result.inner().max_value(), 1.0, epsilon = 1e-10);

    if let Distribution::Function(f) = result.inner() {
      assert_relative_eq!(f.y()[1], 1.0, epsilon = 1e-10);
      assert_relative_eq!(f.y()[0], 0.125, epsilon = 1e-10);
      assert_relative_eq!(f.y()[2], 0.125, epsilon = 1e-10);
    } else {
      panic!("Expected Function distribution");
    }
  }

  #[test]
  fn test_scaled_distribution_multiply_many_functions_underflow_resistance() {
    let small_peak = 1e-50;
    let dists: Vec<ScaledDistribution> = std::iter::repeat_with(|| {
      make_function(array![0.0, 1.0, 2.0], array![small_peak * 0.5, small_peak, small_peak * 0.5])
    })
    .take(10)
    .collect();
    let refs: Vec<&ScaledDistribution> = dists.iter().collect();

    let result = scaled_distribution_multiply_many(&refs).unwrap();

    assert!(!result.is_empty());
    assert!(result.log_scale().is_finite());

    let expected_log_scale = 10.0 * small_peak.ln();
    assert_relative_eq!(result.log_scale(), expected_log_scale, epsilon = 1e-10);
  }

  #[test]
  fn test_scaled_distribution_multiply_many_mixed_types_fallback() {
    let func = make_function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 1.0]);
    let point = make_point(1.0, 3.0);

    let result = scaled_distribution_multiply_many(&[&func, &point]).unwrap();

    assert!(!result.is_empty());
    assert!(result.log_scale().is_finite());
  }

  #[test]
  fn test_scaled_distribution_multiply_many_different_grids_fallback() {
    let a = make_function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 1.0]);
    let b = make_function(array![0.0, 0.5, 1.0, 1.5, 2.0], array![1.0, 1.5, 2.0, 1.5, 1.0]);

    let result = scaled_distribution_multiply_many(&[&a, &b]).unwrap();

    assert!(!result.is_empty());
    assert!(result.log_scale().is_finite());
  }
}
