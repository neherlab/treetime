use crate::distribution::distribution_convolution::distribution_convolution;
use crate::distribution::scaled_distribution::ScaledDistribution;
use eyre::Report;

/// Convolve two scaled distributions.
///
/// Delegates to distribution_convolution for the inner (normalized) distributions,
/// then combines log_scales.
pub fn scaled_distribution_convolution(
  a: &ScaledDistribution,
  b: &ScaledDistribution,
) -> Result<ScaledDistribution, Report> {
  if a.is_empty() || b.is_empty() {
    return Ok(ScaledDistribution::default());
  }

  let conv_inner = distribution_convolution(a.inner(), b.inner())?;

  let max_val = conv_inner.max_value();
  if max_val <= 0.0 || !max_val.is_finite() {
    return Ok(ScaledDistribution::default());
  }

  let combined_log_scale = a.log_scale() + b.log_scale() + max_val.ln();
  let normalized_conv = conv_inner.normalize();

  Ok(ScaledDistribution::from_parts(combined_log_scale, normalized_conv))
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::distribution::distribution::Distribution;
  use crate::distribution::y_axis_policy::Plain;
  use approx::assert_relative_eq;
  use ndarray::{array, Array1};

  fn make_point(t: f64, amplitude: f64) -> ScaledDistribution {
    let dist = Distribution::<Plain>::point(t, amplitude);
    ScaledDistribution::from_plain(&dist)
  }

  fn make_function(x: Array1<f64>, y: Array1<f64>) -> ScaledDistribution {
    let dist = Distribution::<Plain>::function(x, y).unwrap();
    ScaledDistribution::from_plain(&dist)
  }

  #[test]
  fn test_scaled_distribution_convolve_points() {
    let a = make_point(2.0, 3.0);
    let b = make_point(5.0, 4.0);
    let result = scaled_distribution_convolution(&a, &b).unwrap();

    assert_relative_eq!(result.peak_value(), 12.0, epsilon = 1e-10);
  }

  #[test]
  fn test_scaled_distribution_convolve_functions() {
    let a = make_function(array![0.0, 1.0, 2.0], array![1.0, 1.0, 1.0]);
    let b = make_function(array![0.0, 1.0], array![1.0, 1.0]);
    let result = scaled_distribution_convolution(&a, &b).unwrap();

    assert!(!result.is_empty());
    assert!(result.log_scale().is_finite());
  }

  #[test]
  fn test_scaled_distribution_convolve_empty() {
    let a = ScaledDistribution::default();
    let b = make_point(1.0, 2.0);
    let result = scaled_distribution_convolution(&a, &b).unwrap();
    assert!(result.is_empty());
  }

  #[test]
  fn test_scaled_distribution_convolve_preserves_normalization() {
    let a = make_function(array![0.0, 1.0, 2.0], array![10.0, 40.0, 10.0]);
    let b = make_function(array![0.0, 1.0, 2.0], array![5.0, 20.0, 5.0]);
    let result = scaled_distribution_convolution(&a, &b).unwrap();

    assert_relative_eq!(result.inner().max_value(), 1.0, epsilon = 1e-10);
  }

  #[test]
  fn test_scaled_distribution_convolve_symmetric() {
    let a = make_function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 1.0]);
    let b = make_function(array![0.0, 1.0, 2.0], array![2.0, 1.0, 2.0]);

    let result_ab = scaled_distribution_convolution(&a, &b).unwrap();
    let result_ba = scaled_distribution_convolution(&b, &a).unwrap();

    assert_relative_eq!(result_ab.log_scale(), result_ba.log_scale(), epsilon = 1e-10);
  }
}
