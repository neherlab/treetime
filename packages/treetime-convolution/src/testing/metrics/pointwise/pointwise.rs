use crate::testing::metrics::config::PointwiseConfig;
use crate::testing::metrics::pointwise::errors::{PointwiseErrors, compute_pointwise_errors};
use crate::testing::metrics::pointwise::structural::{StructuralErrors, compute_structural_errors};
use crate::testing::metrics::pointwise::tolerance::{ToleranceMetrics, compute_tolerance_metrics};
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use treetime_utils::make_error;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PointwiseMetrics {
  pub total_points: usize,
  pub dx: f64,
  pub errors: PointwiseErrors,
  pub structural: StructuralErrors,
  pub tolerance: ToleranceMetrics,
}

impl PointwiseMetrics {
  pub fn new(x: &Array1<f64>, actual: &Array1<f64>, expected: &Array1<f64>) -> eyre::Result<Self> {
    Self::new_with_config(x, actual, expected, &PointwiseConfig::default())
  }

  pub fn new_with_config(
    x: &Array1<f64>,
    actual: &Array1<f64>,
    expected: &Array1<f64>,
    config: &PointwiseConfig,
  ) -> eyre::Result<Self> {
    validate_inputs(x, actual, expected)?;

    let n = x.len();
    let dx = compute_grid_spacing(x)?;

    let errors = compute_pointwise_errors(actual, expected, config)?;
    let structural = compute_structural_errors(x, actual, expected, dx, config)?;
    let tolerance = compute_tolerance_metrics(actual, expected, config)?;

    Ok(Self {
      total_points: n,
      dx,
      errors,
      structural,
      tolerance,
    })
  }
}

fn validate_inputs(x: &Array1<f64>, actual: &Array1<f64>, expected: &Array1<f64>) -> eyre::Result<()> {
  let n = x.len();
  if n < 2 {
    return make_error!("Input arrays must have at least 2 points, got {n}");
  }
  if actual.len() != n {
    return make_error!("Actual array length {} does not match grid length {n}", actual.len());
  }
  if expected.len() != n {
    return make_error!(
      "Expected array length {} does not match grid length {n}",
      expected.len()
    );
  }
  Ok(())
}

fn compute_grid_spacing(x: &Array1<f64>) -> eyre::Result<f64> {
  if x.len() < 2 {
    return make_error!("Cannot compute grid spacing with fewer than 2 points");
  }
  let dx = x[1] - x[0];
  if !dx.is_finite() || dx <= 0.0 {
    return make_error!("Invalid grid spacing: {dx}");
  }
  Ok(dx)
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_abs_diff_eq;
  use ndarray::array;

  #[test]
  fn test_perfect_agreement() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y = array![1.0, 2.0, 3.0, 2.0, 1.0];
    let result = PointwiseMetrics::new(&x, &y, &y).unwrap();

    assert_eq!(result.total_points, 5);
    assert_abs_diff_eq!(result.errors.summary.abs_max, 0.0, epsilon = 1e-12);
    assert_abs_diff_eq!(result.errors.summary.rel_max, 0.0, epsilon = 1e-12);
    assert_abs_diff_eq!(result.errors.summary.signed_bias, 0.0, epsilon = 1e-12);
  }

  #[test]
  fn test_constant_offset() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let expected = array![1.0, 2.0, 3.0, 2.0, 1.0];
    let actual = &expected + 0.5;
    let result = PointwiseMetrics::new(&x, &actual, &expected).unwrap();

    for i in 0..5 {
      assert_abs_diff_eq!(result.errors.absolute[i], 0.5, epsilon = 1e-12);
      assert_abs_diff_eq!(result.errors.signed[i], 0.5, epsilon = 1e-12);
    }
    assert_abs_diff_eq!(result.errors.summary.signed_bias, 0.5, epsilon = 1e-12);
  }
}
