use crate::make_error;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::BTreeMap;

/// Pure pointwise metrics - one value per evaluation grid point
///
/// These metrics maintain a strict 1:1 correspondence with the evaluation grid,
/// providing error analysis at each individual point without any aggregation
/// or spatial relationships.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PointwiseMetrics {
  /// Total number of evaluation points
  pub total_points: usize,
  /// Grid spacing (assumed uniform)
  pub dx: f64,
  /// Basic pointwise error metrics
  pub errors: PointwiseErrors,
  /// Structural derivative and symmetry analysis
  pub structural: StructuralErrors,
  /// Tolerance compliance at each point
  pub tolerance: ToleranceMetrics,
}

impl PointwiseMetrics {
  /// Creates new pointwise metrics from evaluation grid and function values
  pub fn new(x: &Array1<f64>, actual: &Array1<f64>, expected: &Array1<f64>) -> eyre::Result<Self> {
    Self::new_with_config(x, actual, expected, &PointwiseConfig::default())
  }

  /// Creates new pointwise metrics with custom configuration
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

/// Configuration parameters for pointwise metric computation
#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
pub struct PointwiseConfig {
  /// Floor value for relative error computation to prevent division by zero
  #[default = 1e-15]
  pub epsilon: f64,
  /// Threshold for log error computation (values below threshold are masked)
  #[default = 1e-10]
  pub log_threshold: f64,
  /// Tolerance for monotonicity violation detection
  #[default = 1e-12]
  pub monotonicity_eta: f64,
  /// Absolute tolerance thresholds for pass/fail mask [strict, moderate, loose]
  #[default(_code = "[1e-6, 1e-9, 1e-12]")]
  pub abs_tolerances: [f64; 3],
  /// Relative tolerance thresholds for pass/fail mask [strict, moderate, loose]
  #[default(_code = "[0.01, 0.001, 0.0001]")]
  pub rel_tolerances: [f64; 3],
}

/// Basic pointwise error metrics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PointwiseErrors {
  /// Absolute error: |actual - expected|
  pub absolute: Array1<f64>,
  /// Relative error with floor: |actual - expected| / max(|expected|, ε)
  pub relative: Array1<f64>,
  /// Signed error: actual - expected
  pub signed: Array1<f64>,
  /// Logarithmic error for tail analysis (masked where values below threshold)
  pub logarithmic: Array1<f64>,
  /// Summary statistics
  pub summary: PointwiseErrorSummary,
}

/// Summary statistics for pointwise errors
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PointwiseErrorSummary {
  pub abs_mean: f64,
  pub abs_max: f64,
  pub abs_std: f64,
  pub rel_mean: f64,
  pub rel_max: f64,
  pub rel_median: f64,
  pub signed_bias: f64,
  pub log_valid_count: usize,
  pub log_max: f64,
}

/// Derivative and structural error metrics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StructuralErrors {
  /// First derivative error using central differences
  pub first_derivative: Array1<f64>,
  /// Second derivative error using central differences
  pub second_derivative: Array1<f64>,
  /// Symmetry residual for symmetric distributions
  pub symmetry_residual: Array1<f64>,
  /// Monotonicity violation flags (1.0 where violation detected)
  pub monotonicity_violations: Array1<f64>,
  /// Summary statistics
  pub summary: StructuralErrorSummary,
}

/// Summary statistics for structural errors
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StructuralErrorSummary {
  pub d1_max: f64,
  pub d1_mean: f64,
  pub d2_max: f64,
  pub d2_mean: f64,
  pub symmetry_max: f64,
  pub monotonicity_violation_count: usize,
}

/// Tolerance compliance metrics at each point
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ToleranceMetrics {
  /// Pass/fail masks at three tolerance levels [strict, moderate, loose]
  pub pass_masks: [Array1<f64>; 3],
  /// Support coverage mask (1.0 where support mismatch detected)
  pub support_coverage_mask: Array1<f64>,
  /// Summary statistics
  pub summary: ToleranceSummary,
}

/// Summary statistics for tolerance metrics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ToleranceSummary {
  /// Pass fractions at each tolerance level
  pub pass_fractions: [f64; 3],
  /// Number of support coverage mismatches
  pub support_mismatch_count: usize,
}

// Implementation functions (extracted from pointwise_metrics.rs)
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

fn compute_pointwise_errors(
  actual: &Array1<f64>,
  expected: &Array1<f64>,
  config: &PointwiseConfig,
) -> eyre::Result<PointwiseErrors> {
  let n = actual.len();

  let signed = actual - expected;
  let absolute = signed.mapv(|x| x.abs());

  let mut relative = Array1::zeros(n);
  for i in 0..n {
    let denom = expected[i].abs().max(config.epsilon);
    relative[i] = signed[i].abs() / denom;
  }

  let mut logarithmic = Array1::from_elem(n, f64::NAN);
  let mut log_valid_count = 0;
  let mut log_max = 0.0;
  for i in 0..n {
    if actual[i] > config.log_threshold && expected[i] > config.log_threshold {
      let log_err = (actual[i].ln() - expected[i].ln()).abs();
      logarithmic[i] = log_err;
      log_valid_count += 1;
      if log_err > log_max {
        log_max = log_err;
      }
    }
  }

  let abs_mean = absolute.mean().unwrap_or(0.0);
  let abs_max = absolute.iter().copied().fold(0.0_f64, f64::max);
  let abs_variance = absolute.mapv(|x| (x - abs_mean).powi(2)).mean().unwrap_or(0.0);
  let abs_std = abs_variance.sqrt();

  let rel_mean = relative.mean().unwrap_or(0.0);
  let rel_max = relative.iter().copied().fold(0.0_f64, f64::max);
  let mut rel_sorted: Vec<f64> = relative.iter().copied().collect();
  rel_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
  let rel_median = compute_median(&rel_sorted);

  let signed_bias = signed.mean().unwrap_or(0.0);

  let summary = PointwiseErrorSummary {
    abs_mean,
    abs_max,
    abs_std,
    rel_mean,
    rel_max,
    rel_median,
    signed_bias,
    log_valid_count,
    log_max,
  };

  Ok(PointwiseErrors {
    absolute,
    relative,
    signed,
    logarithmic,
    summary,
  })
}

fn compute_structural_errors(
  x: &Array1<f64>,
  actual: &Array1<f64>,
  expected: &Array1<f64>,
  dx: f64,
  config: &PointwiseConfig,
) -> eyre::Result<StructuralErrors> {
  let n = x.len();

  let d1_actual = compute_first_derivative(actual, dx);
  let d1_expected = compute_first_derivative(expected, dx);
  let first_derivative = (&d1_actual - &d1_expected).mapv(|x| x.abs());

  let d2_actual = compute_second_derivative(actual, dx);
  let d2_expected = compute_second_derivative(expected, dx);
  let second_derivative = (&d2_actual - &d2_expected).mapv(|x| x.abs());

  let symmetry_residual = compute_symmetry_residual(x, actual);

  let mut monotonicity_violations = Array1::zeros(n);
  let mut violation_count = 0;
  for i in 0..n - 1 {
    if actual[i + 1] - actual[i] < -config.monotonicity_eta {
      monotonicity_violations[i] = 1.0;
      violation_count += 1;
    }
  }

  let d1_max = first_derivative.iter().copied().fold(0.0_f64, f64::max);
  let d1_mean = first_derivative.mean().unwrap_or(0.0);
  let d2_max = second_derivative.iter().copied().fold(0.0_f64, f64::max);
  let d2_mean = second_derivative.mean().unwrap_or(0.0);
  let symmetry_max = symmetry_residual.iter().copied().fold(0.0_f64, |a, b| a.max(b.abs()));

  let summary = StructuralErrorSummary {
    d1_max,
    d1_mean,
    d2_max,
    d2_mean,
    symmetry_max,
    monotonicity_violation_count: violation_count,
  };

  Ok(StructuralErrors {
    first_derivative,
    second_derivative,
    symmetry_residual,
    monotonicity_violations,
    summary,
  })
}

fn compute_tolerance_metrics(
  actual: &Array1<f64>,
  expected: &Array1<f64>,
  config: &PointwiseConfig,
) -> eyre::Result<ToleranceMetrics> {
  let n = actual.len();

  // Compute basic errors for tolerance checking
  let abs_errors = (actual - expected).mapv(|x| x.abs());
  let mut rel_errors = Array1::zeros(n);
  for i in 0..n {
    let denom = expected[i].abs().max(config.epsilon);
    rel_errors[i] = abs_errors[i] / denom;
  }

  // Compute tolerance pass masks
  let mut pass_masks = [Array1::zeros(n), Array1::zeros(n), Array1::zeros(n)];
  let mut pass_counts = [0_usize; 3];

  for level in 0..3 {
    let abs_tol = config.abs_tolerances[level];
    let rel_tol = config.rel_tolerances[level];
    for i in 0..n {
      let abs_pass = abs_errors[i] <= abs_tol;
      let rel_pass = rel_errors[i] <= rel_tol;
      if abs_pass || rel_pass {
        pass_masks[level][i] = 1.0;
        pass_counts[level] += 1;
      }
    }
  }

  // Compute support coverage mask
  let mut support_coverage_mask = Array1::zeros(n);
  let mut support_mismatch_count = 0;
  for i in 0..n {
    let actual_significant = actual[i].abs() > config.log_threshold;
    let expected_significant = expected[i].abs() > config.log_threshold;
    if actual_significant != expected_significant {
      support_coverage_mask[i] = 1.0;
      support_mismatch_count += 1;
    }
  }

  let pass_fractions = [
    pass_counts[0] as f64 / n as f64,
    pass_counts[1] as f64 / n as f64,
    pass_counts[2] as f64 / n as f64,
  ];

  let summary = ToleranceSummary {
    pass_fractions,
    support_mismatch_count,
  };

  Ok(ToleranceMetrics {
    pass_masks,
    support_coverage_mask,
    summary,
  })
}

// Helper functions for derivative computation
fn compute_first_derivative(y: &Array1<f64>, dx: f64) -> Array1<f64> {
  let n = y.len();
  let mut dy = Array1::zeros(n);

  if n < 2 {
    return dy;
  }

  dy[0] = (y[1] - y[0]) / dx;
  for i in 1..n - 1 {
    dy[i] = (y[i + 1] - y[i - 1]) / (2.0 * dx);
  }
  dy[n - 1] = (y[n - 1] - y[n - 2]) / dx;

  dy
}

fn compute_second_derivative(y: &Array1<f64>, dx: f64) -> Array1<f64> {
  let n = y.len();
  let mut d2y = Array1::zeros(n);

  if n < 3 {
    return d2y;
  }

  let dx2 = dx * dx;

  d2y[0] = (y[2] - 2.0 * y[1] + y[0]) / dx2;
  for i in 1..n - 1 {
    d2y[i] = (y[i + 1] - 2.0 * y[i] + y[i - 1]) / dx2;
  }
  d2y[n - 1] = (y[n - 1] - 2.0 * y[n - 2] + y[n - 3]) / dx2;

  d2y
}

fn compute_symmetry_residual(x: &Array1<f64>, y: &Array1<f64>) -> Array1<f64> {
  let n = x.len();
  let mut residual = Array1::zeros(n);

  let x_map: BTreeMap<String, (usize, f64)> = x
    .iter()
    .enumerate()
    .map(|(i, &val)| (format!("{val:.12}"), (i, val)))
    .collect();

  for (i, &xi) in x.iter().enumerate() {
    let neg_xi = -xi;
    let key = format!("{neg_xi:.12}");
    if let Some(&(j, _)) = x_map.get(&key) {
      residual[i] = y[i] - y[j];
    }
  }

  residual
}

fn compute_median(sorted: &[f64]) -> f64 {
  if sorted.is_empty() {
    return 0.0;
  }
  let n = sorted.len();
  if n.is_multiple_of(2) {
    f64::midpoint(sorted[n / 2 - 1], sorted[n / 2])
  } else {
    sorted[n / 2]
  }
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
