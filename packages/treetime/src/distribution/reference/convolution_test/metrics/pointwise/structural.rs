use crate::distribution::reference::convolution_test::metrics::config::PointwiseConfig;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StructuralErrors {
  pub first_derivative: Array1<f64>,
  pub second_derivative: Array1<f64>,
  pub symmetry_residual: Array1<f64>,
  pub monotonicity_violations: Array1<f64>,
  pub summary: StructuralErrorSummary,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StructuralErrorSummary {
  pub d1_max: f64,
  pub d1_mean: f64,
  pub d2_max: f64,
  pub d2_mean: f64,
  pub symmetry_max: f64,
  pub monotonicity_violation_count: usize,
}

pub(super) fn compute_structural_errors(
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
