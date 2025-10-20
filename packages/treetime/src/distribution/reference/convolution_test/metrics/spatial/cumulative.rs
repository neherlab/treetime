use ndarray::Array1;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CumulativeMetrics {
  pub cumulative_error: Array1<f64>,
  pub summary: CumulativeSummary,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CumulativeSummary {
  pub final_value: f64,
  pub max_abs: f64,
}

pub(super) fn compute_cumulative_metrics(
  actual: &Array1<f64>,
  expected: &Array1<f64>,
  dx: f64,
) -> eyre::Result<CumulativeMetrics> {
  let n = actual.len();
  let errors = actual - expected;

  let mut cumulative_error = Array1::zeros(n);
  cumulative_error[0] = errors[0] * dx;
  for i in 1..n {
    cumulative_error[i] = cumulative_error[i - 1] + errors[i] * dx;
  }

  let final_value = cumulative_error[n - 1];
  let max_abs = cumulative_error.iter().copied().fold(0.0_f64, |a, b| a.max(b.abs()));

  let summary = CumulativeSummary { final_value, max_abs };

  Ok(CumulativeMetrics {
    cumulative_error,
    summary,
  })
}
