use ndarray::{Array1, Axis};
use serde::{Deserialize, Serialize};
use treetime_utils::ndarray::cumsum_axis;

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
  let errors = actual - expected;
  let scaled_errors = &errors * dx;

  let cumulative_error = cumsum_axis(&scaled_errors, Axis(0));

  let final_value = *cumulative_error.last().unwrap_or(&0.0);
  let max_abs = cumulative_error.iter().copied().fold(0.0_f64, |a, b| a.max(b.abs()));

  let summary = CumulativeSummary { final_value, max_abs };

  Ok(CumulativeMetrics {
    cumulative_error,
    summary,
  })
}
