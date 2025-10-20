use crate::distribution::reference::convolution_test::metrics::config::PointwiseConfig;
use ndarray::Array1;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ToleranceMetrics {
  pub pass_masks: [Array1<f64>; 3],
  pub support_coverage_mask: Array1<f64>,
  pub summary: ToleranceSummary,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ToleranceSummary {
  pub pass_fractions: [f64; 3],
  pub support_mismatch_count: usize,
}

pub(super) fn compute_tolerance_metrics(
  actual: &Array1<f64>,
  expected: &Array1<f64>,
  config: &PointwiseConfig,
) -> eyre::Result<ToleranceMetrics> {
  let n = actual.len();

  let abs_errors = (actual - expected).mapv(|x| x.abs());
  let mut rel_errors = Array1::zeros(n);
  for i in 0..n {
    let denom = expected[i].abs().max(config.epsilon);
    rel_errors[i] = abs_errors[i] / denom;
  }

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
