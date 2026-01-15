use crate::testing::metrics::config::SpatialConfig;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use treetime_utils::serde::{array1_as_vec, array1_from_vec};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RegionalMetrics {
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub peak_region_errors: Array1<f64>,
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub tail_region_errors: Array1<f64>,
  pub summary: RegionalSummary,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RegionalSummary {
  pub peak_region: RegionStats,
  pub tail_region: RegionStats,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RegionStats {
  pub point_count: usize,
  pub mean_error: f64,
  pub max_error: f64,
  pub std_error: f64,
}

pub(super) fn compute_regional_metrics(
  x: &Array1<f64>,
  actual: &Array1<f64>,
  expected: &Array1<f64>,
  config: &SpatialConfig,
) -> eyre::Result<RegionalMetrics> {
  let n = x.len();
  let errors = actual - expected;

  let peak_idx_expected = expected
    .iter()
    .enumerate()
    .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
    .map_or(0, |(i, _)| i);
  let peak_x = x[peak_idx_expected];

  let peak_val = expected.iter().copied().fold(f64::NEG_INFINITY, f64::max);
  let tail_threshold_val = peak_val * config.tail_threshold;

  let mut peak_region_errors = Array1::from_elem(n, f64::NAN);
  let mut peak_errors_valid = Vec::new();
  for i in 0..n {
    if (x[i] - peak_x).abs() <= config.peak_region_radius {
      let rel_err = errors[i].abs() / expected[i].abs().max(config.epsilon);
      peak_region_errors[i] = rel_err;
      peak_errors_valid.push(rel_err);
    }
  }

  let mut tail_region_errors = Array1::from_elem(n, f64::NAN);
  let mut tail_errors_valid = Vec::new();
  for i in 0..n {
    if expected[i] <= tail_threshold_val && actual[i] > config.log_threshold && expected[i] > config.log_threshold {
      let log_err = (actual[i].ln() - expected[i].ln()).abs();
      tail_region_errors[i] = log_err;
      tail_errors_valid.push(log_err);
    }
  }

  let peak_region = compute_region_stats(&peak_errors_valid);
  let tail_region = compute_region_stats(&tail_errors_valid);

  let summary = RegionalSummary {
    peak_region,
    tail_region,
  };

  Ok(RegionalMetrics {
    peak_region_errors,
    tail_region_errors,
    summary,
  })
}

fn compute_region_stats(errors: &[f64]) -> RegionStats {
  if errors.is_empty() {
    return RegionStats {
      point_count: 0,
      mean_error: 0.0,
      max_error: 0.0,
      std_error: 0.0,
    };
  }

  let count = errors.len();
  let mean = errors.iter().sum::<f64>() / count as f64;
  let max = errors.iter().copied().fold(0.0_f64, f64::max);
  let variance = errors.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / count as f64;
  let std = variance.sqrt();

  RegionStats {
    point_count: count,
    mean_error: mean,
    max_error: max,
    std_error: std,
  }
}
