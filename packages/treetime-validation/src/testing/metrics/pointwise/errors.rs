use crate::testing::metrics::config::PointwiseConfig;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use treetime_utils::array::serde::{array1_as_vec, array1_from_vec};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PointwiseErrors {
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub absolute: Array1<f64>,
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub relative: Array1<f64>,
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub signed: Array1<f64>,
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub logarithmic: Array1<f64>,
  pub summary: PointwiseErrorSummary,
}

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

pub(super) fn compute_pointwise_errors(
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
