use crate::testing::metrics::config::DistributionConfig;
use crate::testing::metrics::pointwise::errors::PointwiseErrors;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use treetime_utils::make_error;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HistogramMode {
  Unsigned { use_log_scale: bool },
  Signed,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HistogramMetrics {
  pub abs_error_histogram: ErrorHistogram,
  pub rel_error_histogram: ErrorHistogram,
  pub signed_error_histogram: ErrorHistogram,
  pub summary: HistogramSummary,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ErrorHistogram {
  pub bin_edges: Vec<f64>,
  pub bin_counts: Vec<usize>,
  pub bin_centers: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HistogramSummary {
  pub total_count: usize,
  pub modal_range: (f64, f64),
  pub spread_measure: f64,
}

pub(super) fn compute_histogram_metrics(
  pointwise_errors: &PointwiseErrors,
  config: &DistributionConfig,
) -> eyre::Result<HistogramMetrics> {
  let num_bins = config.histogram_bins;
  let abs_error_histogram = compute_error_histogram(
    &pointwise_errors.absolute,
    num_bins,
    HistogramMode::Unsigned { use_log_scale: true },
  )?;
  let rel_error_histogram = compute_error_histogram(
    &pointwise_errors.relative,
    num_bins,
    HistogramMode::Unsigned { use_log_scale: false },
  )?;
  let signed_error_histogram = compute_error_histogram(&pointwise_errors.signed, num_bins, HistogramMode::Signed)?;

  let total_count = abs_error_histogram.bin_counts.iter().sum();

  let max_bin_idx = abs_error_histogram
    .bin_counts
    .iter()
    .enumerate()
    .max_by_key(|&(_, count)| count)
    .map_or(0, |(idx, _)| idx);

  let modal_range = if max_bin_idx < abs_error_histogram.bin_edges.len() - 1 {
    (
      abs_error_histogram.bin_edges[max_bin_idx],
      abs_error_histogram.bin_edges[max_bin_idx + 1],
    )
  } else {
    (0.0, 0.0)
  };

  let spread_measure = if abs_error_histogram.bin_edges.len() >= 2 {
    abs_error_histogram.bin_edges.last().unwrap() - abs_error_histogram.bin_edges.first().unwrap()
  } else {
    0.0
  };

  let summary = HistogramSummary {
    total_count,
    modal_range,
    spread_measure,
  };

  Ok(HistogramMetrics {
    abs_error_histogram,
    rel_error_histogram,
    signed_error_histogram,
    summary,
  })
}

fn compute_error_histogram(errors: &Array1<f64>, num_bins: usize, mode: HistogramMode) -> eyre::Result<ErrorHistogram> {
  if num_bins == 0 {
    return make_error!("Number of histogram bins must be positive");
  }

  let use_log_scale = matches!(mode, HistogramMode::Unsigned { use_log_scale: true });

  let valid_errors: Vec<f64> = match mode {
    HistogramMode::Unsigned { .. } => errors.iter().copied().filter(|&x| x > 0.0 && x.is_finite()).collect(),
    HistogramMode::Signed => errors.iter().copied().filter(|&x| x.is_finite()).collect(),
  };

  if valid_errors.is_empty() {
    return Ok(ErrorHistogram {
      bin_edges: vec![0.0; num_bins + 1],
      bin_counts: vec![0; num_bins],
      bin_centers: vec![0.0; num_bins],
    });
  }

  let (min_val, max_val) = if use_log_scale {
    let min_log = valid_errors
      .iter()
      .copied()
      .map(|x| x.log10())
      .fold(f64::INFINITY, f64::min);
    let max_log = valid_errors
      .iter()
      .copied()
      .map(|x| x.log10())
      .fold(f64::NEG_INFINITY, f64::max);
    (min_log, max_log)
  } else {
    let min_val = valid_errors.iter().copied().fold(f64::INFINITY, f64::min);
    let max_val = valid_errors.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    (min_val, max_val)
  };

  let range = max_val - min_val;
  let bin_width = if range > 0.0 { range / num_bins as f64 } else { 1.0 };

  let mut bin_edges = Vec::with_capacity(num_bins + 1);
  for i in 0..=num_bins {
    bin_edges.push(min_val + i as f64 * bin_width);
  }

  let mut bin_counts = vec![0_usize; num_bins];
  for &err in &valid_errors {
    let val = if use_log_scale { err.log10() } else { err };
    let bin_idx = ((val - min_val) / bin_width).floor() as usize;
    let bin_idx = bin_idx.min(num_bins - 1);
    bin_counts[bin_idx] += 1;
  }

  let bin_centers: Vec<f64> = (0..num_bins)
    .map(|i| f64::midpoint(bin_edges[i], bin_edges[i + 1]))
    .collect();

  Ok(ErrorHistogram {
    bin_edges,
    bin_counts,
    bin_centers,
  })
}

#[cfg(test)]
mod tests {
  use super::*;
  use ndarray::{Array1, array};

  #[test]
  fn test_histograms_unsigned_log_scale() {
    let errors = Array1::linspace(1e-6, 1e-3, 100);
    let histogram = compute_error_histogram(&errors, 10, HistogramMode::Unsigned { use_log_scale: true }).unwrap();

    assert_eq!(histogram.bin_counts.len(), 10);
    assert_eq!(histogram.bin_edges.len(), 11);
    assert_eq!(histogram.bin_centers.len(), 10);

    let total_count: usize = histogram.bin_counts.iter().sum();
    assert_eq!(total_count, 100);
  }

  #[test]
  fn test_histograms_unsigned_linear_scale() {
    let errors = Array1::linspace(0.01, 1.0, 50);
    let histogram = compute_error_histogram(&errors, 5, HistogramMode::Unsigned { use_log_scale: false }).unwrap();

    assert_eq!(histogram.bin_counts.len(), 5);
    assert_eq!(histogram.bin_edges.len(), 6);

    let total_count: usize = histogram.bin_counts.iter().sum();
    assert_eq!(total_count, 50);
  }

  #[test]
  fn test_histograms_signed() {
    let errors = array![-2.0, -1.0, 0.0, 1.0, 2.0];
    let histogram = compute_error_histogram(&errors, 4, HistogramMode::Signed).unwrap();

    assert_eq!(histogram.bin_counts.len(), 4);
    assert_eq!(histogram.bin_edges.len(), 5);

    let total_count: usize = histogram.bin_counts.iter().sum();
    assert_eq!(total_count, 5);
  }

  #[test]
  fn test_histograms_unsigned_filters_non_positive() {
    let errors = array![-1.0, 0.0, 1.0, 2.0, f64::NAN, f64::INFINITY];
    let histogram = compute_error_histogram(&errors, 2, HistogramMode::Unsigned { use_log_scale: false }).unwrap();

    let total_count: usize = histogram.bin_counts.iter().sum();
    assert_eq!(total_count, 2);
  }

  #[test]
  fn test_histograms_signed_filters_non_finite() {
    let errors = array![-1.0, 0.0, 1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY];
    let histogram = compute_error_histogram(&errors, 3, HistogramMode::Signed).unwrap();

    let total_count: usize = histogram.bin_counts.iter().sum();
    assert_eq!(total_count, 3);
  }

  #[test]
  fn test_histograms_zero_bins_error() {
    let errors = array![1.0, 2.0];
    let result = compute_error_histogram(&errors, 0, HistogramMode::Signed);
    drop(result.unwrap_err());
  }

  #[test]
  fn test_histograms_empty_valid_errors() {
    let errors = array![f64::NAN, f64::INFINITY];
    let histogram = compute_error_histogram(&errors, 5, HistogramMode::Signed).unwrap();

    assert_eq!(histogram.bin_counts, vec![0; 5]);
    assert_eq!(histogram.bin_edges, vec![0.0; 6]);
  }
}
