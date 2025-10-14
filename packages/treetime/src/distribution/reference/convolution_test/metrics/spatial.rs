use ndarray::Array1;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

/// Spatial and regional analysis metrics
///
/// These metrics analyze spatial patterns, regional behavior, and windowed aggregates.
/// They maintain arrays aligned with the evaluation grid but represent spatial
/// relationships and regional properties rather than pure pointwise values.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpatialMetrics {
  /// Total number of evaluation points
  pub total_points: usize,
  /// Grid spacing
  pub dx: f64,
  /// Regional analysis metrics
  pub regional: RegionalMetrics,
  /// Sliding window analysis
  pub windowed: WindowedMetrics,
  /// Cumulative integration metrics
  pub cumulative: CumulativeMetrics,
}

impl SpatialMetrics {
  /// Creates new spatial metrics from evaluation grid and function values
  pub fn new(
    x: &Array1<f64>,
    actual: &Array1<f64>,
    expected: &Array1<f64>,
    dx: f64,
  ) -> eyre::Result<Self> {
    Self::new_with_config(x, actual, expected, dx, &SpatialConfig::default())
  }

  /// Creates new spatial metrics with custom configuration
  pub fn new_with_config(
    x: &Array1<f64>,
    actual: &Array1<f64>,
    expected: &Array1<f64>,
    dx: f64,
    config: &SpatialConfig,
  ) -> eyre::Result<Self> {
    let n = x.len();

    let regional = compute_regional_metrics(x, actual, expected, config)?;
    let windowed = compute_windowed_metrics(actual, expected, config)?;
    let cumulative = compute_cumulative_metrics(actual, expected, dx)?;

    Ok(Self {
      total_points: n,
      dx,
      regional,
      windowed,
      cumulative,
    })
  }
}

/// Configuration parameters for spatial metric computation
#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
pub struct SpatialConfig {
  /// Floor value for relative error computation
  #[default = 1e-15]
  pub epsilon: f64,
  /// Threshold for log error computation
  #[default = 1e-10]
  pub log_threshold: f64,
  /// Peak region radius as multiple of effective width
  #[default = 3.0]
  pub peak_region_radius: f64,
  /// Tail threshold as fraction of peak value
  #[default = 1e-6]
  pub tail_threshold: f64,
  /// Sliding window half-width for spatial metrics (number of points)
  #[default = 5]
  pub window_half_width: usize,
}

/// Regional analysis metrics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RegionalMetrics {
  /// Peak region relative errors (masked outside peak region)
  pub peak_region_errors: Array1<f64>,
  /// Tail region log errors (masked outside tail regions)
  pub tail_region_errors: Array1<f64>,
  /// Summary statistics
  pub summary: RegionalSummary,
}

/// Summary statistics for regional metrics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RegionalSummary {
  /// Peak region statistics
  pub peak_region: RegionStats,
  /// Tail region statistics
  pub tail_region: RegionStats,
}

/// Statistics for a specific region
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RegionStats {
  /// Number of points in region
  pub point_count: usize,
  /// Mean error in region
  pub mean_error: f64,
  /// Maximum error in region
  pub max_error: f64,
  /// Standard deviation of errors in region
  pub std_error: f64,
}

/// Sliding window analysis metrics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WindowedMetrics {
  /// Sliding window RMS error profile
  pub sliding_rms: Array1<f64>,
  /// Sliding window maximum error profile
  pub sliding_max: Array1<f64>,
  /// Summary statistics
  pub summary: WindowedSummary,
}

/// Summary statistics for windowed metrics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WindowedSummary {
  /// Maximum sliding RMS error
  pub sliding_rms_max: f64,
  /// Maximum sliding max error
  pub sliding_max_max: f64,
}

/// Cumulative integration metrics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CumulativeMetrics {
  /// Cumulative error curve: ∫ (actual - expected) dx
  pub cumulative_error: Array1<f64>,
  /// Summary statistics
  pub summary: CumulativeSummary,
}

/// Summary statistics for cumulative metrics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CumulativeSummary {
  /// Final cumulative error value
  pub final_value: f64,
  /// Maximum absolute cumulative error
  pub max_abs: f64,
}

// Implementation functions

fn compute_regional_metrics(
  x: &Array1<f64>,
  actual: &Array1<f64>,
  expected: &Array1<f64>,
  config: &SpatialConfig,
) -> eyre::Result<RegionalMetrics> {
  let n = x.len();
  let errors = actual - expected;

  // Find peak position
  let peak_idx_expected = expected
    .iter()
    .enumerate()
    .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
    .map_or(0, |(i, _)| i);
  let peak_x = x[peak_idx_expected];

  let peak_val = expected.iter().copied().fold(f64::NEG_INFINITY, f64::max);
  let tail_threshold_val = peak_val * config.tail_threshold;

  // Peak region analysis
  let mut peak_region_errors = Array1::from_elem(n, f64::NAN);
  let mut peak_errors_valid = Vec::new();
  for i in 0..n {
    if (x[i] - peak_x).abs() <= config.peak_region_radius {
      let rel_err = errors[i].abs() / expected[i].abs().max(config.epsilon);
      peak_region_errors[i] = rel_err;
      peak_errors_valid.push(rel_err);
    }
  }

  // Tail region analysis
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

fn compute_windowed_metrics(
  actual: &Array1<f64>,
  expected: &Array1<f64>,
  config: &SpatialConfig,
) -> eyre::Result<WindowedMetrics> {
  let errors = actual - expected;
  let abs_errors = errors.mapv(|x| x.abs());

  let sliding_rms = compute_sliding_window_rms(&errors, config.window_half_width);
  let sliding_max = compute_sliding_window_max(&abs_errors, config.window_half_width);

  let sliding_rms_max = sliding_rms.iter().copied().fold(0.0_f64, f64::max);
  let sliding_max_max = sliding_max.iter().copied().fold(0.0_f64, f64::max);

  let summary = WindowedSummary {
    sliding_rms_max,
    sliding_max_max,
  };

  Ok(WindowedMetrics {
    sliding_rms,
    sliding_max,
    summary,
  })
}

fn compute_cumulative_metrics(
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

// Helper functions

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

fn compute_sliding_window_rms(errors: &Array1<f64>, half_width: usize) -> Array1<f64> {
  let n = errors.len();
  let mut rms = Array1::zeros(n);

  for i in 0..n {
    let start = i.saturating_sub(half_width);
    let end = (i + half_width + 1).min(n);
    let window = &errors.slice(ndarray::s![start..end]);
    let sum_sq: f64 = window.iter().map(|&x| x * x).sum();
    let count = window.len() as f64;
    rms[i] = (sum_sq / count).sqrt();
  }

  rms
}

fn compute_sliding_window_max(abs_errors: &Array1<f64>, half_width: usize) -> Array1<f64> {
  let n = abs_errors.len();
  let mut max_err = Array1::zeros(n);

  for i in 0..n {
    let start = i.saturating_sub(half_width);
    let end = (i + half_width + 1).min(n);
    let window = &abs_errors.slice(ndarray::s![start..end]);
    max_err[i] = window.iter().copied().fold(0.0_f64, f64::max);
  }

  max_err
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
    let result = SpatialMetrics::new(&x, &y, &y, 1.0).unwrap();

    assert_eq!(result.total_points, 5);
    assert_abs_diff_eq!(result.cumulative.summary.final_value, 0.0, epsilon = 1e-12);
    assert_abs_diff_eq!(result.cumulative.summary.max_abs, 0.0, epsilon = 1e-12);
  }

  #[test]
  fn test_cumulative_error_systematic_bias() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let expected = array![0.0, 0.0, 0.0, 0.0, 0.0];
    let actual = array![0.1, 0.1, 0.1, 0.1, 0.1];

    let result = SpatialMetrics::new(&x, &actual, &expected, 1.0).unwrap();

    assert!(result.cumulative.summary.final_value > 0.3);
  }
}