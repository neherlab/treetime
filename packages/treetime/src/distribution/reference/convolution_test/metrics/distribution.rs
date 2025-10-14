use crate::make_error;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

/// Statistical distribution analysis metrics
///
/// These metrics analyze the statistical properties of error distributions
/// across the entire domain, providing histogram analysis, percentiles,
/// and distributional characteristics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DistributionMetrics {
  /// Error histogram analysis
  pub histograms: HistogramMetrics,
  /// Statistical percentiles and moments
  pub statistics: StatisticalMetrics,
  /// Distributional properties
  pub properties: DistributionProperties,
}

impl DistributionMetrics {
  /// Creates new distribution metrics from pointwise errors
  pub fn new(pointwise_errors: &PointwiseErrors) -> eyre::Result<Self> {
    Self::new_with_config(pointwise_errors, &DistributionConfig::default())
  }

  /// Creates new distribution metrics with custom configuration
  pub fn new_with_config(
    pointwise_errors: &PointwiseErrors,
    config: &DistributionConfig,
  ) -> eyre::Result<Self> {
    let histograms = compute_histogram_metrics(pointwise_errors, config)?;
    let statistics = compute_statistical_metrics(pointwise_errors)?;
    let properties = compute_distribution_properties(pointwise_errors)?;

    Ok(Self {
      histograms,
      statistics,
      properties,
    })
  }
}

/// Configuration parameters for distribution metric computation
#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
pub struct DistributionConfig {
  /// Number of bins for error histogram
  #[default = 50]
  pub histogram_bins: usize,
  /// Minimum value for log scale histograms
  #[default = 1e-16]
  pub log_min_value: f64,
}

/// Error histogram analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HistogramMetrics {
  /// Histogram of log10(absolute errors)
  pub abs_error_histogram: ErrorHistogram,
  /// Histogram of relative errors
  pub rel_error_histogram: ErrorHistogram,
  /// Histogram of signed errors
  pub signed_error_histogram: ErrorHistogram,
  /// Summary statistics
  pub summary: HistogramSummary,
}

/// Histogram representation for error distributions
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ErrorHistogram {
  /// Bin edges (length = bin_counts.len() + 1)
  pub bin_edges: Vec<f64>,
  /// Count of errors in each bin
  pub bin_counts: Vec<usize>,
  /// Bin centers for plotting
  pub bin_centers: Vec<f64>,
}

/// Summary statistics for histograms
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HistogramSummary {
  /// Total count across all histograms
  pub total_count: usize,
  /// Most frequent error range (mode)
  pub modal_range: (f64, f64),
  /// Histogram spread measure
  pub spread_measure: f64,
}

/// Statistical percentiles and moments
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StatisticalMetrics {
  /// Absolute error statistics
  pub abs_error_stats: ErrorStatistics,
  /// Relative error statistics
  pub rel_error_stats: ErrorStatistics,
  /// Signed error statistics
  pub signed_error_stats: ErrorStatistics,
}

/// Statistical measures for error distributions
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ErrorStatistics {
  /// Mean value
  pub mean: f64,
  /// Standard deviation
  pub std: f64,
  /// Median (50th percentile)
  pub median: f64,
  /// 25th percentile
  pub q25: f64,
  /// 75th percentile
  pub q75: f64,
  /// 95th percentile
  pub q95: f64,
  /// 99th percentile
  pub q99: f64,
  /// Skewness (asymmetry measure)
  pub skewness: f64,
  /// Kurtosis (tail heaviness measure)
  pub kurtosis: f64,
}

/// Distributional properties and characteristics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DistributionProperties {
  /// Effective dynamic range
  pub dynamic_range: f64,
  /// Distribution symmetry measure
  pub symmetry_measure: f64,
  /// Tail behavior classification
  pub tail_behavior: TailBehavior,
  /// Outlier statistics
  pub outlier_stats: OutlierStatistics,
}

/// Classification of tail behavior
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum TailBehavior {
  /// Light tails (exponential-like decay)
  Light,
  /// Heavy tails (polynomial-like decay)
  Heavy,
  /// Extremely heavy tails
  ExtremelyHeavy,
}

/// Statistics about outliers in the distribution
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutlierStatistics {
  /// Number of outliers (beyond 3 standard deviations)
  pub count_3sigma: usize,
  /// Number of extreme outliers (beyond 5 standard deviations)
  pub count_5sigma: usize,
  /// Fraction of outliers
  pub fraction_outliers: f64,
}

// Temporary import from pointwise.rs - this will be resolved when we integrate modules
use super::pointwise::PointwiseErrors;

// Implementation functions

fn compute_histogram_metrics(
  pointwise_errors: &PointwiseErrors,
  config: &DistributionConfig,
) -> eyre::Result<HistogramMetrics> {
  let abs_error_histogram = compute_error_histogram(&pointwise_errors.absolute, config.histogram_bins, true)?;
  let rel_error_histogram = compute_error_histogram(&pointwise_errors.relative, config.histogram_bins, false)?;
  let signed_error_histogram = compute_error_histogram_signed(&pointwise_errors.signed, config.histogram_bins)?;

  let total_count = abs_error_histogram.bin_counts.iter().sum();

  // Find modal range (bin with highest count)
  let max_bin_idx = abs_error_histogram
    .bin_counts
    .iter()
    .enumerate()
    .max_by_key(|&(_, count)| count)
    .map(|(idx, _)| idx)
    .unwrap_or(0);

  let modal_range = if max_bin_idx < abs_error_histogram.bin_edges.len() - 1 {
    (
      abs_error_histogram.bin_edges[max_bin_idx],
      abs_error_histogram.bin_edges[max_bin_idx + 1],
    )
  } else {
    (0.0, 0.0)
  };

  // Compute spread measure (interquartile range in log space)
  let spread_measure = if abs_error_histogram.bin_edges.len() >= 2 {
    let log_range = abs_error_histogram.bin_edges.last().unwrap()
      - abs_error_histogram.bin_edges.first().unwrap();
    log_range
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

fn compute_statistical_metrics(pointwise_errors: &PointwiseErrors) -> eyre::Result<StatisticalMetrics> {
  let abs_error_stats = compute_error_statistics(&pointwise_errors.absolute)?;
  let rel_error_stats = compute_error_statistics(&pointwise_errors.relative)?;
  let signed_error_stats = compute_error_statistics(&pointwise_errors.signed)?;

  Ok(StatisticalMetrics {
    abs_error_stats,
    rel_error_stats,
    signed_error_stats,
  })
}

fn compute_distribution_properties(pointwise_errors: &PointwiseErrors) -> eyre::Result<DistributionProperties> {
  let abs_errors = &pointwise_errors.absolute;

  // Dynamic range
  let min_val = abs_errors.iter().copied().fold(f64::INFINITY, f64::min);
  let max_val = abs_errors.iter().copied().fold(0.0_f64, f64::max);
  let dynamic_range = if min_val > 0.0 && min_val.is_finite() {
    (max_val / min_val).log10()
  } else {
    f64::INFINITY
  };

  // Symmetry measure (based on signed errors)
  let signed_errors = &pointwise_errors.signed;
  let mean_signed = signed_errors.mean().unwrap_or(0.0);
  let abs_mean_signed = mean_signed.abs();
  let std_signed = compute_std(signed_errors);
  let symmetry_measure = if std_signed > 0.0 {
    abs_mean_signed / std_signed
  } else {
    0.0
  };

  // Tail behavior classification
  let tail_behavior = classify_tail_behavior(abs_errors);

  // Outlier statistics
  let outlier_stats = compute_outlier_statistics(abs_errors);

  Ok(DistributionProperties {
    dynamic_range,
    symmetry_measure,
    tail_behavior,
    outlier_stats,
  })
}

// Helper functions

fn compute_error_histogram(errors: &Array1<f64>, num_bins: usize, use_log_scale: bool) -> eyre::Result<ErrorHistogram> {
  if num_bins == 0 {
    return make_error!("Number of histogram bins must be positive");
  }

  let valid_errors: Vec<f64> = errors.iter().copied().filter(|&x| x > 0.0 && x.is_finite()).collect();

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

fn compute_error_histogram_signed(errors: &Array1<f64>, num_bins: usize) -> eyre::Result<ErrorHistogram> {
  if num_bins == 0 {
    return make_error!("Number of histogram bins must be positive");
  }

  let valid_errors: Vec<f64> = errors.iter().copied().filter(|&x| x.is_finite()).collect();

  if valid_errors.is_empty() {
    return Ok(ErrorHistogram {
      bin_edges: vec![0.0; num_bins + 1],
      bin_counts: vec![0; num_bins],
      bin_centers: vec![0.0; num_bins],
    });
  }

  let min_val = valid_errors.iter().copied().fold(f64::INFINITY, f64::min);
  let max_val = valid_errors.iter().copied().fold(f64::NEG_INFINITY, f64::max);

  let range = max_val - min_val;
  let bin_width = if range > 0.0 { range / num_bins as f64 } else { 1.0 };

  let mut bin_edges = Vec::with_capacity(num_bins + 1);
  for i in 0..=num_bins {
    bin_edges.push(min_val + i as f64 * bin_width);
  }

  let mut bin_counts = vec![0_usize; num_bins];
  for &err in &valid_errors {
    let bin_idx = ((err - min_val) / bin_width).floor() as usize;
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

fn compute_error_statistics(errors: &Array1<f64>) -> eyre::Result<ErrorStatistics> {
  let mut sorted_errors: Vec<f64> = errors.iter().copied().filter(|x| x.is_finite()).collect();
  sorted_errors.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

  if sorted_errors.is_empty() {
    return Ok(ErrorStatistics {
      mean: 0.0,
      std: 0.0,
      median: 0.0,
      q25: 0.0,
      q75: 0.0,
      q95: 0.0,
      q99: 0.0,
      skewness: 0.0,
      kurtosis: 0.0,
    });
  }

  let mean = sorted_errors.iter().sum::<f64>() / sorted_errors.len() as f64;
  let std = compute_std_from_sorted(&sorted_errors, mean);

  let median = compute_quantile(&sorted_errors, 0.5);
  let q25 = compute_quantile(&sorted_errors, 0.25);
  let q75 = compute_quantile(&sorted_errors, 0.75);
  let q95 = compute_quantile(&sorted_errors, 0.95);
  let q99 = compute_quantile(&sorted_errors, 0.99);

  let skewness = compute_skewness(&sorted_errors, mean, std);
  let kurtosis = compute_kurtosis(&sorted_errors, mean, std);

  Ok(ErrorStatistics {
    mean,
    std,
    median,
    q25,
    q75,
    q95,
    q99,
    skewness,
    kurtosis,
  })
}

fn classify_tail_behavior(errors: &Array1<f64>) -> TailBehavior {
  let mut sorted_errors: Vec<f64> = errors.iter().copied().filter(|&x| x > 0.0 && x.is_finite()).collect();
  sorted_errors.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

  if sorted_errors.len() < 10 {
    return TailBehavior::Light;
  }

  // Compare 95th vs 99th percentile ratio
  let q95 = compute_quantile(&sorted_errors, 0.95);
  let q99 = compute_quantile(&sorted_errors, 0.99);

  if q95 <= 0.0 {
    return TailBehavior::Light;
  }

  let tail_ratio = q99 / q95;

  if tail_ratio > 100.0 {
    TailBehavior::ExtremelyHeavy
  } else if tail_ratio > 10.0 {
    TailBehavior::Heavy
  } else {
    TailBehavior::Light
  }
}

fn compute_outlier_statistics(errors: &Array1<f64>) -> OutlierStatistics {
  let mean = errors.mean().unwrap_or(0.0);
  let std = compute_std(errors);

  if std <= 0.0 {
    return OutlierStatistics {
      count_3sigma: 0,
      count_5sigma: 0,
      fraction_outliers: 0.0,
    };
  }

  let threshold_3sigma = mean + 3.0 * std;
  let threshold_5sigma = mean + 5.0 * std;

  let count_3sigma = errors.iter().filter(|&&x| x.abs() > threshold_3sigma).count();
  let count_5sigma = errors.iter().filter(|&&x| x.abs() > threshold_5sigma).count();

  let fraction_outliers = count_3sigma as f64 / errors.len() as f64;

  OutlierStatistics {
    count_3sigma,
    count_5sigma,
    fraction_outliers,
  }
}

// Statistical helper functions

fn compute_std(data: &Array1<f64>) -> f64 {
  let mean = data.mean().unwrap_or(0.0);
  compute_std_from_sorted(&data.iter().copied().collect::<Vec<_>>(), mean)
}

fn compute_std_from_sorted(sorted_data: &[f64], mean: f64) -> f64 {
  if sorted_data.len() <= 1 {
    return 0.0;
  }
  let variance = sorted_data.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (sorted_data.len() - 1) as f64;
  variance.sqrt()
}

fn compute_quantile(sorted: &[f64], q: f64) -> f64 {
  if sorted.is_empty() {
    return 0.0;
  }
  let n = sorted.len();
  let index = (q * (n - 1) as f64).round() as usize;
  sorted[index.min(n - 1)]
}

fn compute_skewness(sorted_data: &[f64], mean: f64, std: f64) -> f64 {
  if std <= 0.0 || sorted_data.len() < 3 {
    return 0.0;
  }
  let n = sorted_data.len() as f64;
  let skew_sum = sorted_data.iter().map(|&x| ((x - mean) / std).powi(3)).sum::<f64>();
  (n / ((n - 1.0) * (n - 2.0))) * skew_sum
}

fn compute_kurtosis(sorted_data: &[f64], mean: f64, std: f64) -> f64 {
  if std <= 0.0 || sorted_data.len() < 4 {
    return 0.0;
  }
  let n = sorted_data.len() as f64;
  let kurt_sum = sorted_data.iter().map(|&x| ((x - mean) / std).powi(4)).sum::<f64>();
  let kurtosis_raw = (n * (n + 1.0) / ((n - 1.0) * (n - 2.0) * (n - 3.0))) * kurt_sum;
  let correction = 3.0 * (n - 1.0).powi(2) / ((n - 2.0) * (n - 3.0));
  kurtosis_raw - correction // Excess kurtosis (subtract 3 for normal distribution)
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_abs_diff_eq;
  use ndarray::array;

  #[test]
  fn test_histogram_creation() {
    let errors = Array1::linspace(1e-6, 1e-3, 100);
    let histogram = compute_error_histogram(&errors, 10, true).unwrap();

    assert_eq!(histogram.bin_counts.len(), 10);
    assert_eq!(histogram.bin_edges.len(), 11);
    assert_eq!(histogram.bin_centers.len(), 10);

    let total_count: usize = histogram.bin_counts.iter().sum();
    assert_eq!(total_count, 100);
  }

  #[test]
  fn test_statistical_metrics() {
    let errors = array![1.0, 2.0, 3.0, 4.0, 5.0];
    let stats = compute_error_statistics(&errors).unwrap();

    assert_abs_diff_eq!(stats.mean, 3.0, epsilon = 1e-12);
    assert_abs_diff_eq!(stats.median, 3.0, epsilon = 1e-12);
    assert!(stats.std > 0.0);
  }
}