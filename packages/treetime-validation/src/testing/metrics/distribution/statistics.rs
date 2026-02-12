use crate::testing::metrics::pointwise::errors::PointwiseErrors;
use ndarray::Array1;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StatisticalMetrics {
  pub abs_error_stats: ErrorStatistics,
  pub rel_error_stats: ErrorStatistics,
  pub signed_error_stats: ErrorStatistics,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ErrorStatistics {
  pub mean: f64,
  pub std: f64,
  pub median: f64,
  pub q25: f64,
  pub q75: f64,
  pub q95: f64,
  pub q99: f64,
  pub skewness: f64,
  pub kurtosis: f64,
}

pub(super) fn compute_statistical_metrics(pointwise_errors: &PointwiseErrors) -> eyre::Result<StatisticalMetrics> {
  let abs_error_stats = compute_error_statistics(&pointwise_errors.absolute)?;
  let rel_error_stats = compute_error_statistics(&pointwise_errors.relative)?;
  let signed_error_stats = compute_error_statistics(&pointwise_errors.signed)?;

  Ok(StatisticalMetrics {
    abs_error_stats,
    rel_error_stats,
    signed_error_stats,
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

pub(super) fn compute_std(data: &Array1<f64>) -> f64 {
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

pub(super) fn compute_quantile(sorted: &[f64], q: f64) -> f64 {
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
  kurtosis_raw - correction
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_ulps_eq;
  use ndarray::array;

  #[test]
  fn test_statistical_metrics() {
    let errors = array![1.0, 2.0, 3.0, 4.0, 5.0];
    let stats = compute_error_statistics(&errors).unwrap();

    assert_ulps_eq!(stats.mean, 3.0, max_ulps = 4);
    assert_ulps_eq!(stats.median, 3.0, max_ulps = 4);
    assert!(stats.std > 0.0);
  }
}
