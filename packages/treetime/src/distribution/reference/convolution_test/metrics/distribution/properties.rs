use crate::distribution::reference::convolution_test::metrics::pointwise::PointwiseErrors;
use ndarray::Array1;
use serde::{Deserialize, Serialize};

use super::statistics::{compute_quantile, compute_std};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DistributionProperties {
  pub dynamic_range: f64,
  pub symmetry_measure: f64,
  pub tail_behavior: TailBehavior,
  pub outlier_stats: OutlierStatistics,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum TailBehavior {
  Light,
  Heavy,
  ExtremelyHeavy,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutlierStatistics {
  pub count_3sigma: usize,
  pub count_5sigma: usize,
  pub fraction_outliers: f64,
}

pub(super) fn compute_distribution_properties(
  pointwise_errors: &PointwiseErrors,
) -> eyre::Result<DistributionProperties> {
  let abs_errors = &pointwise_errors.absolute;

  let min_val = abs_errors.iter().copied().fold(f64::INFINITY, f64::min);
  let max_val = abs_errors.iter().copied().fold(0.0_f64, f64::max);
  let dynamic_range = if min_val > 0.0 && min_val.is_finite() {
    (max_val / min_val).log10()
  } else {
    f64::INFINITY
  };

  let signed_errors = &pointwise_errors.signed;
  let mean_signed = signed_errors.mean().unwrap_or(0.0);
  let abs_mean_signed = mean_signed.abs();
  let std_signed = compute_std(signed_errors);
  let symmetry_measure = if std_signed > 0.0 {
    abs_mean_signed / std_signed
  } else {
    0.0
  };

  let tail_behavior = classify_tail_behavior(abs_errors);

  let outlier_stats = compute_outlier_statistics(abs_errors);

  Ok(DistributionProperties {
    dynamic_range,
    symmetry_measure,
    tail_behavior,
    outlier_stats,
  })
}

fn classify_tail_behavior(errors: &Array1<f64>) -> TailBehavior {
  let mut sorted_errors: Vec<f64> = errors.iter().copied().filter(|&x| x > 0.0 && x.is_finite()).collect();
  sorted_errors.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

  if sorted_errors.len() < 10 {
    return TailBehavior::Light;
  }

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
