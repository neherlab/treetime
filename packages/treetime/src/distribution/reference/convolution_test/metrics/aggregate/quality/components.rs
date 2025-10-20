use serde::{Deserialize, Serialize};

use crate::distribution::reference::convolution_test::metrics::aggregate::domain_agreement::DomainAgreementMetrics;
use crate::distribution::reference::convolution_test::metrics::aggregate::performance::PerformanceMetrics;

/// Quality score components
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityComponents {
  pub accuracy_score: f64,
  pub precision_score: f64,
  pub stability_score: f64,
  pub robustness_score: f64,
  pub overall_score: f64,
}

pub fn compute_accuracy_score(domain_agreement: &DomainAgreementMetrics) -> f64 {
  let r_squared = domain_agreement.quality_metrics.r_squared;
  (r_squared * 100.0).clamp(0.0, 100.0)
}

pub fn compute_precision_score(performance: &PerformanceMetrics) -> f64 {
  let snr = performance.signal_to_noise_ratio;
  if snr > 60.0 {
    100.0
  } else if snr > 0.0 {
    (snr / 60.0 * 100.0).max(0.0)
  } else {
    0.0
  }
}

pub fn compute_stability_score(domain_agreement: &DomainAgreementMetrics) -> f64 {
  let max_error = domain_agreement.abs_error_stats.max;
  let mean_error = domain_agreement.abs_error_stats.mean;

  if mean_error > 0.0 {
    let stability_ratio = 1.0 - (max_error / (mean_error * 10.0)).min(1.0);
    (stability_ratio * 100.0).max(0.0_f64)
  } else {
    100.0
  }
}

pub fn compute_robustness_score(performance: &PerformanceMetrics) -> f64 {
  let mape = performance.mean_absolute_percentage_error;
  if mape < 1.0 {
    100.0
  } else if mape < 100.0 {
    (100.0 - mape).max(0.0)
  } else {
    0.0
  }
}
