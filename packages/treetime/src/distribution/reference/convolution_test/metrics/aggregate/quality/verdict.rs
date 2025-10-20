use serde::{Deserialize, Serialize};

use crate::distribution::reference::convolution_test::metrics::aggregate::domain_agreement::DomainAgreementMetrics;
use crate::distribution::reference::convolution_test::metrics::aggregate::performance::PerformanceMetrics;
use crate::distribution::reference::convolution_test::metrics::aggregate::quality::components::QualityComponents;
use crate::o;

/// Quality verdict with strengths, weaknesses, and recommendations
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityVerdict {
  pub strengths: Vec<String>,
  pub weaknesses: Vec<String>,
  pub recommendations: Vec<String>,
  pub critical_issues: Vec<String>,
}

pub fn generate_quality_verdict(
  domain_agreement: &DomainAgreementMetrics,
  performance: &PerformanceMetrics,
  components: &QualityComponents,
) -> QualityVerdict {
  let mut strengths = Vec::new();
  let mut weaknesses = Vec::new();
  let mut recommendations = Vec::new();
  let mut critical_issues = Vec::new();

  // Analyze strengths
  if components.accuracy_score > 95.0 {
    strengths.push(o!("Excellent accuracy"));
  }
  if performance.signal_to_noise_ratio > 60.0 {
    strengths.push(o!("Very high signal-to-noise ratio"));
  }
  if domain_agreement.quality_metrics.r_squared > 0.99 {
    strengths.push(o!("Near-perfect correlation"));
  }

  // Analyze weaknesses
  if components.accuracy_score < 70.0 {
    weaknesses.push(o!("Poor accuracy performance"));
  }
  if performance.normalized_rmse > 0.1 {
    weaknesses.push(o!("High normalized error"));
  }
  if domain_agreement.quality_metrics.r_squared < 0.9 {
    weaknesses.push(o!("Low correlation with reference"));
  }

  // Generate recommendations
  if components.overall_score > 90.0 {
    recommendations.push(o!("Suitable for production use"));
    recommendations.push(o!("Can be used for high-precision applications"));
  } else if components.overall_score > 70.0 {
    recommendations.push(o!(
      "Suitable for most applications with moderate precision requirements"
    ));
  } else {
    recommendations.push(o!("Requires improvement before production use"));
    recommendations.push(o!("Consider parameter tuning or algorithm alternatives"));
  }

  // Identify critical issues
  if domain_agreement.quality_metrics.r_squared < 0.5 {
    critical_issues.push(o!("Very poor correlation - fundamental algorithm issues"));
  }
  if performance.mean_absolute_percentage_error > 50.0 {
    critical_issues.push(o!("Extremely high percentage errors"));
  }

  QualityVerdict {
    strengths,
    weaknesses,
    recommendations,
    critical_issues,
  }
}
