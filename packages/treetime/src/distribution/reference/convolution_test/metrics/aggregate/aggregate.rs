use crate::distribution::reference::convolution_test::metrics::aggregate::domain_agreement::domain_agreement::DomainAgreementMetrics;
use crate::distribution::reference::convolution_test::metrics::aggregate::efficiency::metrics::{
  EfficiencyMetrics, compute_efficiency_metrics,
};
use crate::distribution::reference::convolution_test::metrics::aggregate::performance::performance::{
  PerformanceMetrics, compute_performance_metrics,
};
use crate::distribution::reference::convolution_test::metrics::aggregate::quality::assessment::{
  QualityAssessment, compute_quality_assessment,
};
use crate::distribution::reference::convolution_test::metrics::config::ToleranceThresholds;
use eyre::Result;
use ndarray::Array1;
use serde::{Deserialize, Serialize};

/// Enhanced aggregate metrics combining domain agreement with comprehensive summaries
///
/// These metrics provide domain-wide scalar summaries that give an overall
/// assessment of algorithm performance across all evaluation points.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AggregateMetrics {
  /// Original domain agreement metrics (for backward compatibility)
  pub domain_agreement: DomainAgreementMetrics,
  /// Enhanced performance assessment
  pub performance: PerformanceMetrics,
  /// Quality assessment based on multiple criteria
  pub quality: QualityAssessment,
  /// Algorithm efficiency metrics
  pub efficiency: EfficiencyMetrics,
}

impl AggregateMetrics {
  pub fn new(x: &Array1<f64>, actual: &Array1<f64>, expected: &Array1<f64>, execution_time_ms: f64) -> Result<Self> {
    let domain_agreement = DomainAgreementMetrics::new(x, actual, expected)?;
    let performance = compute_performance_metrics(actual, expected)?;
    let quality = compute_quality_assessment(&domain_agreement, &performance)?;
    let efficiency = compute_efficiency_metrics(execution_time_ms, x.len())?;

    Ok(Self {
      domain_agreement,
      performance,
      quality,
      efficiency,
    })
  }

  pub fn new_with_thresholds(
    x: &Array1<f64>,
    actual: &Array1<f64>,
    expected: &Array1<f64>,
    execution_time_ms: f64,
    thresholds: &ToleranceThresholds,
  ) -> Result<Self> {
    let domain_agreement = DomainAgreementMetrics::new_with_thresholds(x, actual, expected, thresholds)?;
    let performance = compute_performance_metrics(actual, expected)?;
    let quality = compute_quality_assessment(&domain_agreement, &performance)?;
    let efficiency = compute_efficiency_metrics(execution_time_ms, x.len())?;

    Ok(Self {
      domain_agreement,
      performance,
      quality,
      efficiency,
    })
  }
}
