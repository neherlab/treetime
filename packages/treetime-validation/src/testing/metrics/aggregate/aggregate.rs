use crate::testing::metrics::aggregate::domain_agreement::domain_agreement::DomainAgreementMetrics;
use crate::testing::metrics::aggregate::performance::performance::{PerformanceMetrics, compute_performance_metrics};
use crate::testing::metrics::config::ToleranceThresholds;
use eyre::Result;
use ndarray::Array1;
use serde::{Deserialize, Serialize};

/// Aggregate metrics combining domain agreement with performance assessment
///
/// These metrics provide domain-wide scalar summaries that give an overall
/// assessment of algorithm accuracy across all evaluation points.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AggregateMetrics {
  /// Domain agreement metrics
  pub domain_agreement: DomainAgreementMetrics,
  /// Performance assessment
  pub performance: PerformanceMetrics,
  /// Execution time in milliseconds
  pub execution_time_ms: f64,
}

impl AggregateMetrics {
  pub fn new(x: &Array1<f64>, actual: &Array1<f64>, expected: &Array1<f64>, execution_time_ms: f64) -> Result<Self> {
    let domain_agreement = DomainAgreementMetrics::new(x, actual, expected)?;
    let performance = compute_performance_metrics(actual, expected)?;

    Ok(Self {
      domain_agreement,
      performance,
      execution_time_ms,
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

    Ok(Self {
      domain_agreement,
      performance,
      execution_time_ms,
    })
  }
}
