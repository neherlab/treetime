use crate::testing::metrics::aggregate::domain_agreement::domain_agreement::DomainAgreementMetrics;
use crate::testing::metrics::aggregate::performance::performance::PerformanceMetrics;
use crate::testing::metrics::aggregate::quality::components::{
  QualityComponents, compute_accuracy_score, compute_precision_score, compute_robustness_score, compute_stability_score,
};
use crate::testing::metrics::aggregate::quality::grade::QualityGrade;
use crate::testing::metrics::aggregate::quality::verdict::{QualityVerdict, generate_quality_verdict};
use eyre::Result;
use serde::{Deserialize, Serialize};

/// Quality assessment with overall score, grade, and detailed analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityAssessment {
  pub overall_score: f64,
  pub grade: QualityGrade,
  pub components: QualityComponents,
  pub verdict: QualityVerdict,
}

pub fn compute_quality_assessment(
  domain_agreement: &DomainAgreementMetrics,
  performance: &PerformanceMetrics,
) -> Result<QualityAssessment> {
  // Compute individual component scores
  let accuracy_score = compute_accuracy_score(domain_agreement);
  let precision_score = compute_precision_score(performance);
  let stability_score = compute_stability_score(domain_agreement);
  let robustness_score = compute_robustness_score(performance);

  // Overall score is weighted average
  let overall_score = 0.4 * accuracy_score + 0.3 * precision_score + 0.2 * stability_score + 0.1 * robustness_score;

  let components = QualityComponents {
    accuracy_score,
    precision_score,
    stability_score,
    robustness_score,
    overall_score,
  };

  let grade = match overall_score {
    s if s >= 90.0 => QualityGrade::A,
    s if s >= 80.0 => QualityGrade::B,
    s if s >= 70.0 => QualityGrade::C,
    s if s >= 60.0 => QualityGrade::D,
    _ => QualityGrade::F,
  };

  let verdict = generate_quality_verdict(domain_agreement, performance, &components);

  Ok(QualityAssessment {
    overall_score,
    grade,
    components,
    verdict,
  })
}
