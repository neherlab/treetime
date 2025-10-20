pub mod aggregate;
pub mod domain_agreement;
pub mod efficiency;
pub mod performance;
pub mod quality;

pub use aggregate::AggregateMetrics;
pub use domain_agreement::{AgreementAssessment, DomainAgreementMetrics};
pub use efficiency::{
  EfficiencyMetrics, EfficiencyRating, MemoryEfficiency, MemoryRating, ScalabilityAssessment,
  compute_efficiency_metrics,
};
pub use performance::{PerformanceMetrics, compute_performance_metrics};
pub use quality::{
  QualityAssessment, QualityComponents, QualityGrade, QualityVerdict, compute_accuracy_score, compute_precision_score,
  compute_quality_assessment, compute_robustness_score, compute_stability_score, generate_quality_verdict,
};
