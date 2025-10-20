use serde::{Deserialize, Serialize};

/// Algorithm scalability characteristics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScalabilityAssessment {
  /// Expected scaling behavior with problem size
  pub scaling_behavior: String,
  /// Predicted performance at 10x problem size
  pub performance_10x: f64,
  /// Predicted performance at 100x problem size
  pub performance_100x: f64,
}
