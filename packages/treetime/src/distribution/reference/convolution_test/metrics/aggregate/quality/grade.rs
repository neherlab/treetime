use serde::{Deserialize, Serialize};

/// Quality grade classification
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum QualityGrade {
  /// Excellent quality (90-100)
  A,
  /// Good quality (80-89)
  B,
  /// Acceptable quality (70-79)
  C,
  /// Poor quality (60-69)
  D,
  /// Unacceptable quality (<60)
  F,
}
