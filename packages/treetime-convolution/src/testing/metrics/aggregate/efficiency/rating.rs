use serde::{Deserialize, Serialize};

/// Computational efficiency classification
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum EfficiencyRating {
  /// Very fast execution
  VeryFast,
  /// Fast execution
  Fast,
  /// Moderate execution time
  Moderate,
  /// Slow execution
  Slow,
  /// Very slow execution
  VerySlow,
}
