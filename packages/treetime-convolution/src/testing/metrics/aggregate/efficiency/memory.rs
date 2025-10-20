use serde::{Deserialize, Serialize};

/// Memory usage efficiency assessment
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MemoryEfficiency {
  /// Estimated memory complexity (e.g., "O(n)", "O(n log n)")
  pub complexity_class: String,
  /// Efficiency rating
  pub rating: MemoryRating,
}

/// Memory efficiency classification
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MemoryRating {
  /// Very memory efficient
  Excellent,
  /// Good memory usage
  Good,
  /// Moderate memory usage
  Moderate,
  /// High memory usage
  Poor,
  /// Excessive memory usage
  Critical,
}
