use eyre::Result;
use serde::{Deserialize, Serialize};

use crate::distribution::reference::convolution_test::metrics::aggregate::efficiency::memory::{
  MemoryEfficiency, MemoryRating,
};
use crate::distribution::reference::convolution_test::metrics::aggregate::efficiency::rating::EfficiencyRating;
use crate::distribution::reference::convolution_test::metrics::aggregate::efficiency::scalability::ScalabilityAssessment;
use crate::o;

/// Algorithm efficiency and performance characteristics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EfficiencyMetrics {
  /// Execution time in milliseconds
  pub execution_time_ms: f64,
  /// Time per evaluation point (ms/point)
  pub time_per_point_ms: f64,
  /// Computational efficiency rating
  pub efficiency_rating: EfficiencyRating,
  /// Memory efficiency estimate
  pub memory_efficiency: MemoryEfficiency,
  /// Scalability assessment
  pub scalability: ScalabilityAssessment,
}

pub fn compute_efficiency_metrics(execution_time_ms: f64, num_points: usize) -> Result<EfficiencyMetrics> {
  let time_per_point_ms = execution_time_ms / num_points as f64;

  let efficiency_rating = match time_per_point_ms {
    t if t < 0.001 => EfficiencyRating::VeryFast,
    t if t < 0.01 => EfficiencyRating::Fast,
    t if t < 0.1 => EfficiencyRating::Moderate,
    t if t < 1.0 => EfficiencyRating::Slow,
    _ => EfficiencyRating::VerySlow,
  };

  // Estimate memory complexity (simplified heuristic)
  let memory_efficiency = if num_points < 1000 {
    MemoryEfficiency {
      complexity_class: o!("O(n)"),
      rating: MemoryRating::Excellent,
    }
  } else if num_points < 10000 {
    MemoryEfficiency {
      complexity_class: o!("O(n log n)"),
      rating: MemoryRating::Good,
    }
  } else {
    MemoryEfficiency {
      complexity_class: o!("O(n²)"),
      rating: MemoryRating::Moderate,
    }
  };

  // Scalability prediction (simplified)
  let performance_10x = execution_time_ms * 10.0; // Linear scaling assumption
  let performance_100x = execution_time_ms * 100.0;

  let scalability = ScalabilityAssessment {
    scaling_behavior: o!("Linear"),
    performance_10x,
    performance_100x,
  };

  Ok(EfficiencyMetrics {
    execution_time_ms,
    time_per_point_ms,
    efficiency_rating,
    memory_efficiency,
    scalability,
  })
}
