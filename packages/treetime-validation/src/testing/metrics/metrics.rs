#[cfg(test)]
mod __tests__;

use crate::testing::metrics::aggregate::aggregate::AggregateMetrics;
use crate::testing::metrics::config::MetricsConfig;
use crate::testing::metrics::distribution::distribution::DistributionMetrics;
use crate::testing::metrics::pointwise::pointwise::PointwiseMetrics;
use crate::testing::metrics::spatial::spatial::SpatialMetrics;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use treetime_utils::make_error;

#[cfg(test)]
use approx::assert_ulps_eq;

/// Comprehensive convolution algorithm metrics
///
/// This structure cleanly separates different types of analysis:
/// - Aggregate metrics provide overall assessment and quality scores
/// - Pointwise metrics give detailed per-point error analysis
/// - Spatial metrics analyze regional patterns and windowed behavior
/// - Distribution metrics examine statistical properties and histograms
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationMetrics {
  /// Domain-wide aggregate accuracy and quality metrics
  pub aggregate: AggregateMetrics,

  /// Point-by-point error analysis (1:1 with evaluation grid)
  pub pointwise: PointwiseMetrics,

  /// Spatial and regional analysis metrics
  pub spatial: SpatialMetrics,

  /// Statistical distribution analysis
  pub distribution: DistributionMetrics,
}

impl ValidationMetrics {
  /// Creates comprehensive metrics from evaluation data
  pub fn new(
    x: &Array1<f64>,
    actual: &Array1<f64>,
    expected: &Array1<f64>,
    execution_time_ms: f64,
  ) -> eyre::Result<Self> {
    Self::new_with_config(x, actual, expected, execution_time_ms, &MetricsConfig::default())
  }

  /// Creates comprehensive metrics with custom configuration
  pub fn new_with_config(
    x: &Array1<f64>,
    actual: &Array1<f64>,
    expected: &Array1<f64>,
    execution_time_ms: f64,
    config: &MetricsConfig,
  ) -> eyre::Result<Self> {
    let dx = compute_grid_spacing(x)?;

    // Compute each metric type independently
    let aggregate = AggregateMetrics::new(x, actual, expected, execution_time_ms)?;
    let pointwise = PointwiseMetrics::new_with_config(x, actual, expected, &config.pointwise)?;
    let spatial = SpatialMetrics::new_with_config(x, actual, expected, dx, &config.spatial)?;
    let distribution = DistributionMetrics::new_with_config(&pointwise.errors, &config.distribution)?;

    Ok(Self {
      aggregate,
      pointwise,
      spatial,
      distribution,
    })
  }
}

/// Helper function to compute grid spacing
fn compute_grid_spacing(x: &Array1<f64>) -> eyre::Result<f64> {
  if x.len() < 2 {
    return make_error!("Cannot compute grid spacing with fewer than 2 points");
  }
  let dx = x[1] - x[0];
  if !dx.is_finite() || dx <= 0.0 {
    return make_error!("Invalid grid spacing: {dx}");
  }
  Ok(dx)
}
