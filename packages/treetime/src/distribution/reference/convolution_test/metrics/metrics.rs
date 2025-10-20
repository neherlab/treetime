use crate::distribution::reference::convolution_test::metrics::aggregate::aggregate::AggregateMetrics;
use crate::distribution::reference::convolution_test::metrics::config::MetricsConfig;
use crate::distribution::reference::convolution_test::metrics::distribution::distribution::DistributionMetrics;
use crate::distribution::reference::convolution_test::metrics::pointwise::pointwise::PointwiseMetrics;
use crate::distribution::reference::convolution_test::metrics::spatial::spatial::SpatialMetrics;
use crate::make_error;
use ndarray::Array1;
use serde::{Deserialize, Serialize};

/// Comprehensive convolution algorithm metrics
///
/// This structure cleanly separates different types of analysis:
/// - Aggregate metrics provide overall assessment and quality scores
/// - Pointwise metrics give detailed per-point error analysis
/// - Spatial metrics analyze regional patterns and windowed behavior
/// - Distribution metrics examine statistical properties and histograms
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConvolutionMetrics {
  /// Domain-wide aggregate accuracy and quality metrics
  pub aggregate: AggregateMetrics,

  /// Point-by-point error analysis (1:1 with evaluation grid)
  pub pointwise: PointwiseMetrics,

  /// Spatial and regional analysis metrics
  pub spatial: SpatialMetrics,

  /// Statistical distribution analysis
  pub distribution: DistributionMetrics,
}

impl ConvolutionMetrics {
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

#[cfg(test)]
mod tests {
  use super::*;
  use ndarray::array;

  #[test]
  fn test_comprehensive_metrics_creation() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let expected = array![1.0, 2.0, 3.0, 2.0, 1.0];
    let actual = &expected + 0.1; // Small systematic error

    let metrics = ConvolutionMetrics::new(&x, &actual, &expected, 100.0).unwrap();

    // Check that all metric types are computed
    assert_eq!(metrics.pointwise.total_points, 5);
    assert_eq!(metrics.spatial.total_points, 5);
    assert!(metrics.aggregate.domain_agreement.quality_metrics.r_squared > 0.9);
    assert!(metrics.distribution.statistics.abs_error_stats.mean > 0.0);
  }

  #[test]
  fn test_metrics_with_custom_config() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y = array![1.0, 2.0, 3.0, 2.0, 1.0];

    let mut config = MetricsConfig::default();
    config.distribution.histogram_bins = 20;
    config.spatial.window_half_width = 2;

    let metrics = ConvolutionMetrics::new_with_config(&x, &y, &y, 50.0, &config).unwrap();

    assert_eq!(metrics.distribution.histograms.abs_error_histogram.bin_counts.len(), 20);
    // TODO: Test quality grade when PerformanceMetrics, QualityAssessment, EfficiencyMetrics are added
  }
}
