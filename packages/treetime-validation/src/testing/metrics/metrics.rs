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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationMetrics {
  pub aggregate: AggregateMetrics,

  pub pointwise: PointwiseMetrics,

  pub spatial: SpatialMetrics,

  pub distribution: DistributionMetrics,
}

impl ValidationMetrics {
  pub fn new(
    x: &Array1<f64>,
    actual: &Array1<f64>,
    expected: &Array1<f64>,
    execution_time_ms: f64,
  ) -> eyre::Result<Self> {
    Self::new_with_config(x, actual, expected, execution_time_ms, &MetricsConfig::default())
  }

  pub fn new_with_config(
    x: &Array1<f64>,
    actual: &Array1<f64>,
    expected: &Array1<f64>,
    execution_time_ms: f64,
    config: &MetricsConfig,
  ) -> eyre::Result<Self> {
    let dx = compute_grid_spacing(x)?;

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
