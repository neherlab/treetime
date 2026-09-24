use crate::testing::metrics::config::SpatialConfig;
use crate::testing::metrics::spatial::cumulative::{CumulativeMetrics, compute_cumulative_metrics};
use crate::testing::metrics::spatial::regional::{RegionalMetrics, compute_regional_metrics};
use crate::testing::metrics::spatial::windowed::{WindowedMetrics, compute_windowed_metrics};
use ndarray::Array1;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpatialMetrics {
  pub(crate) total_points: usize,
  pub(crate) dx: f64,
  pub(crate) regional: RegionalMetrics,
  pub(crate) windowed: WindowedMetrics,
  pub(crate) cumulative: CumulativeMetrics,
}

impl SpatialMetrics {
  pub(crate) fn new_with_config(
    x: &Array1<f64>,
    actual: &Array1<f64>,
    expected: &Array1<f64>,
    dx: f64,
    config: &SpatialConfig,
  ) -> eyre::Result<Self> {
    let n = x.len();

    let regional = compute_regional_metrics(x, actual, expected, config)?;
    let windowed = compute_windowed_metrics(actual, expected, config)?;
    let cumulative = compute_cumulative_metrics(actual, expected, dx)?;

    Ok(Self {
      total_points: n,
      dx,
      regional,
      windowed,
      cumulative,
    })
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_ulps_eq;
  use ndarray::array;

  #[test]
  fn test_perfect_agreement() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y = array![1.0, 2.0, 3.0, 2.0, 1.0];
    let result = SpatialMetrics::new_with_config(&x, &y, &y, 1.0, &SpatialConfig::default()).unwrap();

    assert_eq!(5, result.total_points);
    assert_ulps_eq!(result.cumulative.summary.final_value, 0.0, max_ulps = 4);
    assert_ulps_eq!(result.cumulative.summary.max_abs, 0.0, max_ulps = 4);
  }

  #[test]
  fn test_cumulative_error_systematic_bias() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let expected = array![0.0, 0.0, 0.0, 0.0, 0.0];
    let actual = array![0.1, 0.1, 0.1, 0.1, 0.1];

    let result = SpatialMetrics::new_with_config(&x, &actual, &expected, 1.0, &SpatialConfig::default()).unwrap();

    assert_ulps_eq!(result.cumulative.summary.final_value, 0.5, max_ulps = 4);
    assert_ulps_eq!(result.cumulative.summary.max_abs, 0.5, max_ulps = 4);
  }
}
