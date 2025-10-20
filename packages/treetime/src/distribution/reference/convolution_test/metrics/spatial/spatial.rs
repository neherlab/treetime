use crate::distribution::reference::convolution_test::metrics::config::SpatialConfig;
use ndarray::Array1;
use serde::{Deserialize, Serialize};

use super::cumulative::{CumulativeMetrics, compute_cumulative_metrics};
use super::regional::{RegionalMetrics, compute_regional_metrics};
use super::windowed::{WindowedMetrics, compute_windowed_metrics};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpatialMetrics {
  pub total_points: usize,
  pub dx: f64,
  pub regional: RegionalMetrics,
  pub windowed: WindowedMetrics,
  pub cumulative: CumulativeMetrics,
}

impl SpatialMetrics {
  pub fn new(x: &Array1<f64>, actual: &Array1<f64>, expected: &Array1<f64>, dx: f64) -> eyre::Result<Self> {
    Self::new_with_config(x, actual, expected, dx, &SpatialConfig::default())
  }

  pub fn new_with_config(
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
  use approx::assert_abs_diff_eq;
  use ndarray::array;

  #[test]
  fn test_perfect_agreement() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y = array![1.0, 2.0, 3.0, 2.0, 1.0];
    let result = SpatialMetrics::new(&x, &y, &y, 1.0).unwrap();

    assert_eq!(result.total_points, 5);
    assert_abs_diff_eq!(result.cumulative.summary.final_value, 0.0, epsilon = 1e-12);
    assert_abs_diff_eq!(result.cumulative.summary.max_abs, 0.0, epsilon = 1e-12);
  }

  #[test]
  fn test_cumulative_error_systematic_bias() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let expected = array![0.0, 0.0, 0.0, 0.0, 0.0];
    let actual = array![0.1, 0.1, 0.1, 0.1, 0.1];

    let result = SpatialMetrics::new(&x, &actual, &expected, 1.0).unwrap();

    assert!(result.cumulative.summary.final_value > 0.3);
  }
}
