use ndarray::Array1;
use serde::{Deserialize, Serialize};

/// Performance metrics for numerical algorithm comparison
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerformanceMetrics {
  /// Signal-to-noise ratio in dB
  pub signal_to_noise_ratio: f64,
  /// Normalized root mean square error
  pub normalized_rmse: f64,
  /// Coefficient of determination (R²)
  pub coefficient_of_determination: f64,
  /// Nash-Sutcliffe efficiency coefficient
  pub nash_sutcliffe_efficiency: f64,
}

pub fn compute_performance_metrics(actual: &Array1<f64>, expected: &Array1<f64>) -> eyre::Result<PerformanceMetrics> {
  let errors = actual - expected;
  let squared_errors = errors.mapv(|x| x * x);

  let mse = squared_errors.mean().unwrap_or(0.0);
  let rmse = mse.sqrt();

  let signal_power = expected.mapv(|x| x * x).mean().unwrap_or(0.0);
  let noise_power = mse;

  let signal_to_noise_ratio = if noise_power > 0.0 {
    10.0 * (signal_power / noise_power).log10()
  } else {
    f64::INFINITY
  };

  let expected_range =
    expected.iter().copied().fold(0.0_f64, f64::max) - expected.iter().copied().fold(f64::INFINITY, f64::min);
  let normalized_rmse = if expected_range > 0.0 {
    rmse / expected_range
  } else {
    0.0
  };

  let expected_mean = expected.mean().unwrap_or(0.0);
  let total_sum_squares = expected.mapv(|x| (x - expected_mean).powi(2)).sum();
  let residual_sum_squares = squared_errors.sum();

  let coefficient_of_determination = if total_sum_squares > 0.0 {
    1.0 - (residual_sum_squares / total_sum_squares)
  } else {
    0.0
  };

  let nash_sutcliffe_efficiency = coefficient_of_determination;

  Ok(PerformanceMetrics {
    signal_to_noise_ratio,
    normalized_rmse,
    coefficient_of_determination,
    nash_sutcliffe_efficiency,
  })
}
