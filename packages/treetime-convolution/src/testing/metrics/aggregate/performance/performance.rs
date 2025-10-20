use ndarray::Array1;
use serde::{Deserialize, Serialize};

/// Enhanced performance metrics beyond basic domain agreement
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerformanceMetrics {
  /// Signal-to-noise ratio
  pub signal_to_noise_ratio: f64,
  /// Peak signal-to-noise ratio (for image/signal processing)
  pub peak_snr: f64,
  /// Normalized root mean square error
  pub normalized_rmse: f64,
  /// Mean absolute percentage error
  pub mean_absolute_percentage_error: f64,
  /// Symmetric mean absolute percentage error
  pub symmetric_mape: f64,
  /// Coefficient of determination (R²) alternative calculation
  pub coefficient_of_determination: f64,
  /// Nash-Sutcliffe efficiency
  pub nash_sutcliffe_efficiency: f64,
  /// Index of agreement (Willmott's d)
  pub index_of_agreement: f64,
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

  let max_expected = expected.iter().copied().fold(0.0_f64, f64::max);
  let peak_snr = if mse > 0.0 && max_expected > 0.0 {
    20.0 * (max_expected / rmse).log10()
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

  let mut mape_sum = 0.0;
  let mut smape_sum = 0.0;
  let mut valid_count = 0;

  for i in 0..actual.len() {
    let exp_val = expected[i];
    let act_val = actual[i];
    let abs_error = (act_val - exp_val).abs();

    if exp_val.abs() > 1e-10 {
      mape_sum += abs_error / exp_val.abs();
      smape_sum += abs_error / f64::midpoint(exp_val.abs(), act_val.abs());
      valid_count += 1;
    }
  }

  let mean_absolute_percentage_error = if valid_count > 0 {
    (mape_sum / valid_count as f64) * 100.0
  } else {
    0.0
  };

  let symmetric_mape = if valid_count > 0 {
    (smape_sum / valid_count as f64) * 100.0
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

  let potential_error_sum = actual
    .iter()
    .zip(expected.iter())
    .map(|(&a, &e)| (a - expected_mean).abs() + (e - expected_mean).abs())
    .sum::<f64>();

  let index_of_agreement = if potential_error_sum > 0.0 {
    1.0 - (residual_sum_squares / potential_error_sum.powi(2))
  } else {
    1.0
  };

  Ok(PerformanceMetrics {
    signal_to_noise_ratio,
    peak_snr,
    normalized_rmse,
    mean_absolute_percentage_error,
    symmetric_mape,
    coefficient_of_determination,
    nash_sutcliffe_efficiency,
    index_of_agreement,
  })
}
