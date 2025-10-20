use crate::distribution::reference::convolution_test::metrics::config::SpatialConfig;
use ndarray::Array1;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WindowedMetrics {
  pub sliding_rms: Array1<f64>,
  pub sliding_max: Array1<f64>,
  pub summary: WindowedSummary,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WindowedSummary {
  pub sliding_rms_max: f64,
  pub sliding_max_max: f64,
}

pub(super) fn compute_windowed_metrics(
  actual: &Array1<f64>,
  expected: &Array1<f64>,
  config: &SpatialConfig,
) -> eyre::Result<WindowedMetrics> {
  let errors = actual - expected;
  let abs_errors = errors.mapv(|x| x.abs());

  let sliding_rms = compute_sliding_window_rms(&errors, config.window_half_width);
  let sliding_max = compute_sliding_window_max(&abs_errors, config.window_half_width);

  let sliding_rms_max = sliding_rms.iter().copied().fold(0.0_f64, f64::max);
  let sliding_max_max = sliding_max.iter().copied().fold(0.0_f64, f64::max);

  let summary = WindowedSummary {
    sliding_rms_max,
    sliding_max_max,
  };

  Ok(WindowedMetrics {
    sliding_rms,
    sliding_max,
    summary,
  })
}

fn compute_sliding_window_rms(errors: &Array1<f64>, half_width: usize) -> Array1<f64> {
  let n = errors.len();
  let mut rms = Array1::zeros(n);

  for i in 0..n {
    let start = i.saturating_sub(half_width);
    let end = (i + half_width + 1).min(n);
    let window = &errors.slice(ndarray::s![start..end]);
    let sum_sq: f64 = window.iter().map(|&x| x * x).sum();
    let count = window.len() as f64;
    rms[i] = (sum_sq / count).sqrt();
  }

  rms
}

fn compute_sliding_window_max(abs_errors: &Array1<f64>, half_width: usize) -> Array1<f64> {
  let n = abs_errors.len();
  let mut max_err = Array1::zeros(n);

  for i in 0..n {
    let start = i.saturating_sub(half_width);
    let end = (i + half_width + 1).min(n);
    let window = &abs_errors.slice(ndarray::s![start..end]);
    max_err[i] = window.iter().copied().fold(0.0_f64, f64::max);
  }

  max_err
}
