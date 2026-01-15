use approx::ulps_eq;
use ndarray::Array1;
use ndarray_stats::QuantileExt;
use treetime_utils::make_error;

/// Peak-related accuracy metrics for distribution analysis
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PeakMetrics {
  /// Peak value error: $\frac{\left|\max_i y_{\text{num}}(x_i) - \max_i y_{\text{ref}}(x_i)\right|}{\max_i y_{\text{ref}}(x_i)}$
  /// Relative error in maximum amplitude, critical for probability distributions
  pub value_error: f64,
  /// Peak location error: $\left|x_{\operatorname{argmax}(y_{\text{num}})} - x_{\operatorname{argmax}(y_{\text{ref}})}\right|$
  /// Difference in location of the peak, important for temporal accuracy
  pub location_error: f64,
}

/// Compute peak-related accuracy metrics for distribution analysis
/// Analyzes both peak amplitude and location accuracy
pub fn compute_peak_metrics(
  x: &Array1<f64>,
  actual: &Array1<f64>,
  expected: &Array1<f64>,
) -> eyre::Result<PeakMetrics> {
  let actual_peak = *actual.max().unwrap_or(&0.0);
  let expected_peak = *expected.max().unwrap_or(&0.0);

  let actual_peak_idx = actual.argmax().unwrap_or(0);
  let expected_peak_idx = expected.argmax().unwrap_or(0);

  if ulps_eq!(expected_peak, 0.0, max_ulps = 3) {
    return make_error!("Expected peak value too close to zero: {:.2e}", expected_peak);
  }
  let value_error = (actual_peak - expected_peak).abs() / expected_peak.abs();

  let location_error = (x[actual_peak_idx] - x[expected_peak_idx]).abs();

  Ok(PeakMetrics {
    value_error,
    location_error,
  })
}
