use approx::ulps_eq;
use eyre::WrapErr;
use ndarray::Array1;
use ndarray_stats::QuantileExt;
use treetime_utils::array::ndarray::argmax_first;
use treetime_utils::make_error;

pub(crate) fn compute_peak_metrics(
  x: &Array1<f64>,
  actual: &Array1<f64>,
  expected: &Array1<f64>,
) -> eyre::Result<PeakMetrics> {
  let actual_peak = *actual.max().wrap_err("When finding the peak of the actual values")?;
  let expected_peak = *expected
    .max()
    .wrap_err("When finding the peak of the expected values")?;

  let actual_peak_idx = argmax_first(&actual.view()).unwrap_or(0);
  let expected_peak_idx = argmax_first(&expected.view()).unwrap_or(0);

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

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PeakMetrics {
  pub value_error: f64,
  pub location_error: f64,
}
