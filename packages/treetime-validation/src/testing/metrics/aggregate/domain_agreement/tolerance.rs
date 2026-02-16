use crate::testing::metrics::config::ToleranceThresholds;
use itertools::izip;
use ndarray::Array1;
use ordered_float::OrderedFloat;

/// Tolerance counts for different threshold levels
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ToleranceCounts {
  /// Count of points within absolute tolerance thresholds [strict, moderate, loose]
  pub within_abs_tolerances: [usize; 3],
  /// Count of points within relative tolerance thresholds [strict, moderate, loose]
  pub within_rel_tolerances: [usize; 3],
}

/// Maximum error location information
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct MaxErrorLocation {
  /// Array index where maximum absolute error occurs
  pub idx: usize,
  /// X-coordinate value where maximum absolute error occurs
  pub x_value: f64,
}

/// Count points within specified tolerance thresholds
pub fn compute_tolerance_counts(
  actual: &Array1<f64>,
  expected: &Array1<f64>,
  thresholds: &ToleranceThresholds,
) -> ToleranceCounts {
  let mut within_abs_tolerances = [0_usize; 3];
  let mut within_rel_tolerances = [0_usize; 3];

  for (&a, &e) in izip!(actual.iter(), expected.iter()) {
    let abs_error = (a - e).abs();
    for (i, &threshold) in thresholds.abs_tolerances.iter().enumerate() {
      if abs_error < threshold {
        within_abs_tolerances[i] += 1;
      }
    }

    let rel_error = if e.abs() > 1e-15 { abs_error / e.abs() } else { 0.0 };
    for (i, &threshold) in thresholds.rel_tolerances.iter().enumerate() {
      if rel_error < threshold {
        within_rel_tolerances[i] += 1;
      }
    }
  }

  ToleranceCounts {
    within_abs_tolerances,
    within_rel_tolerances,
  }
}

/// Find location of maximum absolute error
pub fn find_max_error_location(x: &Array1<f64>, actual: &Array1<f64>, expected: &Array1<f64>) -> MaxErrorLocation {
  let abs_errors: Vec<f64> = (actual - expected).mapv(|x| x.abs()).to_vec();
  let max_idx = abs_errors
    .iter()
    .enumerate()
    .max_by_key(|&(_, v)| OrderedFloat(*v))
    .map_or(0, |(idx, _)| idx);

  MaxErrorLocation {
    idx: max_idx,
    x_value: x[max_idx],
  }
}
