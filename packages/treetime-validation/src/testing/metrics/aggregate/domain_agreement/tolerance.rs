use crate::testing::metrics::config::ToleranceThresholds;
use itertools::izip;
use ndarray::Array1;
use ordered_float::OrderedFloat;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ToleranceCounts {
  pub within_abs_tolerances: [usize; 3],
  pub within_rel_tolerances: [usize; 3],
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct MaxErrorLocation {
  pub idx: usize,
  pub x_value: f64,
}

pub(crate) fn compute_tolerance_counts(
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

pub(crate) fn find_max_error_location(
  x: &Array1<f64>,
  actual: &Array1<f64>,
  expected: &Array1<f64>,
) -> MaxErrorLocation {
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
