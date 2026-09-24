use approx::ulps_eq;
use itertools::izip;
use ndarray::Array1;
use ordered_float::OrderedFloat;

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub(crate) fn compute_absolute_error_statistics(actual: &Array1<f64>, expected: &Array1<f64>) -> AbsoluteErrorStats {
  let abs_errors: Vec<f64> = (actual - expected).mapv(|x| x.abs()).to_vec();
  let signed_errors: Vec<f64> = (actual - expected).to_vec();

  let mean = abs_errors.iter().sum::<f64>() / abs_errors.len() as f64;
  let max = abs_errors.iter().map(|&x| OrderedFloat(x)).max().map_or(0.0, |x| x.0);
  let bias = signed_errors.iter().sum::<f64>() / signed_errors.len() as f64;

  let variance = abs_errors.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / abs_errors.len() as f64;
  let std = variance.sqrt();

  AbsoluteErrorStats { mean, max, std, bias }
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct AbsoluteErrorStats {
  pub(crate) mean: f64,
  pub(crate) max: f64,
  pub(crate) std: f64,
  pub(crate) bias: f64,
}

#[allow(
  clippy::as_conversions,
  clippy::integer_division,
  reason = "count/index numeric cast is exact for the domain range; integer division is the intended floor division"
)]
pub(crate) fn compute_relative_error_statistics(actual: &Array1<f64>, expected: &Array1<f64>) -> RelativeErrorStats {
  let mut rel_errors = Vec::new();
  let mut abs_rel_errors = Vec::new();

  for (&a, &e) in izip!(actual.iter(), expected.iter()) {
    if !ulps_eq!(e, 0.0, max_ulps = 3) {
      let rel_error = (a - e) / e;
      rel_errors.push(rel_error);
      abs_rel_errors.push(rel_error.abs());
    }
  }

  if rel_errors.is_empty() {
    return RelativeErrorStats {
      mean: 0.0,
      max: 0.0,
      mape: 0.0,
      median: 0.0,
    };
  }

  let mean = rel_errors.iter().sum::<f64>() / rel_errors.len() as f64;
  let max = abs_rel_errors
    .iter()
    .map(|&x| OrderedFloat(x))
    .max()
    .map_or(0.0, |x| x.0);
  let mape = (abs_rel_errors.iter().sum::<f64>() / abs_rel_errors.len() as f64) * 100.0;

  let mut sorted_abs_errors = abs_rel_errors.clone();
  sorted_abs_errors.sort_by_key(|&x| OrderedFloat(x));
  let median = if sorted_abs_errors.len() % 2 == 0 {
    let mid = sorted_abs_errors.len() / 2;
    f64::midpoint(sorted_abs_errors[mid - 1], sorted_abs_errors[mid])
  } else {
    sorted_abs_errors[sorted_abs_errors.len() / 2]
  };

  RelativeErrorStats {
    mean,
    max,
    mape,
    median,
  }
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct RelativeErrorStats {
  pub(crate) mean: f64,
  pub(crate) max: f64,
  pub(crate) mape: f64,
  pub(crate) median: f64,
}
