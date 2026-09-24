use itertools::izip;
use ndarray::Array1;
use ordered_float::OrderedFloat;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct QualityMetrics {
  pub(crate) rmse: f64,
  pub(crate) r_squared: f64,
  pub(crate) correlation: f64,
  pub(crate) mass_error: f64,
  pub(crate) rel_l2_error: f64,
  pub(crate) rel_l1_error: f64,
  pub(crate) rel_linf_error: f64,
  pub(crate) max_log_error: f64,
  pub(crate) symmetry_error: f64,
  pub(crate) quantile_95_error: f64,
}

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub(crate) fn compute_rmse(actual: &Array1<f64>, expected: &Array1<f64>) -> f64 {
  let squared_errors: f64 = (actual - expected).mapv(|x| x * x).sum();
  (squared_errors / actual.len() as f64).sqrt()
}

pub(crate) fn compute_r_squared(actual: &Array1<f64>, expected: &Array1<f64>) -> f64 {
  let mean_expected = expected.mean().unwrap_or(0.0);
  let ss_res: f64 = (actual - expected).mapv(|x| x * x).sum();
  let ss_tot: f64 = expected.mapv(|x| (x - mean_expected).powi(2)).sum();

  if ss_tot.abs() < 1e-15 {
    1.0
  } else {
    (1.0 - ss_res / ss_tot).max(0.0)
  }
}

pub(crate) fn compute_correlation(actual: &Array1<f64>, expected: &Array1<f64>) -> f64 {
  let mean_actual = actual.mean().unwrap_or(0.0);
  let mean_expected = expected.mean().unwrap_or(0.0);

  let cov: f64 = izip!(actual.iter(), expected.iter())
    .map(|(&a, &e)| (a - mean_actual) * (e - mean_expected))
    .sum();

  let var_actual: f64 = actual.iter().map(|&a| (a - mean_actual).powi(2)).sum();
  let var_expected: f64 = expected.iter().map(|&e| (e - mean_expected).powi(2)).sum();

  let denominator = (var_actual * var_expected).sqrt();
  if denominator < 1e-15 { 1.0 } else { cov / denominator }
}

pub(crate) fn compute_mass_error(x: &Array1<f64>, actual: &Array1<f64>, expected: &Array1<f64>) -> f64 {
  let compute_integral = |values: &Array1<f64>| -> f64 {
    let mut integral = 0.0;
    for i in 0..values.len() - 1 {
      let dx = x[i + 1] - x[i];
      integral += f64::midpoint(values[i], values[i + 1]) * dx;
    }
    integral
  };

  let actual_integral = compute_integral(actual);
  let expected_integral = compute_integral(expected);

  (actual_integral - expected_integral).abs()
}

pub(crate) fn compute_relative_l2_norm_error(actual: &Array1<f64>, expected: &Array1<f64>) -> f64 {
  let error_norm = (actual - expected).mapv(|x| x * x).sum().sqrt();
  let expected_norm = expected.mapv(|x| x * x).sum().sqrt();

  if expected_norm < 1e-15 {
    0.0
  } else {
    error_norm / expected_norm
  }
}

pub(crate) fn compute_relative_l1_norm_error(actual: &Array1<f64>, expected: &Array1<f64>) -> f64 {
  let error_norm: f64 = (actual - expected).mapv(|x| x.abs()).sum();
  let expected_norm: f64 = expected.mapv(|x| x.abs()).sum();

  if expected_norm < 1e-15 {
    0.0
  } else {
    error_norm / expected_norm
  }
}

pub(crate) fn compute_relative_linf_norm_error(actual: &Array1<f64>, expected: &Array1<f64>) -> f64 {
  let error_norm = (actual - expected)
    .iter()
    .map(|&x| OrderedFloat(x.abs()))
    .max()
    .map_or(0.0, |x| x.0);

  let expected_norm = expected
    .iter()
    .map(|&x| OrderedFloat(x.abs()))
    .max()
    .map_or(0.0, |x| x.0);

  if expected_norm < 1e-15 {
    0.0
  } else {
    error_norm / expected_norm
  }
}

pub(crate) fn compute_max_log_error(actual: &Array1<f64>, expected: &Array1<f64>, threshold: f64) -> f64 {
  let mut max_log_error = 0.0_f64;
  for (&a, &e) in izip!(actual.iter(), expected.iter()) {
    if e > threshold && a > threshold {
      let log_error = (a.ln() - e.ln()).abs();
      max_log_error = max_log_error.max(log_error);
    }
  }
  max_log_error
}

pub(crate) fn compute_symmetry_error(x: &Array1<f64>, actual: &Array1<f64>) -> f64 {
  let mut max_symmetry_error = 0.0_f64;
  let mut found_symmetric_pairs = false;

  for (i, &xi) in x.iter().enumerate() {
    let neg_xi = -xi;
    for (j, &xj) in x.iter().enumerate() {
      if (xj - neg_xi).abs() < 1e-10 {
        let symmetry_error = (actual[i] - actual[j]).abs();
        max_symmetry_error = max_symmetry_error.max(symmetry_error);
        found_symmetric_pairs = true;
      }
    }
  }

  if !found_symmetric_pairs {
    0.0
  } else {
    max_symmetry_error
  }
}

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub(crate) fn compute_quantile_error(actual: &Array1<f64>, expected: &Array1<f64>, quantile: f64) -> f64 {
  let mut abs_errors: Vec<f64> = (actual - expected).mapv(|x| x.abs()).to_vec();
  abs_errors.sort_by_key(|&x| OrderedFloat(x));

  let index = ((quantile * abs_errors.len() as f64).ceil() as usize).saturating_sub(1);
  abs_errors.get(index).copied().unwrap_or(0.0)
}
