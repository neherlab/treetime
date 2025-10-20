use itertools::izip;
use ndarray::Array1;
use ordered_float::OrderedFloat;

/// Quality metrics for agreement assessment including conservation properties
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct QualityMetrics {
  /// Root Mean Square Error: $\sqrt{\frac{1}{N}\sum_i (y_{\text{num}}(x_i) - y_{\text{ref}}(x_i))^2}$
  /// Global measure of deviation magnitude, sensitive to large errors
  pub rmse: f64,
  /// Coefficient of determination: $R^2 = 1 - \frac{SS_{\text{res}}}{SS_{\text{tot}}}$
  /// Fraction of variance explained
  pub r_squared: f64,
  /// Pearson correlation coefficient: measures linear relationship strength
  pub correlation: f64,
  /// Integral (mass) conservation error: $\left|\sum_i y_{\text{num}}(x_i)\Delta x - \sum_i y_{\text{ref}}(x_i)\Delta x\right|$
  /// Critical for probability distributions that must integrate to 1
  pub mass_error: f64,
  /// Relative L2 norm error: $\frac{\left\|y_{\text{num}} - y_{\text{ref}}\right\|_2}{\left\|y_{\text{ref}}\right\|_2}$
  /// Scale-invariant measure of total deviation
  pub rel_l2_error: f64,
  /// Relative L1 norm error: $\frac{\sum_i \left|y_{\text{num}}(x_i) - y_{\text{ref}}(x_i)\right|}{\sum_i \left|y_{\text{ref}}(x_i)\right|}$
  /// Robust to outliers, global deviation measure
  pub rel_l1_error: f64,
  /// Relative L∞ norm error: $\frac{\max_i \left|y_{\text{num}}(x_i) - y_{\text{ref}}(x_i)\right|}{\max_i \left|y_{\text{ref}}(x_i)\right|}$
  /// Supremum error normalized by global peak
  pub rel_linf_error: f64,
  /// Maximum log error: $\max_{i: y_{\text{ref}}(x_i) > \tau} \left|\log y_{\text{num}}(x_i) - \log y_{\text{ref}}(x_i)\right|$
  /// Sensitive to tail accuracy for small values
  pub max_log_error: f64,
  /// Symmetry error: $\max_i \left|y_{\text{num}}(x_i) - y_{\text{num}}(-x_i)\right|$
  /// For symmetric kernels, measures deviation from symmetry
  pub symmetry_error: f64,
  /// 95th percentile error threshold
  /// Value below which 95% of absolute errors fall
  pub quantile_95_error: f64,
}

/// Compute root-mean-square error
pub fn compute_rmse(actual: &Array1<f64>, expected: &Array1<f64>) -> f64 {
  let squared_errors: f64 = (actual - expected).mapv(|x| x * x).sum();
  (squared_errors / actual.len() as f64).sqrt()
}

/// Compute R-squared
pub fn compute_r_squared(actual: &Array1<f64>, expected: &Array1<f64>) -> f64 {
  let mean_expected = expected.mean().unwrap_or(0.0);
  let ss_res: f64 = (actual - expected).mapv(|x| x * x).sum();
  let ss_tot: f64 = expected.mapv(|x| (x - mean_expected).powi(2)).sum();

  if ss_tot.abs() < 1e-15 {
    1.0
  } else {
    (1.0 - ss_res / ss_tot).max(0.0)
  }
}

/// Compute Pearson correlation coefficient
pub fn compute_correlation(actual: &Array1<f64>, expected: &Array1<f64>) -> f64 {
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

/// Compute integral mass error using trapezoidal integration
pub fn compute_mass_error(x: &Array1<f64>, actual: &Array1<f64>, expected: &Array1<f64>) -> f64 {
  let compute_integral = |values: &Array1<f64>| -> f64 {
    let mut integral = 0.0;
    for i in 0..values.len() - 1 {
      let dx = x[i + 1] - x[i];
      integral += 0.5 * (values[i] + values[i + 1]) * dx;
    }
    integral
  };

  let actual_integral = compute_integral(actual);
  let expected_integral = compute_integral(expected);

  (actual_integral - expected_integral).abs()
}

/// Compute relative L2 norm error
pub fn compute_relative_l2_norm_error(actual: &Array1<f64>, expected: &Array1<f64>) -> f64 {
  let error_norm = (actual - expected).mapv(|x| x * x).sum().sqrt();
  let expected_norm = expected.mapv(|x| x * x).sum().sqrt();

  if expected_norm < 1e-15 {
    0.0
  } else {
    error_norm / expected_norm
  }
}

/// Compute relative L1 norm error
pub fn compute_relative_l1_norm_error(actual: &Array1<f64>, expected: &Array1<f64>) -> f64 {
  let error_norm: f64 = (actual - expected).mapv(|x| x.abs()).sum();
  let expected_norm: f64 = expected.mapv(|x| x.abs()).sum();

  if expected_norm < 1e-15 {
    0.0
  } else {
    error_norm / expected_norm
  }
}

/// Compute relative L-infinity norm error
pub fn compute_relative_linf_norm_error(actual: &Array1<f64>, expected: &Array1<f64>) -> f64 {
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

/// Compute maximum logarithmic error
pub fn compute_max_log_error(actual: &Array1<f64>, expected: &Array1<f64>, threshold: f64) -> f64 {
  let mut max_log_error = 0.0_f64;
  for (&a, &e) in izip!(actual.iter(), expected.iter()) {
    if e > threshold && a > threshold {
      let log_error = (a.ln() - e.ln()).abs();
      max_log_error = max_log_error.max(log_error);
    }
  }
  max_log_error
}

/// Compute symmetry error for symmetric distributions
pub fn compute_symmetry_error(x: &Array1<f64>, actual: &Array1<f64>) -> f64 {
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

/// Compute quantile error
pub fn compute_quantile_error(actual: &Array1<f64>, expected: &Array1<f64>, quantile: f64) -> f64 {
  let mut abs_errors: Vec<f64> = (actual - expected).mapv(|x| x.abs()).to_vec();
  abs_errors.sort_by_key(|&x| OrderedFloat(x));

  let index = ((quantile * abs_errors.len() as f64).ceil() as usize).saturating_sub(1);
  abs_errors.get(index).copied().unwrap_or(0.0)
}
