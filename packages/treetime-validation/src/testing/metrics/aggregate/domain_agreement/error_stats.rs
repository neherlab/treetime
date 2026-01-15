use approx::ulps_eq;
use itertools::izip;
use ndarray::Array1;
use ordered_float::OrderedFloat;

/// Absolute error statistics including bias detection
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct AbsoluteErrorStats {
  /// Mean absolute error: $\frac{1}{N}\sum_i \left|y_{\text{num}}(x_i) - y_{\text{ref}}(x_i)\right|$
  /// Average magnitude of deviations across all points
  pub mean: f64,
  /// Maximum absolute error: $\max_i \left|y_{\text{num}}(x_i) - y_{\text{ref}}(x_i)\right|$
  /// Largest deviation in absolute terms, indicates worst-case accuracy
  pub max: f64,
  /// Standard deviation of absolute errors
  /// Measures variability in error magnitudes across the domain
  pub std: f64,
  /// Signed error bias: $\frac{1}{N}\sum_i (y_{\text{num}}(x_i) - y_{\text{ref}}(x_i))$
  /// Detects systematic positive/negative bias in numerical methods
  pub bias: f64,
}

/// Relative error statistics with robust central tendency
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct RelativeErrorStats {
  /// Mean relative error: $\frac{1}{N}\sum_i \frac{y_{\text{num}}(x_i) - y_{\text{ref}}(x_i)}{y_{\text{ref}}(x_i)}$
  /// Average signed relative deviation, can indicate systematic scaling bias
  pub mean: f64,
  /// Maximum relative error: $\max_i \frac{\left|y_{\text{num}}(x_i) - y_{\text{ref}}(x_i)\right|}{\left|y_{\text{ref}}(x_i)\right|}$
  /// Worst-case relative deviation, normalized by reference magnitude
  pub max: f64,
  /// Mean Absolute Percentage Error: $\frac{100}{N}\sum_i \frac{\left|y_{\text{num}}(x_i) - y_{\text{ref}}(x_i)\right|}{\left|y_{\text{ref}}(x_i)\right|}$
  /// Average relative error magnitude as percentage, commonly used metric
  pub mape: f64,
  /// Median relative error: more robust to outliers than mean
  /// Uses $\frac{\left|y_{\text{num}}(x_i) - y_{\text{ref}}(x_i)\right|}{\max(\left|y_{\text{ref}}(x_i)\right|, \epsilon)}$
  pub median: f64,
}

/// Compute absolute error statistics with bias detection
/// Calculates mean, maximum, standard deviation, and signed bias of absolute errors
pub fn compute_absolute_error_statistics(actual: &Array1<f64>, expected: &Array1<f64>) -> AbsoluteErrorStats {
  let abs_errors: Vec<f64> = (actual - expected).mapv(|x| x.abs()).to_vec();
  let signed_errors: Vec<f64> = (actual - expected).to_vec();

  let mean = abs_errors.iter().sum::<f64>() / abs_errors.len() as f64;
  let max = abs_errors.iter().map(|&x| OrderedFloat(x)).max().map_or(0.0, |x| x.0);
  let bias = signed_errors.iter().sum::<f64>() / signed_errors.len() as f64;

  let variance = abs_errors.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / abs_errors.len() as f64;
  let std = variance.sqrt();

  AbsoluteErrorStats { mean, max, std, bias }
}

/// Compute relative error statistics
/// Returns mean, max, MAPE, and median relative errors
pub fn compute_relative_error_statistics(actual: &Array1<f64>, expected: &Array1<f64>) -> RelativeErrorStats {
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
