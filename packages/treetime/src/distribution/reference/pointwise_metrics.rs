use crate::make_error;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::BTreeMap;

/// Comprehensive pointwise and spatial accuracy metrics for convolution validation
///
/// Unlike aggregate metrics which produce single scalars, pointwise metrics produce arrays
/// aligned with the evaluation grid, enabling visualization of error spatial distribution,
/// regional analysis, and identification of problematic domain regions.
///
/// All metric arrays have the same length as the input evaluation grid.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PointwiseMetrics {
  /// Total number of evaluation points
  pub total_points: usize,
  /// Grid spacing (assumed uniform)
  pub dx: f64,
  /// Pointwise error metrics (arrays aligned with evaluation grid)
  pub pointwise_errors: PointwiseErrors,
  /// Derivative and structural error metrics
  pub structural_errors: StructuralErrors,
  /// Regional and spatial metrics with summary statistics
  pub spatial_metrics: SpatialMetrics,
  /// Error distribution analysis
  pub error_distributions: ErrorDistributions,
}

impl PointwiseMetrics {
  /// Creates new pointwise metrics from evaluation grid and function values
  ///
  /// All input arrays must have the same length and represent values on a uniform grid.
  pub fn new(x: &Array1<f64>, actual: &Array1<f64>, expected: &Array1<f64>) -> eyre::Result<Self> {
    Self::new_with_config(x, actual, expected, &PointwiseConfig::default())
  }

  /// Creates new pointwise metrics with custom configuration
  pub fn new_with_config(
    x: &Array1<f64>,
    actual: &Array1<f64>,
    expected: &Array1<f64>,
    config: &PointwiseConfig,
  ) -> eyre::Result<Self> {
    validate_inputs(x, actual, expected)?;

    let n = x.len();
    let dx = compute_grid_spacing(x)?;

    let pointwise_errors = compute_pointwise_errors(actual, expected, config)?;
    let structural_errors = compute_structural_errors(x, actual, expected, dx, config)?;
    let spatial_metrics = compute_spatial_metrics(x, actual, expected, dx, config)?;
    let error_distributions = compute_error_distributions(x, actual, expected, &pointwise_errors, config)?;

    Ok(Self {
      total_points: n,
      dx,
      pointwise_errors,
      structural_errors,
      spatial_metrics,
      error_distributions,
    })
  }
}

/// Configuration parameters for pointwise metric computation
#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
pub struct PointwiseConfig {
  /// Floor value for relative error computation to prevent division by zero
  #[default = 1e-15]
  pub epsilon: f64,
  /// Threshold for log error computation (values below threshold are masked)
  #[default = 1e-10]
  pub log_threshold: f64,
  /// Tolerance for monotonicity violation detection
  #[default = 1e-12]
  pub monotonicity_eta: f64,
  /// Peak region radius as multiple of effective width
  #[default = 3.0]
  pub peak_region_radius: f64,
  /// Tail threshold as fraction of peak value
  #[default = 1e-6]
  pub tail_threshold: f64,
  /// Sliding window half-width for spatial metrics (number of points)
  #[default = 5]
  pub window_half_width: usize,
  /// Absolute tolerance thresholds for pass/fail mask [strict, moderate, loose]
  #[default(_code = "[1e-6, 1e-9, 1e-12]")]
  pub abs_tolerances: [f64; 3],
  /// Relative tolerance thresholds for pass/fail mask [strict, moderate, loose]
  #[default(_code = "[0.01, 0.001, 0.0001]")]
  pub rel_tolerances: [f64; 3],
  /// Number of bins for error histogram
  #[default = 50]
  pub histogram_bins: usize,
}

/// Pointwise error metrics producing arrays aligned with evaluation grid
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PointwiseErrors {
  /// Absolute error: $e_{\text{abs}}(x_i) = |y_{\text{num}}(x_i) - y_{\text{ref}}(x_i)|$
  pub absolute: Array1<f64>,
  /// Relative error with floor: $e_{\text{rel}}(x_i) = \frac{|y_{\text{num}}(x_i) - y_{\text{ref}}(x_i)|}{\max(|y_{\text{ref}}(x_i)|, \epsilon)}$
  pub relative: Array1<f64>,
  /// Signed error: $e_{\text{sgn}}(x_i) = y_{\text{num}}(x_i) - y_{\text{ref}}(x_i)$
  pub signed: Array1<f64>,
  /// Logarithmic error for tail analysis (masked where values below threshold)
  /// $e_{\log}(x_i) = \log y_{\text{num}}(x_i) - \log y_{\text{ref}}(x_i)$ where both values $> \tau$
  pub logarithmic: Array1<f64>,
  /// Summary statistics for quick reference
  pub summary: PointwiseErrorSummary,
}

/// Summary statistics for pointwise errors
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PointwiseErrorSummary {
  /// Mean absolute error
  pub abs_mean: f64,
  /// Maximum absolute error
  pub abs_max: f64,
  /// Standard deviation of absolute errors
  pub abs_std: f64,
  /// Mean relative error
  pub rel_mean: f64,
  /// Maximum relative error
  pub rel_max: f64,
  /// Median relative error
  pub rel_median: f64,
  /// Signed error bias
  pub signed_bias: f64,
  /// Number of valid logarithmic error points
  pub log_valid_count: usize,
  /// Maximum logarithmic error (among valid points)
  pub log_max: f64,
}

/// Derivative and structural error metrics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StructuralErrors {
  /// First derivative error: $e_{D1}(x_i) = |D y_{\text{num}}(x_i) - D y_{\text{ref}}(x_i)|$
  /// Uses central differences: $D y(x_i) = \frac{y(x_{i+1}) - y(x_{i-1})}{2\Delta x}$
  pub first_derivative: Array1<f64>,
  /// Second derivative error: $e_{D2}(x_i) = |D^2 y_{\text{num}}(x_i) - D^2 y_{\text{ref}}(x_i)|$
  /// Uses central differences: $D^2 y(x_i) = \frac{y(x_{i+1}) - 2y(x_i) + y(x_{i-1})}{\Delta x^2}$
  pub second_derivative: Array1<f64>,
  /// Symmetry residual for symmetric distributions
  /// $r_{\text{sym}}(x_i) = y_{\text{num}}(x_i) - y_{\text{num}}(-x_i)$
  pub symmetry_residual: Array1<f64>,
  /// Monotonicity violation flags (1.0 where violation detected, 0.0 otherwise)
  /// Flags points where $y_{\text{num}}(x_{i+1}) - y_{\text{num}}(x_i) < -\eta$
  pub monotonicity_violations: Array1<f64>,
  /// Summary statistics
  pub summary: StructuralErrorSummary,
}

/// Summary statistics for structural errors
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StructuralErrorSummary {
  /// Maximum first derivative error
  pub d1_max: f64,
  /// Mean first derivative error
  pub d1_mean: f64,
  /// Maximum second derivative error
  pub d2_max: f64,
  /// Mean second derivative error
  pub d2_mean: f64,
  /// Maximum symmetry residual
  pub symmetry_max: f64,
  /// Number of monotonicity violations
  pub monotonicity_violation_count: usize,
}

/// Regional and spatial metrics with sliding window analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpatialMetrics {
  /// Peak region relative errors (masked outside peak region)
  pub peak_region_errors: Array1<f64>,
  /// Tail region log errors (masked outside tail regions)
  pub tail_region_errors: Array1<f64>,
  /// Sliding window RMS error profile
  /// $E_{\text{rms}}(x_i) = \sqrt{\frac{1}{|W_i|}\sum_{x_j \in W_i} (y_{\text{num}}(x_j)-y_{\text{ref}}(x_j))^2}$
  pub sliding_rms: Array1<f64>,
  /// Sliding window maximum error profile
  /// $E_{\max}(x_i) = \max_{x_j \in W_i} |y_{\text{num}}(x_j) - y_{\text{ref}}(x_j)|$
  pub sliding_max: Array1<f64>,
  /// Cumulative error curve: $C(x_i) = \sum_{j \le i} (y_{\text{num}}(x_j) - y_{\text{ref}}(x_j))\Delta x$
  pub cumulative_error: Array1<f64>,
  /// Summary statistics
  pub summary: SpatialMetricsSummary,
}

/// Summary statistics for spatial metrics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpatialMetricsSummary {
  /// Peak region metrics
  pub peak_region: RegionSummary,
  /// Tail region metrics
  pub tail_region: RegionSummary,
  /// Maximum sliding RMS error
  pub sliding_rms_max: f64,
  /// Maximum sliding max error
  pub sliding_max_max: f64,
  /// Final cumulative error value
  pub cumulative_final: f64,
  /// Maximum absolute cumulative error
  pub cumulative_max_abs: f64,
}

/// Summary statistics for a specific region
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RegionSummary {
  /// Number of points in region
  pub point_count: usize,
  /// Mean error in region
  pub mean_error: f64,
  /// Maximum error in region
  pub max_error: f64,
  /// Standard deviation of errors in region
  pub std_error: f64,
}

/// Error distribution analysis including histograms and masks
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ErrorDistributions {
  /// Histogram of log10(absolute errors)
  pub abs_error_histogram: ErrorHistogram,
  /// Support coverage mask (1.0 where support mismatch detected, 0.0 otherwise)
  /// Flags where one output has significant value but other is near zero
  pub support_coverage_mask: Array1<f64>,
  /// Tolerance pass mask at three levels (1.0 = pass, 0.0 = fail)
  /// Point passes if it meets absolute OR relative tolerance at that level
  pub tolerance_pass_masks: [Array1<f64>; 3],
  /// Summary statistics
  pub summary: ErrorDistributionSummary,
}

/// Histogram representation for error distributions
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ErrorHistogram {
  /// Bin edges (length = bin_counts.len() + 1)
  pub bin_edges: Vec<f64>,
  /// Count of errors in each bin
  pub bin_counts: Vec<usize>,
  /// Bin centers for plotting
  pub bin_centers: Vec<f64>,
}

/// Summary statistics for error distributions
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ErrorDistributionSummary {
  /// Histogram statistics
  pub histogram_total_count: usize,
  /// Number of points with support mismatch
  pub support_mismatch_count: usize,
  /// Fractions passing tolerance at each level [strict, moderate, loose]
  pub tolerance_pass_fractions: [f64; 3],
  /// Statistical summary
  pub stats: DistributionStats,
}

/// Statistical summary of error distribution
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DistributionStats {
  /// Median absolute error
  pub median: f64,
  /// 25th percentile
  pub q25: f64,
  /// 75th percentile
  pub q75: f64,
  /// 95th percentile
  pub q95: f64,
  /// 99th percentile
  pub q99: f64,
}

/// Validates that all input arrays have compatible dimensions
fn validate_inputs(x: &Array1<f64>, actual: &Array1<f64>, expected: &Array1<f64>) -> eyre::Result<()> {
  let n = x.len();
  if n < 2 {
    return make_error!("Input arrays must have at least 2 points, got {n}");
  }
  if actual.len() != n {
    return make_error!("Actual array length {} does not match grid length {n}", actual.len());
  }
  if expected.len() != n {
    return make_error!(
      "Expected array length {} does not match grid length {n}",
      expected.len()
    );
  }
  Ok(())
}

/// Computes grid spacing (assumes uniform grid)
fn compute_grid_spacing(x: &Array1<f64>) -> eyre::Result<f64> {
  if x.len() < 2 {
    return make_error!("Cannot compute grid spacing with fewer than 2 points");
  }
  let dx = x[1] - x[0];
  if !dx.is_finite() || dx <= 0.0 {
    return make_error!("Invalid grid spacing: {dx}");
  }
  Ok(dx)
}

/// Computes all pointwise error metrics
fn compute_pointwise_errors(
  actual: &Array1<f64>,
  expected: &Array1<f64>,
  config: &PointwiseConfig,
) -> eyre::Result<PointwiseErrors> {
  let n = actual.len();

  let signed = actual - expected;
  let absolute = signed.mapv(|x| x.abs());

  let mut relative = Array1::zeros(n);
  for i in 0..n {
    let denom = expected[i].abs().max(config.epsilon);
    relative[i] = signed[i].abs() / denom;
  }

  let mut logarithmic = Array1::from_elem(n, f64::NAN);
  let mut log_valid_count = 0;
  let mut log_max = 0.0;
  for i in 0..n {
    if actual[i] > config.log_threshold && expected[i] > config.log_threshold {
      let log_err = (actual[i].ln() - expected[i].ln()).abs();
      logarithmic[i] = log_err;
      log_valid_count += 1;
      if log_err > log_max {
        log_max = log_err;
      }
    }
  }

  let abs_mean = absolute.mean().unwrap_or(0.0);
  let abs_max = absolute.iter().copied().fold(0.0_f64, f64::max);
  let abs_variance = absolute.mapv(|x| (x - abs_mean).powi(2)).mean().unwrap_or(0.0);
  let abs_std = abs_variance.sqrt();

  let rel_mean = relative.mean().unwrap_or(0.0);
  let rel_max = relative.iter().copied().fold(0.0_f64, f64::max);
  let mut rel_sorted: Vec<f64> = relative.iter().copied().collect();
  rel_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
  let rel_median = compute_median(&rel_sorted);

  let signed_bias = signed.mean().unwrap_or(0.0);

  let summary = PointwiseErrorSummary {
    abs_mean,
    abs_max,
    abs_std,
    rel_mean,
    rel_max,
    rel_median,
    signed_bias,
    log_valid_count,
    log_max,
  };

  Ok(PointwiseErrors {
    absolute,
    relative,
    signed,
    logarithmic,
    summary,
  })
}

/// Computes structural error metrics including derivatives and symmetry
fn compute_structural_errors(
  x: &Array1<f64>,
  actual: &Array1<f64>,
  expected: &Array1<f64>,
  dx: f64,
  config: &PointwiseConfig,
) -> eyre::Result<StructuralErrors> {
  let n = x.len();

  let d1_actual = compute_first_derivative(actual, dx);
  let d1_expected = compute_first_derivative(expected, dx);
  let first_derivative = (&d1_actual - &d1_expected).mapv(|x| x.abs());

  let d2_actual = compute_second_derivative(actual, dx);
  let d2_expected = compute_second_derivative(expected, dx);
  let second_derivative = (&d2_actual - &d2_expected).mapv(|x| x.abs());

  let symmetry_residual = compute_symmetry_residual(x, actual);

  let mut monotonicity_violations = Array1::zeros(n);
  let mut violation_count = 0;
  for i in 0..n - 1 {
    if actual[i + 1] - actual[i] < -config.monotonicity_eta {
      monotonicity_violations[i] = 1.0;
      violation_count += 1;
    }
  }

  let d1_max = first_derivative.iter().copied().fold(0.0_f64, f64::max);
  let d1_mean = first_derivative.mean().unwrap_or(0.0);
  let d2_max = second_derivative.iter().copied().fold(0.0_f64, f64::max);
  let d2_mean = second_derivative.mean().unwrap_or(0.0);
  let symmetry_max = symmetry_residual.iter().copied().fold(0.0_f64, |a, b| a.max(b.abs()));

  let summary = StructuralErrorSummary {
    d1_max,
    d1_mean,
    d2_max,
    d2_mean,
    symmetry_max,
    monotonicity_violation_count: violation_count,
  };

  Ok(StructuralErrors {
    first_derivative,
    second_derivative,
    symmetry_residual,
    monotonicity_violations,
    summary,
  })
}

/// Computes first derivative using central differences with one-sided at boundaries
fn compute_first_derivative(y: &Array1<f64>, dx: f64) -> Array1<f64> {
  let n = y.len();
  let mut dy = Array1::zeros(n);

  if n < 2 {
    return dy;
  }

  dy[0] = (y[1] - y[0]) / dx;
  for i in 1..n - 1 {
    dy[i] = (y[i + 1] - y[i - 1]) / (2.0 * dx);
  }
  dy[n - 1] = (y[n - 1] - y[n - 2]) / dx;

  dy
}

/// Computes second derivative using central differences with one-sided at boundaries
fn compute_second_derivative(y: &Array1<f64>, dx: f64) -> Array1<f64> {
  let n = y.len();
  let mut d2y = Array1::zeros(n);

  if n < 3 {
    return d2y;
  }

  let dx2 = dx * dx;

  d2y[0] = (y[2] - 2.0 * y[1] + y[0]) / dx2;
  for i in 1..n - 1 {
    d2y[i] = (y[i + 1] - 2.0 * y[i] + y[i - 1]) / dx2;
  }
  d2y[n - 1] = (y[n - 1] - 2.0 * y[n - 2] + y[n - 3]) / dx2;

  d2y
}

/// Computes symmetry residual for values around x=0
fn compute_symmetry_residual(x: &Array1<f64>, y: &Array1<f64>) -> Array1<f64> {
  let n = x.len();
  let mut residual = Array1::zeros(n);

  let x_map: BTreeMap<String, (usize, f64)> = x
    .iter()
    .enumerate()
    .map(|(i, &val)| (format!("{val:.12}"), (i, val)))
    .collect();

  for (i, &xi) in x.iter().enumerate() {
    let neg_xi = -xi;
    let key = format!("{neg_xi:.12}");
    if let Some(&(j, _)) = x_map.get(&key) {
      residual[i] = y[i] - y[j];
    }
  }

  residual
}

/// Computes spatial and regional metrics
fn compute_spatial_metrics(
  x: &Array1<f64>,
  actual: &Array1<f64>,
  expected: &Array1<f64>,
  dx: f64,
  config: &PointwiseConfig,
) -> eyre::Result<SpatialMetrics> {
  let n = x.len();
  let errors = actual - expected;

  let peak_idx_expected = expected
    .iter()
    .enumerate()
    .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
    .map_or(0, |(i, _)| i);
  let peak_x = x[peak_idx_expected];

  let peak_val = expected.iter().copied().fold(f64::NEG_INFINITY, f64::max);
  let tail_threshold_val = peak_val * config.tail_threshold;

  let mut peak_region_errors = Array1::from_elem(n, f64::NAN);
  let mut peak_errors_valid = Vec::new();
  for i in 0..n {
    if (x[i] - peak_x).abs() <= config.peak_region_radius {
      let rel_err = errors[i].abs() / expected[i].abs().max(config.epsilon);
      peak_region_errors[i] = rel_err;
      peak_errors_valid.push(rel_err);
    }
  }

  let mut tail_region_errors = Array1::from_elem(n, f64::NAN);
  let mut tail_errors_valid = Vec::new();
  for i in 0..n {
    if expected[i] <= tail_threshold_val && actual[i] > config.log_threshold && expected[i] > config.log_threshold {
      let log_err = (actual[i].ln() - expected[i].ln()).abs();
      tail_region_errors[i] = log_err;
      tail_errors_valid.push(log_err);
    }
  }

  let sliding_rms = compute_sliding_window_rms(&errors, config.window_half_width);
  let sliding_max = compute_sliding_window_max(&errors.mapv(|x| x.abs()), config.window_half_width);

  let mut cumulative_error = Array1::zeros(n);
  cumulative_error[0] = errors[0] * dx;
  for i in 1..n {
    cumulative_error[i] = cumulative_error[i - 1] + errors[i] * dx;
  }

  let peak_region = compute_region_summary(&peak_errors_valid);
  let tail_region = compute_region_summary(&tail_errors_valid);

  let sliding_rms_max = sliding_rms.iter().copied().fold(0.0_f64, f64::max);
  let sliding_max_max = sliding_max.iter().copied().fold(0.0_f64, f64::max);
  let cumulative_final = cumulative_error[n - 1];
  let cumulative_max_abs = cumulative_error.iter().copied().fold(0.0_f64, |a, b| a.max(b.abs()));

  let summary = SpatialMetricsSummary {
    peak_region,
    tail_region,
    sliding_rms_max,
    sliding_max_max,
    cumulative_final,
    cumulative_max_abs,
  };

  Ok(SpatialMetrics {
    peak_region_errors,
    tail_region_errors,
    sliding_rms,
    sliding_max,
    cumulative_error,
    summary,
  })
}

/// Computes sliding window RMS error
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

/// Computes sliding window maximum error
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

/// Computes summary statistics for a region
fn compute_region_summary(errors: &[f64]) -> RegionSummary {
  if errors.is_empty() {
    return RegionSummary {
      point_count: 0,
      mean_error: 0.0,
      max_error: 0.0,
      std_error: 0.0,
    };
  }

  let count = errors.len();
  let mean = errors.iter().sum::<f64>() / count as f64;
  let max = errors.iter().copied().fold(0.0_f64, f64::max);
  let variance = errors.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / count as f64;
  let std = variance.sqrt();

  RegionSummary {
    point_count: count,
    mean_error: mean,
    max_error: max,
    std_error: std,
  }
}

/// Computes error distribution metrics including histograms and masks
fn compute_error_distributions(
  x: &Array1<f64>,
  actual: &Array1<f64>,
  expected: &Array1<f64>,
  pointwise_errors: &PointwiseErrors,
  config: &PointwiseConfig,
) -> eyre::Result<ErrorDistributions> {
  let n = x.len();

  let abs_error_histogram = compute_error_histogram(&pointwise_errors.absolute, config.histogram_bins)?;

  let mut support_coverage_mask = Array1::zeros(n);
  let mut support_mismatch_count = 0;
  for i in 0..n {
    let actual_significant = actual[i].abs() > config.log_threshold;
    let expected_significant = expected[i].abs() > config.log_threshold;
    if actual_significant != expected_significant {
      support_coverage_mask[i] = 1.0;
      support_mismatch_count += 1;
    }
  }

  let mut tolerance_pass_masks = [Array1::zeros(n), Array1::zeros(n), Array1::zeros(n)];
  let mut pass_counts = [0_usize; 3];

  for level in 0..3 {
    let abs_tol = config.abs_tolerances[level];
    let rel_tol = config.rel_tolerances[level];
    for i in 0..n {
      let abs_pass = pointwise_errors.absolute[i] <= abs_tol;
      let rel_pass = pointwise_errors.relative[i] <= rel_tol;
      if abs_pass || rel_pass {
        tolerance_pass_masks[level][i] = 1.0;
        pass_counts[level] += 1;
      }
    }
  }

  let tolerance_pass_fractions = [
    pass_counts[0] as f64 / n as f64,
    pass_counts[1] as f64 / n as f64,
    pass_counts[2] as f64 / n as f64,
  ];

  let mut abs_errors_sorted: Vec<f64> = pointwise_errors.absolute.iter().copied().collect();
  abs_errors_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

  let stats = DistributionStats {
    median: compute_median(&abs_errors_sorted),
    q25: compute_quantile(&abs_errors_sorted, 0.25),
    q75: compute_quantile(&abs_errors_sorted, 0.75),
    q95: compute_quantile(&abs_errors_sorted, 0.95),
    q99: compute_quantile(&abs_errors_sorted, 0.99),
  };

  let summary = ErrorDistributionSummary {
    histogram_total_count: abs_error_histogram.bin_counts.iter().sum(),
    support_mismatch_count,
    tolerance_pass_fractions,
    stats,
  };

  Ok(ErrorDistributions {
    abs_error_histogram,
    support_coverage_mask,
    tolerance_pass_masks,
    summary,
  })
}

/// Computes histogram of errors using log10 scale
fn compute_error_histogram(errors: &Array1<f64>, num_bins: usize) -> eyre::Result<ErrorHistogram> {
  if num_bins == 0 {
    return make_error!("Number of histogram bins must be positive");
  }

  let valid_errors: Vec<f64> = errors.iter().copied().filter(|&x| x > 0.0 && x.is_finite()).collect();

  if valid_errors.is_empty() {
    return Ok(ErrorHistogram {
      bin_edges: vec![0.0; num_bins + 1],
      bin_counts: vec![0; num_bins],
      bin_centers: vec![0.0; num_bins],
    });
  }

  let min_log = valid_errors
    .iter()
    .copied()
    .map(|x| x.log10())
    .fold(f64::INFINITY, f64::min);
  let max_log = valid_errors
    .iter()
    .copied()
    .map(|x| x.log10())
    .fold(f64::NEG_INFINITY, f64::max);

  let range = max_log - min_log;
  let bin_width = if range > 0.0 { range / num_bins as f64 } else { 1.0 };

  let mut bin_edges = Vec::with_capacity(num_bins + 1);
  for i in 0..=num_bins {
    bin_edges.push(min_log + i as f64 * bin_width);
  }

  let mut bin_counts = vec![0_usize; num_bins];
  for &err in &valid_errors {
    let log_err = err.log10();
    let bin_idx = ((log_err - min_log) / bin_width).floor() as usize;
    let bin_idx = bin_idx.min(num_bins - 1);
    bin_counts[bin_idx] += 1;
  }

  let bin_centers: Vec<f64> = (0..num_bins)
    .map(|i| f64::midpoint(bin_edges[i], bin_edges[i + 1]))
    .collect();

  Ok(ErrorHistogram {
    bin_edges,
    bin_counts,
    bin_centers,
  })
}

/// Computes median of a sorted array
fn compute_median(sorted: &[f64]) -> f64 {
  if sorted.is_empty() {
    return 0.0;
  }
  let n = sorted.len();
  if n.is_multiple_of(2) {
    f64::midpoint(sorted[n / 2 - 1], sorted[n / 2])
  } else {
    sorted[n / 2]
  }
}

/// Computes quantile of a sorted array
fn compute_quantile(sorted: &[f64], q: f64) -> f64 {
  if sorted.is_empty() {
    return 0.0;
  }
  let n = sorted.len();
  let index = (q * (n - 1) as f64).round() as usize;
  sorted[index.min(n - 1)]
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_abs_diff_eq;
  use ndarray::array;

  #[test]
  fn test_perfect_agreement() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y = array![1.0, 2.0, 3.0, 2.0, 1.0];
    let result = PointwiseMetrics::new(&x, &y, &y).unwrap();

    assert_eq!(result.total_points, 5);
    assert_abs_diff_eq!(result.pointwise_errors.summary.abs_max, 0.0, epsilon = 1e-12);
    assert_abs_diff_eq!(result.pointwise_errors.summary.rel_max, 0.0, epsilon = 1e-12);
    assert_abs_diff_eq!(result.pointwise_errors.summary.signed_bias, 0.0, epsilon = 1e-12);
  }

  #[test]
  fn test_constant_offset() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let expected = array![1.0, 2.0, 3.0, 2.0, 1.0];
    let actual = &expected + 0.5;
    let result = PointwiseMetrics::new(&x, &actual, &expected).unwrap();

    for i in 0..5 {
      assert_abs_diff_eq!(result.pointwise_errors.absolute[i], 0.5, epsilon = 1e-12);
      assert_abs_diff_eq!(result.pointwise_errors.signed[i], 0.5, epsilon = 1e-12);
    }
    assert_abs_diff_eq!(result.pointwise_errors.summary.signed_bias, 0.5, epsilon = 1e-12);
  }

  #[test]
  fn test_scale_factor() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let expected = array![1.0, 2.0, 3.0, 2.0, 1.0];
    let actual = &expected * 1.1;
    let result = PointwiseMetrics::new(&x, &actual, &expected).unwrap();

    for i in 0..5 {
      let expected_rel = 0.1;
      assert_abs_diff_eq!(result.pointwise_errors.relative[i], expected_rel, epsilon = 1e-10);
    }
  }

  #[test]
  fn test_grid_spacing_uniform() {
    let x = array![0.0, 0.5, 1.0, 1.5, 2.0];
    let dx = compute_grid_spacing(&x).unwrap();
    assert_abs_diff_eq!(dx, 0.5, epsilon = 1e-12);
  }

  #[test]
  fn test_first_derivative_linear() {
    let y = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let dx = 1.0;
    let dy = compute_first_derivative(&y, dx);

    for i in 0..5 {
      assert_abs_diff_eq!(dy[i], 1.0, epsilon = 1e-10);
    }
  }

  #[test]
  fn test_second_derivative_quadratic() {
    let y = array![0.0, 1.0, 4.0, 9.0, 16.0];
    let dx = 1.0;
    let d2y = compute_second_derivative(&y, dx);

    for i in 1..4 {
      assert_abs_diff_eq!(d2y[i], 2.0, epsilon = 1e-10);
    }
  }

  #[test]
  fn test_symmetry_perfect() {
    let x = array![-2.0, -1.0, 0.0, 1.0, 2.0];
    let y = array![4.0, 1.0, 0.0, 1.0, 4.0];
    let residual = compute_symmetry_residual(&x, &y);

    for i in 0..5 {
      assert_abs_diff_eq!(residual[i], 0.0, epsilon = 1e-12);
    }
  }

  #[test]
  fn test_symmetry_asymmetric() {
    let x = array![-2.0, -1.0, 0.0, 1.0, 2.0];
    let y = array![3.0, 1.0, 0.0, 1.0, 4.0];
    let residual = compute_symmetry_residual(&x, &y);

    assert_abs_diff_eq!(residual[0].abs(), 1.0, epsilon = 1e-12);
    assert_abs_diff_eq!(residual[4].abs(), 1.0, epsilon = 1e-12);
  }

  #[test]
  fn test_monotonicity_violations() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let actual = array![1.0, 2.0, 1.5, 3.0, 2.5];
    let expected = array![1.0, 2.0, 3.0, 4.0, 5.0];

    let config = PointwiseConfig::default();
    let result = compute_structural_errors(&x, &actual, &expected, 1.0, &config).unwrap();

    assert!(result.summary.monotonicity_violation_count > 0);
  }

  #[test]
  fn test_cumulative_error_cancellation() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let expected = array![0.0, 0.0, 0.0, 0.0, 0.0];
    let actual = array![1.0, -1.0, 1.0, -1.0, 0.0];

    let result = PointwiseMetrics::new(&x, &actual, &expected).unwrap();

    assert!(result.spatial_metrics.summary.cumulative_max_abs < 2.0);
  }

  #[test]
  fn test_cumulative_error_systematic_bias() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let expected = array![0.0, 0.0, 0.0, 0.0, 0.0];
    let actual = array![0.1, 0.1, 0.1, 0.1, 0.1];

    let result = PointwiseMetrics::new(&x, &actual, &expected).unwrap();

    assert!(result.spatial_metrics.summary.cumulative_final > 0.3);
  }

  #[test]
  fn test_tolerance_masks() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let expected = array![1.0, 1.0, 1.0, 1.0, 1.0];
    let actual = array![1.0, 1.0 + 1e-7, 1.0 + 1e-10, 1.0 + 1e-13, 1.0];

    let result = PointwiseMetrics::new(&x, &actual, &expected).unwrap();

    assert_abs_diff_eq!(
      result.error_distributions.tolerance_pass_masks[2].sum(),
      5.0,
      epsilon = 1e-12
    );
    assert!(result.error_distributions.tolerance_pass_masks[1].sum() >= 3.0);
    assert!(result.error_distributions.tolerance_pass_masks[0].sum() >= 2.0);
  }

  #[test]
  fn test_support_coverage_mismatch() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let expected = array![1.0, 1e-12, 1.0, 1e-12, 1.0];
    let actual = array![1.0, 1.0, 1e-12, 1e-12, 1.0];

    let result = PointwiseMetrics::new(&x, &actual, &expected).unwrap();

    assert_eq!(result.error_distributions.summary.support_mismatch_count, 2);
  }

  #[test]
  fn test_sliding_window_rms() {
    let errors = array![1.0, 0.0, 0.0, 0.0, 1.0];
    let rms = compute_sliding_window_rms(&errors, 1);

    assert!(rms[0] > rms[2]);
    assert!(rms[4] > rms[2]);
  }

  #[test]
  fn test_error_histogram_uniform() {
    let errors = Array1::linspace(1e-6, 1e-3, 100);
    let histogram = compute_error_histogram(&errors, 10).unwrap();

    assert_eq!(histogram.bin_counts.len(), 10);
    assert_eq!(histogram.bin_edges.len(), 11);
    assert_eq!(histogram.bin_centers.len(), 10);

    let total_count: usize = histogram.bin_counts.iter().sum();
    assert_eq!(total_count, 100);
  }

  #[test]
  fn test_quantile_computation() {
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    assert_abs_diff_eq!(compute_median(&data), 3.0, epsilon = 1e-12);
    assert_abs_diff_eq!(compute_quantile(&data, 0.0), 1.0, epsilon = 1e-12);
    assert_abs_diff_eq!(compute_quantile(&data, 1.0), 5.0, epsilon = 1e-12);
    assert_abs_diff_eq!(compute_quantile(&data, 0.5), 3.0, epsilon = 1e-12);
  }

  #[test]
  fn test_validation_mismatched_lengths() {
    let x = array![0.0, 1.0, 2.0];
    let actual = array![1.0, 2.0];
    let expected = array![1.0, 2.0, 3.0];

    let result = PointwiseMetrics::new(&x, &actual, &expected);
    drop(result.unwrap_err());
  }

  #[test]
  fn test_validation_insufficient_points() {
    let x = array![0.0];
    let y = array![1.0];

    let result = PointwiseMetrics::new(&x, &y, &y);
    drop(result.unwrap_err());
  }

  #[test]
  fn test_empty_region_summary() {
    let summary = compute_region_summary(&[]);
    assert_eq!(summary.point_count, 0);
    assert_abs_diff_eq!(summary.mean_error, 0.0, epsilon = 1e-12);
    assert_abs_diff_eq!(summary.max_error, 0.0, epsilon = 1e-12);
    assert_abs_diff_eq!(summary.std_error, 0.0, epsilon = 1e-12);
  }

  #[test]
  fn test_log_error_with_zeros() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let expected = array![0.0, 1.0, 2.0, 1.0, 0.0];
    let actual = array![0.0, 1.1, 2.2, 1.1, 0.0];

    let result = PointwiseMetrics::new(&x, &actual, &expected).unwrap();

    assert!(result.pointwise_errors.summary.log_valid_count < 5);
    assert!(result.pointwise_errors.summary.log_max.is_finite());
  }

  #[test]
  fn test_custom_config() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let expected = array![1.0, 2.0, 3.0, 2.0, 1.0];
    let actual = &expected + 0.01;

    let config = PointwiseConfig {
      abs_tolerances: [1e-3, 1e-6, 1e-9],
      ..Default::default()
    };

    let result = PointwiseMetrics::new_with_config(&x, &actual, &expected, &config).unwrap();

    assert!(result.error_distributions.summary.tolerance_pass_fractions[0] > 0.5);
  }
}
