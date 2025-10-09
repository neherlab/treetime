use crate::make_error;
use crate::utils::float_fmt::float_to_digits;
use itertools::{Itertools, izip};
use ndarray::Array1;
use ndarray_stats::QuantileExt;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::cmp::Ordering;
use std::fmt;

/// Comprehensive domain-wide agreement metrics between actual and expected solutions
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DomainAgreementMetrics {
  /// Total number of evaluation points
  pub total_points: usize,
  /// Absolute error statistics including bias detection
  pub abs_error_stats: AbsoluteErrorStats,
  /// Relative error statistics with robust central tendency measures
  pub rel_error_stats: RelativeErrorStats,
  /// Global quality metrics including conservation properties
  pub quality_metrics: QualityMetrics,
  /// Peak-specific accuracy metrics for distribution analysis
  pub peak_metrics: PeakMetrics,
  /// Location information for the maximum absolute error
  pub max_error_location: MaxErrorLocation,
  /// Fractions of points within absolute tolerance thresholds [0.0, 1.0]
  pub abs_tolerance_fractions: [f64; 3],
  /// Fractions of points within relative tolerance thresholds [0.0, 1.0]
  pub rel_tolerance_fractions: [f64; 3],
  /// Overall assessment based on R² value
  pub overall_assessment: AgreementAssessment,
}

impl DomainAgreementMetrics {
  /// Creates new domain agreement metrics from actual and expected values
  pub fn new(x: &Array1<f64>, actual: &Array1<f64>, expected: &Array1<f64>) -> eyre::Result<Self> {
    Self::new_with_thresholds(x, actual, expected, &ToleranceThresholds::default())
  }

  /// Creates new domain agreement metrics with custom tolerance thresholds
  pub fn new_with_thresholds(
    x: &Array1<f64>,
    actual: &Array1<f64>,
    expected: &Array1<f64>,
    thresholds: &ToleranceThresholds,
  ) -> eyre::Result<Self> {
    let total_points = actual.len();

    if actual.len() != expected.len() || actual.len() != x.len() {
      return make_error!(
        "Input arrays must have same length: actual={}, expected={}, x={}",
        actual.len(),
        expected.len(),
        x.len()
      );
    }

    if total_points == 0 {
      return make_error!("Input arrays cannot be empty");
    }

    let abs_error_stats = compute_absolute_error_statistics(actual, expected);
    let rel_error_stats = compute_relative_error_statistics(actual, expected)?;
    let rmse = compute_rmse(actual, expected);
    let r_squared = compute_r_squared(actual, expected)?;
    let correlation = compute_correlation(actual, expected)?;
    let mass_error = compute_mass_error(x, actual, expected);
    let rel_l2_error = compute_relative_l2_norm_error(actual, expected)?;
    let rel_l1_error = compute_relative_l1_norm_error(actual, expected)?;
    let rel_linf_error = compute_relative_linf_norm_error(actual, expected)?;
    let max_log_error = compute_max_log_error(actual, expected, 1e-300)?;
    let symmetry_error = compute_symmetry_error(x, actual);
    let quantile_95_error = compute_quantile_error(actual, expected, 0.95);
    let peak_metrics = compute_peak_metrics(x, actual, expected)?;

    let quality_metrics = QualityMetrics {
      rmse,
      r_squared,
      correlation,
      mass_error,
      rel_l2_error,
      rel_l1_error,
      rel_linf_error,
      max_log_error,
      symmetry_error,
      quantile_95_error,
    };

    let tolerance_counts = compute_tolerance_counts(actual, expected, thresholds);
    let max_error_location = find_max_error_location(x, actual, expected);

    let overall_assessment = compute_overall_assessment(quality_metrics.r_squared, &thresholds.r2_thresholds);

    let abs_tolerance_fractions = [
      tolerance_counts.within_abs_tolerances[0] as f64 / total_points as f64,
      tolerance_counts.within_abs_tolerances[1] as f64 / total_points as f64,
      tolerance_counts.within_abs_tolerances[2] as f64 / total_points as f64,
    ];

    let rel_tolerance_fractions = [
      tolerance_counts.within_rel_tolerances[0] as f64 / total_points as f64,
      tolerance_counts.within_rel_tolerances[1] as f64 / total_points as f64,
      tolerance_counts.within_rel_tolerances[2] as f64 / total_points as f64,
    ];

    Ok(Self {
      total_points,
      abs_error_stats,
      rel_error_stats,
      quality_metrics,
      peak_metrics,
      max_error_location,
      abs_tolerance_fractions,
      rel_tolerance_fractions,
      overall_assessment,
    })
  }

  /// Returns the precomputed overall assessment based on R² value
  pub fn overall_assessment(&self) -> AgreementAssessment {
    self.overall_assessment
  }

  /// Returns fraction of points within absolute tolerance at given level (0-2)
  pub fn abs_tolerance_fraction(&self, level: usize) -> f64 {
    if level >= 3 {
      return 0.0;
    }
    self.abs_tolerance_fractions[level]
  }

  /// Returns fraction of points within relative tolerance at given level (0-2)
  pub fn rel_tolerance_fraction(&self, level: usize) -> f64 {
    if level >= 3 {
      return 0.0;
    }
    self.rel_tolerance_fractions[level]
  }
}

/// Assessment levels for domain agreement quality
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AgreementAssessment {
  Excellent,
  VeryGood,
  Good,
  Poor,
}

impl fmt::Display for AgreementAssessment {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    match self {
      AgreementAssessment::Excellent => write!(f, "🟢 EXCELLENT: Near-perfect agreement"),
      AgreementAssessment::VeryGood => write!(f, "🟡 VERY GOOD: High agreement"),
      AgreementAssessment::Good => write!(f, "🟠 GOOD: Reasonable agreement"),
      AgreementAssessment::Poor => write!(f, "🔴 POOR: Low agreement"),
    }
  }
}

/// Absolute error statistics including bias detection
#[derive(Debug, Clone, Serialize, Deserialize)]
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
#[derive(Debug, Clone, Serialize, Deserialize)]
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

/// Quality metrics for agreement assessment including conservation properties
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityMetrics {
  /// Root Mean Square Error: $\sqrt{\frac{1}{N}\sum_i (y_{\text{num}}(x_i) - y_{\text{ref}}(x_i))^2}$
  /// Global measure of deviation magnitude, sensitive to large errors
  pub rmse: f64,
  /// Coefficient of determination: $R^2 = 1 - \frac{SS_{\text{res}}}{SS_{\text{tot}}}$
  /// Fraction of variance explained, with $SS_{\text{res}} = \sum_i (y_{\text{num}}(x_i) - y_{\text{ref}}(x_i))^2$
  pub r_squared: f64,
  /// Pearson correlation coefficient: measures linear relationship strength
  /// $r = \frac{\sum_i (y_{\text{num}}(x_i) - \bar{y}_{\text{num}})(y_{\text{ref}}(x_i) - \bar{y}_{\text{ref}})}{\sqrt{\sum_i (y_{\text{num}}(x_i) - \bar{y}_{\text{num}})^2 \sum_i (y_{\text{ref}}(x_i) - \bar{y}_{\text{ref}})^2}}$
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
  /// Sensitive to tail accuracy for small values, uses threshold τ to avoid log(0)
  pub max_log_error: f64,
  /// Symmetry error: $\max_i \left|y_{\text{num}}(x_i) - y_{\text{num}}(-x_i)\right|$
  /// For symmetric kernels (e.g., Gaussians), measures deviation from symmetry
  pub symmetry_error: f64,
  /// 95th percentile error threshold: value below which 95% of absolute errors fall
  /// Quantile-based robustness measure for error distribution analysis
  pub quantile_95_error: f64,
}

/// Peak-related accuracy metrics for distribution analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PeakMetrics {
  /// Peak value error: $\frac{\left|\max_i y_{\text{num}}(x_i) - \max_i y_{\text{ref}}(x_i)\right|}{\max_i y_{\text{ref}}(x_i)}$
  /// Relative error in maximum amplitude, critical for probability distributions
  pub value_error: f64,
  /// Peak location error: $\left|x_{\operatorname{argmax}(y_{\text{num}})} - x_{\operatorname{argmax}(y_{\text{ref}})}\right|$
  /// Difference in location of the peak, important for temporal accuracy
  pub location_error: f64,
}

/// Tolerance counts for different threshold levels
/// Provides pass/fail statistics at multiple precision requirements
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ToleranceCounts {
  /// Count of points within absolute tolerance thresholds [strict, moderate, loose]
  /// Points where $\left|y_{\text{num}}(x_i) - y_{\text{ref}}(x_i)\right| < \text{threshold}$
  pub within_abs_tolerances: [usize; 3],
  /// Count of points within relative tolerance thresholds [strict, moderate, loose]
  /// Points where $\frac{\left|y_{\text{num}}(x_i) - y_{\text{ref}}(x_i)\right|}{\left|y_{\text{ref}}(x_i)\right|} < \text{threshold}$
  pub within_rel_tolerances: [usize; 3],
}

/// Maximum error location information for error analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MaxErrorLocation {
  /// Array index where maximum absolute error occurs
  pub idx: usize,
  /// X-coordinate value where maximum absolute error occurs
  /// Useful for identifying problematic regions in the domain
  pub x_value: f64,
}

/// Tolerance thresholds for agreement assessment at multiple precision levels
#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
pub struct ToleranceThresholds {
  /// Absolute tolerance thresholds: [strict, moderate, loose]
  /// Default: [1e-6, 1e-9, 1e-12] for high-precision numerical analysis
  #[default(_code = "[1e-6, 1e-9, 1e-12]")]
  pub abs_tolerances: [f64; 3],
  /// Relative tolerance thresholds: [strict, moderate, loose]
  /// Default: [1%, 0.1%, 0.01%] for percentage-based accuracy requirements
  #[default(_code = "[0.01, 0.001, 0.0001]")] // 1%, 0.1%, 0.01%
  pub rel_tolerances: [f64; 3],
  /// R² thresholds for overall assessment: [excellent, very_good, good]
  /// Default: [0.999999, 0.9999, 0.99] for very high correlation requirements
  #[default(_code = "[0.999999, 0.9999, 0.99]")] // [excellent, very_good, good]
  pub r2_thresholds: [f64; 3],
}

/// Compute overall assessment based on R² value and thresholds
fn compute_overall_assessment(r2: f64, thresholds: &[f64; 3]) -> AgreementAssessment {
  if r2 >= thresholds[0] {
    AgreementAssessment::Excellent
  } else if r2 >= thresholds[1] {
    AgreementAssessment::VeryGood
  } else if r2 >= thresholds[2] {
    AgreementAssessment::Good
  } else {
    AgreementAssessment::Poor
  }
}

impl fmt::Display for DomainAgreementMetrics {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    let Self {
      total_points,
      abs_error_stats:
        AbsoluteErrorStats {
          mean: mean_abs_error,
          max: max_abs_error,
          std: std_abs_error,
          bias,
        },
      rel_error_stats:
        RelativeErrorStats {
          mean: mean_rel_error,
          max: max_rel_error,
          median: median_rel_error,
          mape,
        },
      quality_metrics:
        QualityMetrics {
          rmse,
          mass_error,
          rel_l2_error,
          rel_l1_error,
          rel_linf_error,
          max_log_error,
          symmetry_error,
          quantile_95_error,
          r_squared,
          correlation,
        },
      peak_metrics:
        PeakMetrics {
          value_error: peak_value_error,
          location_error: peak_location_error,
        },
      max_error_location: MaxErrorLocation {
        x_value: max_error_x_value,
        idx: max_error_idx,
      },
      abs_tolerance_fractions,
      rel_tolerance_fractions,
      overall_assessment,
    } = self;

    let mean_rel_error_pct = mean_rel_error * 100.0;
    let max_rel_error_pct = max_rel_error * 100.0;
    let median_rel_error_pct = median_rel_error * 100.0;
    let peak_value_error_pct = peak_value_error * 100.0;

    // Recreate tolerance thresholds for display
    let abs_tolerances = [1e-6, 1e-9, 1e-12];
    let rel_tolerances = [0.01, 0.001, 0.0001];

    let abs_tolerance_lines = (0..3)
      .map(|i| {
        let count = (abs_tolerance_fractions[i] * *total_points as f64).round() as usize;
        let percentage = abs_tolerance_fractions[i] * 100.0;
        format!(
          "  < {:.0e}: {:3}/{total_points} ({:5.1}%)",
          abs_tolerances[i], count, percentage
        )
      })
      .join("\n");

    let rel_tolerance_lines = (0..3)
      .map(|i| {
        let count = (rel_tolerance_fractions[i] * *total_points as f64).round() as usize;
        let percentage_threshold = rel_tolerances[i] * 100.0;
        let formatted_percentage = float_to_digits(percentage_threshold, Some(3), None);
        let fraction_percentage = rel_tolerance_fractions[i] * 100.0;
        format!("  < {formatted_percentage:>6}%: {count:3}/{total_points} ({fraction_percentage:5.1}%)")
      })
      .join("\n");

    write!(
      f,
      r#"=== Domain-Wide Agreement Metrics ===
Total evaluation points: {total_points}

Absolute Error Statistics:
  Mean absolute error:    {mean_abs_error:.6e}
  Max absolute error:     {max_abs_error:.6e}
  Std absolute error:     {std_abs_error:.6e}
  Signed error bias:      {bias:.6e}
  RMSE:                   {rmse:.6e}
  Max error at x:         {max_error_x_value:.3} (index: {max_error_idx})

Relative Error Statistics:
  Mean relative error:    {mean_rel_error:.6e} ({mean_rel_error_pct:.4}%)
  Max relative error:     {max_rel_error:.6e} ({max_rel_error_pct:.4}%)
  Median relative error:  {median_rel_error:.6e} ({median_rel_error_pct:.4}%)
  MAPE:                   {mape:.6}%

Conservation & Global Quality:
  Mass (integral) error:             {mass_error:.6e}
  Relative L2 norm error:            {rel_l2_error:.6e}
  Relative L1 norm error:            {rel_l1_error:.6e}
  Relative L∞ norm error:            {rel_linf_error:.6e}
  Maximum log error:                 {max_log_error:.6e}
  Symmetry error:                    {symmetry_error:.6e}
  95th percentile error:             {quantile_95_error:.6e}
  R² (coefficient of determination): {r_squared:.12}
  Correlation coefficient:           {correlation:.12}

Peak Accuracy Metrics:
  Peak value error:       {peak_value_error:.6e} ({peak_value_error_pct:.4}%)
  Peak location error:    {peak_location_error:.6e}

Points within absolute tolerance:
{abs_tolerance_lines}

Points within relative tolerance:
{rel_tolerance_lines}

{overall_assessment}"#
    )
  }
}

/// Compute absolute error statistics with bias detection
/// Calculates mean, maximum, standard deviation, and signed bias of absolute errors
/// - Mean: $\frac{1}{N}\sum_i \left|y_{\text{num}}(x_i) - y_{\text{ref}}(x_i)\right|$
/// - Max: $\max_i \left|y_{\text{num}}(x_i) - y_{\text{ref}}(x_i)\right|$
/// - Bias: $\frac{1}{N}\sum_i (y_{\text{num}}(x_i) - y_{\text{ref}}(x_i))$ (signed, detects systematic errors)
fn compute_absolute_error_statistics(actual: &Array1<f64>, expected: &Array1<f64>) -> AbsoluteErrorStats {
  let errors = actual - expected;
  let abs_errors = errors.mapv(|x| x.abs());

  let mean = abs_errors.mean().unwrap_or(0.0);
  let max = *abs_errors.max().unwrap_or(&0.0);
  let bias = errors.mean().unwrap_or(0.0);

  let variance = abs_errors.mapv(|x| (x - mean).powi(2)).mean().unwrap_or(0.0);
  let std = variance.sqrt();

  AbsoluteErrorStats { mean, max, std, bias }
}

/// Compute relative error statistics with robust central tendency measures
/// Calculates mean, maximum, MAPE, and median relative errors with zero-protection
/// - Mean relative: $\frac{1}{N}\sum_i \frac{y_{\text{num}}(x_i) - y_{\text{ref}}(x_i)}{y_{\text{ref}}(x_i)}$
/// - Max relative: $\max_i \frac{\left|y_{\text{num}}(x_i) - y_{\text{ref}}(x_i)\right|}{\left|y_{\text{ref}}(x_i)\right|}$
/// - MAPE: Mean Absolute Percentage Error as percentage
/// - Median: More robust to outliers than mean relative error
fn compute_relative_error_statistics(actual: &Array1<f64>, expected: &Array1<f64>) -> eyre::Result<RelativeErrorStats> {
  let mut rel_errors = Vec::new();
  let mut abs_percentage_errors = Vec::new();
  let mut abs_rel_errors = Vec::new();

  for (&a, &e) in izip!(actual, expected) {
    if e.abs() < f64::EPSILON {
      return make_error!("Expected value too close to zero: {:.2e}", e);
    }

    let rel_error = (a - e) / e;
    let abs_rel_error = rel_error.abs();
    let abs_percentage_error = abs_rel_error * 100.0;

    rel_errors.push(rel_error);
    abs_percentage_errors.push(abs_percentage_error);
    abs_rel_errors.push(abs_rel_error);
  }

  let mean = rel_errors.iter().sum::<f64>() / rel_errors.len() as f64;
  let max = rel_errors.iter().map(|x| x.abs()).fold(0.0_f64, |acc, x| acc.max(x));
  let mape = abs_percentage_errors.iter().sum::<f64>() / abs_percentage_errors.len() as f64;

  // Compute median relative error
  abs_rel_errors.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
  let median = if abs_rel_errors.len() % 2 == 0 {
    let mid = abs_rel_errors.len() / 2;
    f64::midpoint(abs_rel_errors[mid - 1], abs_rel_errors[mid])
  } else {
    abs_rel_errors[abs_rel_errors.len() / 2]
  };

  Ok(RelativeErrorStats {
    mean,
    max,
    mape,
    median,
  })
}

/// Compute Root Mean Square Error
/// Formula: $\mathrm{RMSE} = \sqrt{\frac{1}{N}\sum_i (y_{\text{num}}(x_i) - y_{\text{ref}}(x_i))^2}$
/// Global measure sensitive to large deviations
fn compute_rmse(actual: &Array1<f64>, expected: &Array1<f64>) -> f64 {
  let squared_errors = (actual - expected).mapv(|x| x.powi(2));
  let mse = squared_errors.mean().unwrap_or(0.0);
  mse.sqrt()
}

/// Compute coefficient of determination (R²)
/// Formula: $R^2 = 1 - \frac{SS_{\text{res}}}{SS_{\text{tot}}}$ where:
/// - $SS_{\mathrm{res}} = \sum_i (y_{\text{num}}(x_i) - y_{\text{ref}}(x_i))^2$ (residual sum of squares)
/// - $SS_{\mathrm{tot}} = \sum_i (y_{\text{ref}}(x_i) - \bar{y}_{\text{ref}})^2$ (total sum of squares)
fn compute_r_squared(actual: &Array1<f64>, expected: &Array1<f64>) -> eyre::Result<f64> {
  let expected_mean = expected.mean().unwrap_or(0.0);

  let ss_res = (actual - expected).mapv(|x| x.powi(2)).sum();
  let ss_tot = expected.mapv(|x| (x - expected_mean).powi(2)).sum();

  if ss_tot.abs() < f64::EPSILON {
    return make_error!("Total sum of squares is zero - cannot compute R²");
  }

  Ok(1.0 - (ss_res / ss_tot))
}

/// Compute Pearson correlation coefficient
/// Formula: $r = \frac{\sum_i (y_{\text{num}}(x_i) - \bar{y}_{\text{num}})(y_{\text{ref}}(x_i) - \bar{y}_{\text{ref}})}{\sqrt{\sum_i (y_{\text{num}}(x_i) - \bar{y}_{\text{num}})^2 \sum_i (y_{\text{ref}}(x_i) - \bar{y}_{\text{ref}})^2}}$
/// Measures strength of linear relationship between numerical and reference values
fn compute_correlation(actual: &Array1<f64>, expected: &Array1<f64>) -> eyre::Result<f64> {
  let actual_mean = actual.mean().unwrap_or(0.0);
  let expected_mean = expected.mean().unwrap_or(0.0);

  let numerator = izip!(actual, expected)
    .map(|(&a, &e)| (a - actual_mean) * (e - expected_mean))
    .sum::<f64>();

  let actual_var = actual.mapv(|x| (x - actual_mean).powi(2)).sum();
  let expected_var = expected.mapv(|x| (x - expected_mean).powi(2)).sum();

  let denominator = (actual_var * expected_var).sqrt();

  if denominator.abs() < f64::EPSILON {
    return make_error!("Cannot compute correlation - zero variance detected");
  }

  Ok(numerator / denominator)
}

/// Compute tolerance counts for both absolute and relative thresholds
/// Counts points meeting precision requirements at three levels: strict, moderate, loose
fn compute_tolerance_counts(
  actual: &Array1<f64>,
  expected: &Array1<f64>,
  thresholds: &ToleranceThresholds,
) -> ToleranceCounts {
  let mut within_abs_tolerances = [0; 3];
  let mut within_rel_tolerances = [0; 3];

  for (&a, &e) in izip!(actual, expected) {
    let abs_error = (a - e).abs();

    // Count absolute tolerance compliance
    for (i, &threshold) in thresholds.abs_tolerances.iter().enumerate() {
      if abs_error < threshold {
        within_abs_tolerances[i] += 1;
      }
    }

    // Count relative tolerance compliance (skip if expected is near zero)
    if e.abs() > f64::EPSILON {
      let rel_error = abs_error / e.abs();
      for (i, &threshold) in thresholds.rel_tolerances.iter().enumerate() {
        if rel_error < threshold {
          within_rel_tolerances[i] += 1;
        }
      }
    }
  }

  ToleranceCounts {
    within_abs_tolerances,
    within_rel_tolerances,
  }
}

/// Find location of maximum absolute error for diagnostic purposes
/// Identifies the domain point with worst accuracy for targeted analysis
fn find_max_error_location(x: &Array1<f64>, actual: &Array1<f64>, expected: &Array1<f64>) -> MaxErrorLocation {
  let abs_errors = (actual - expected).mapv(|x| x.abs());

  let max_idx = abs_errors.argmax().unwrap_or(0);

  MaxErrorLocation {
    idx: max_idx,
    x_value: x[max_idx],
  }
}

/// Compute integral (mass) conservation error
/// Formula: $\left|\sum_i y_{\text{num}}(x_i)\Delta x - \sum_i y_{\text{ref}}(x_i)\Delta x\right|$
/// Critical for probability distributions that must integrate to 1
fn compute_mass_error(x: &Array1<f64>, actual: &Array1<f64>, expected: &Array1<f64>) -> f64 {
  // Estimate dx from x values - use mean spacing for non-uniform grids
  let dx = if x.len() > 1 {
    let total_range = x[x.len() - 1] - x[0];
    total_range / (x.len() - 1) as f64
  } else {
    1.0 // fallback for single point
  };

  let actual_mass = actual.sum() * dx;
  let expected_mass = expected.sum() * dx;
  (actual_mass - expected_mass).abs()
}

/// Compute relative L2 norm error for scale-invariant global assessment
/// Formula: $\frac{\left\|y_{\text{num}} - y_{\text{ref}}\right\|_2}{\left\|y_{\text{ref}}\right\|_2}$ where $\left\|v\right\|_2 = \sqrt{\sum_i v_i^2}$
/// Provides normalized measure of total deviation that's independent of signal magnitude
fn compute_relative_l2_norm_error(actual: &Array1<f64>, expected: &Array1<f64>) -> eyre::Result<f64> {
  let diff_norm = (actual - expected).mapv(|x| x.powi(2)).sum().sqrt();
  let expected_norm = expected.mapv(|x| x.powi(2)).sum().sqrt();

  if expected_norm.abs() < f64::EPSILON {
    return make_error!("Expected values have zero L2 norm - cannot compute relative L2 error");
  }

  Ok(diff_norm / expected_norm)
}

/// Compute relative L1 norm error for robust global assessment
/// Formula: $\frac{\sum_i \left|y_{\text{num}}(x_i) - y_{\text{ref}}(x_i)\right|}{\sum_i \left|y_{\text{ref}}(x_i)\right|}$
/// Robust to outliers compared to L2 norm, measures total absolute deviation
fn compute_relative_l1_norm_error(actual: &Array1<f64>, expected: &Array1<f64>) -> eyre::Result<f64> {
  let diff_l1_norm = (actual - expected).mapv(|x| x.abs()).sum();
  let expected_l1_norm = expected.mapv(|x| x.abs()).sum();

  if expected_l1_norm.abs() < f64::EPSILON {
    return make_error!("Expected values have zero L1 norm - cannot compute relative L1 error");
  }

  Ok(diff_l1_norm / expected_l1_norm)
}

/// Compute relative L∞ norm error (supremum error normalized by global peak)
/// Formula: $\frac{\max_i \left|y_{\text{num}}(x_i) - y_{\text{ref}}(x_i)\right|}{\max_i \left|y_{\text{ref}}(x_i)\right|}$
/// Worst-case error normalized by maximum reference value
fn compute_relative_linf_norm_error(actual: &Array1<f64>, expected: &Array1<f64>) -> eyre::Result<f64> {
  let max_abs_error = *(actual - expected).mapv(|x| x.abs()).max().unwrap_or(&0.0);
  let max_expected = *expected.mapv(|x| x.abs()).max().unwrap_or(&0.0);

  if max_expected.abs() < f64::EPSILON {
    return make_error!("Expected values have zero maximum - cannot compute relative L∞ error");
  }

  Ok(max_abs_error / max_expected)
}

/// Compute maximum log error for tail accuracy assessment
/// Formula: $\max_{i: y_{\text{ref}}(x_i) > \tau} \left|\log y_{\text{num}}(x_i) - \log y_{\text{ref}}(x_i)\right|$
/// Sensitive to accuracy for small values, uses threshold τ to avoid log(0)
fn compute_max_log_error(actual: &Array1<f64>, expected: &Array1<f64>, threshold: f64) -> eyre::Result<f64> {
  let mut max_log_error = 0.0_f64;
  let mut valid_points = 0;

  for (&a, &e) in izip!(actual, expected) {
    if e > threshold && a > 0.0 {
      let log_error = (a.ln() - e.ln()).abs();
      max_log_error = max_log_error.max(log_error);
      valid_points += 1;
    }
  }

  if valid_points == 0 {
    return make_error!(
      "No valid points above threshold {:.2e} for log error computation",
      threshold
    );
  }

  Ok(max_log_error)
}

/// Compute symmetry error for symmetric distributions
/// Formula: $\max_i \left|y_{\text{num}}(x_i) - y_{\text{num}}(-x_i)\right|$
/// Measures deviation from symmetry about x=0, useful for Gaussian-like kernels
/// Returns 0.0 if no symmetric pairs exist in the grid
fn compute_symmetry_error(x: &Array1<f64>, actual: &Array1<f64>) -> f64 {
  let mut max_symmetry_error = 0.0_f64;
  let mut found_symmetric_pairs = false;

  for (i, &xi) in x.iter().enumerate() {
    // Find corresponding point at -xi
    let neg_xi = -xi;

    // Find closest point to -xi in the grid
    if let Some((j, _)) = x.iter().enumerate().min_by(|(_, a), (_, b)| {
      (**a - neg_xi)
        .abs()
        .partial_cmp(&(**b - neg_xi).abs())
        .unwrap_or(Ordering::Equal)
    }) {
      // Only consider as symmetric if the points are reasonably close to being symmetric
      let closest_x = x[j];
      if (closest_x - neg_xi).abs() < 1e-10 {
        let symmetry_error = (actual[i] - actual[j]).abs();
        max_symmetry_error = max_symmetry_error.max(symmetry_error);
        found_symmetric_pairs = true;
      }
    }
  }

  // Return 0.0 if no symmetric pairs found (not a symmetric grid)
  if !found_symmetric_pairs {
    0.0
  } else {
    max_symmetry_error
  }
}

/// Compute quantile error (e.g., 95th percentile)
/// Returns error threshold below which the specified quantile of points fall
/// Provides robust measure of error distribution characteristics
fn compute_quantile_error(actual: &Array1<f64>, expected: &Array1<f64>, quantile: f64) -> f64 {
  let mut abs_errors: Vec<f64> = (actual - expected).mapv(|x| x.abs()).to_vec();
  abs_errors.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));

  let index = ((quantile * abs_errors.len() as f64).ceil() as usize).saturating_sub(1);
  abs_errors.get(index).copied().unwrap_or(0.0)
}

/// Compute peak-related accuracy metrics for distribution analysis
/// Analyzes both peak amplitude and location accuracy, critical for probability distributions
/// - Peak value error: $\frac{\left|\max_i y_{\text{num}}(x_i) - \max_i y_{\text{ref}}(x_i)\right|}{\left|\max_i y_{\text{ref}}(x_i)\right|}$
/// - Peak location error: $\left|x_{\operatorname{argmax}(y_{\text{num}})} - x_{\operatorname{argmax}(y_{\text{ref}})}\right|$
fn compute_peak_metrics(x: &Array1<f64>, actual: &Array1<f64>, expected: &Array1<f64>) -> eyre::Result<PeakMetrics> {
  // Find peak values
  let actual_peak = *actual.max().unwrap_or(&0.0);
  let expected_peak = *expected.max().unwrap_or(&0.0);

  // Find peak locations
  let actual_peak_idx = actual.argmax().unwrap_or(0);
  let expected_peak_idx = expected.argmax().unwrap_or(0);

  // Compute peak value error
  if expected_peak.abs() < f64::EPSILON {
    return make_error!("Expected peak value too close to zero: {:.2e}", expected_peak);
  }
  let value_error = (actual_peak - expected_peak).abs() / expected_peak.abs();

  // Compute peak location error
  let location_error = (x[actual_peak_idx] - x[expected_peak_idx]).abs();

  Ok(PeakMetrics {
    value_error,
    location_error,
  })
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_ulps_eq;
  use ndarray::array;

  #[test]
  fn test_perfect_agreement() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let values = array![1.0, 2.0, 3.0, 4.0, 5.0];
    let expected = values.clone();

    let metrics = DomainAgreementMetrics::new(&x, &values, &expected).unwrap();

    assert_ulps_eq!(metrics.abs_error_stats.mean, 0.0);
    assert_ulps_eq!(metrics.abs_error_stats.max, 0.0);
    assert_ulps_eq!(metrics.abs_error_stats.std, 0.0);
    assert_ulps_eq!(metrics.abs_error_stats.bias, 0.0);
    assert_ulps_eq!(metrics.rel_error_stats.mean, 0.0);
    assert_ulps_eq!(metrics.rel_error_stats.max, 0.0);
    assert_ulps_eq!(metrics.rel_error_stats.mape, 0.0);
    assert_ulps_eq!(metrics.rel_error_stats.median, 0.0);
    assert_ulps_eq!(metrics.quality_metrics.rmse, 0.0);
    assert_ulps_eq!(metrics.quality_metrics.r_squared, 1.0);
    assert_ulps_eq!(metrics.quality_metrics.correlation, 1.0);
    assert_ulps_eq!(metrics.quality_metrics.mass_error, 0.0);
    assert_ulps_eq!(metrics.quality_metrics.rel_l2_error, 0.0);
    assert_ulps_eq!(metrics.quality_metrics.rel_l1_error, 0.0);
    assert_ulps_eq!(metrics.quality_metrics.rel_linf_error, 0.0);
    assert_ulps_eq!(metrics.quality_metrics.max_log_error, 0.0);
    assert_ulps_eq!(metrics.quality_metrics.symmetry_error, 0.0);
    assert_ulps_eq!(metrics.quality_metrics.quantile_95_error, 0.0);
    assert_ulps_eq!(metrics.peak_metrics.value_error, 0.0);
    assert_ulps_eq!(metrics.peak_metrics.location_error, 0.0);

    for i in 0..3 {
      assert_ulps_eq!(metrics.abs_tolerance_fraction(i), 1.0);
      assert_ulps_eq!(metrics.rel_tolerance_fraction(i), 1.0);
    }
  }

  #[test]
  fn test_small_differences() {
    let x = array![0.0, 1.0, 2.0];
    let actual = array![1.0, 2.0, 3.0];
    let expected = array![1.001, 2.001, 3.001];

    let metrics = DomainAgreementMetrics::new(&x, &actual, &expected).unwrap();

    assert!(metrics.abs_error_stats.mean < 0.01);
    assert!(metrics.quality_metrics.r_squared > 0.99);
    // Small differences might not reach "Excellent" threshold due to R² calculation
    assert!(matches!(
      metrics.overall_assessment(),
      AgreementAssessment::Excellent | AgreementAssessment::VeryGood
    ));
  }

  #[test]
  fn test_large_differences() {
    let x = array![0.0, 1.0, 2.0];
    let actual = array![1.0, 2.0, 3.0];
    let expected = array![2.0, 4.0, 6.0];

    let metrics = DomainAgreementMetrics::new(&x, &actual, &expected).unwrap();

    assert!(metrics.abs_error_stats.mean > 1.0);
    assert!(metrics.quality_metrics.r_squared < 0.99);
    assert!(matches!(metrics.overall_assessment(), AgreementAssessment::Poor));
  }

  #[test]
  fn test_custom_thresholds() {
    let x = array![0.0, 1.0, 2.0];
    let actual = array![1.0, 2.0, 3.0];
    let expected = array![1.1, 2.1, 3.1];

    let custom_thresholds = ToleranceThresholds {
      abs_tolerances: [0.2, 0.1, 0.05],
      rel_tolerances: [0.2, 0.1, 0.05],
      r2_thresholds: [0.9, 0.8, 0.7],
    };

    let metrics = DomainAgreementMetrics::new_with_thresholds(&x, &actual, &expected, &custom_thresholds).unwrap();

    assert_ulps_eq!(metrics.abs_tolerance_fraction(0), 1.0);
    assert!(metrics.abs_tolerance_fraction(1) < 1.0);
  }

  #[test]
  #[allow(unused_must_use)]
  fn test_error_cases() {
    let x = array![0.0, 1.0];
    let actual = array![1.0, 2.0];
    let expected_wrong_size = array![1.0];

    DomainAgreementMetrics::new(&x, &actual, &expected_wrong_size).unwrap_err();

    let empty_x = array![];
    let empty_actual = array![];
    let empty_expected = array![];

    DomainAgreementMetrics::new(&empty_x, &empty_actual, &empty_expected).unwrap_err();
  }

  #[test]
  fn test_tolerance_fractions() {
    let x = array![0.0, 1.0, 2.0, 3.0];
    let actual = array![1.0, 2.0, 3.0, 4.0];
    let expected = array![1.0, 2.0, 3.0, 4.0];

    let metrics = DomainAgreementMetrics::new(&x, &actual, &expected).unwrap();

    assert_ulps_eq!(metrics.abs_tolerance_fraction(0), 1.0);
    assert_ulps_eq!(metrics.rel_tolerance_fraction(0), 1.0);
    assert_ulps_eq!(metrics.abs_tolerance_fraction(3), 0.0); // Invalid level
  }

  #[test]
  fn test_accessor_methods() {
    let x = array![0.0, 1.0, 2.0];
    let actual = array![1.0, 2.0, 3.0];
    let expected = array![1.0, 2.0, 3.0];

    let metrics = DomainAgreementMetrics::new(&x, &actual, &expected).unwrap();

    assert_ulps_eq!(metrics.abs_tolerance_fraction(0), 1.0);
    assert_ulps_eq!(metrics.rel_tolerance_fraction(0), 1.0);
  }

  #[test]
  fn test_new_metrics() {
    let x = array![0.0, 1.0, 2.0, 3.0];
    let actual = array![1.0, 2.0, 4.0, 5.0]; // Peak at x=3
    let expected = array![1.0, 2.0, 3.0, 6.0]; // Peak at x=3

    let metrics = DomainAgreementMetrics::new(&x, &actual, &expected).unwrap();

    // Test bias (systematic error)
    assert!(metrics.abs_error_stats.bias.abs() < 1.0);

    // Test mass conservation
    assert!(metrics.quality_metrics.mass_error.abs() < 2.0);

    // Test relative L2 error
    assert!(metrics.quality_metrics.rel_l2_error > 0.0);

    // Test peak metrics
    assert!(metrics.peak_metrics.value_error > 0.0); // Different peak values
    assert_ulps_eq!(metrics.peak_metrics.location_error, 0.0); // Same peak location

    // Test median relative error
    assert!(metrics.rel_error_stats.median >= 0.0);
  }

  #[test]
  fn test_peak_location_difference() {
    let x = array![0.0, 1.0, 2.0, 3.0];
    let actual = array![1.0, 5.0, 2.0, 1.0]; // Peak at x=1
    let expected = array![1.0, 2.0, 2.0, 4.0]; // Peak at x=3

    let metrics = DomainAgreementMetrics::new(&x, &actual, &expected).unwrap();

    // Peak locations differ by 2.0 (x=1 vs x=3)
    assert_ulps_eq!(metrics.peak_metrics.location_error, 2.0);
    assert!(metrics.peak_metrics.value_error > 0.0);
  }

  #[test]
  fn test_comprehensive_metrics_example() {
    // Create test data with known characteristics for comprehensive metric validation
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];

    // Reference: Gaussian-like distribution with peak at x=5
    let expected = array![0.1, 0.3, 0.7, 1.2, 1.8, 2.0, 1.7, 1.1, 0.6, 0.2];

    // Actual: Similar but with systematic bias, peak shift, and some noise
    let actual = array![0.15, 0.28, 0.75, 1.15, 1.85, 1.9, 1.75, 1.08, 0.58, 0.25];

    let metrics = DomainAgreementMetrics::new(&x, &actual, &expected).unwrap();

    // Test all error statistics
    assert!(metrics.abs_error_stats.mean > 0.0);
    assert!(metrics.abs_error_stats.max > 0.0);
    assert!(metrics.abs_error_stats.std >= 0.0);
    assert!(metrics.abs_error_stats.bias.abs() < 0.1); // Small systematic bias

    // Test relative error statistics
    assert!(metrics.rel_error_stats.mean.abs() < 0.5); // Reasonable relative error
    assert!(metrics.rel_error_stats.max > 0.0);
    assert!(metrics.rel_error_stats.mape > 0.0);
    assert!(metrics.rel_error_stats.median >= 0.0);

    // Test quality metrics
    assert!(metrics.quality_metrics.rmse > 0.0);
    assert!(metrics.quality_metrics.r_squared > 0.0);
    assert!(metrics.quality_metrics.correlation > 0.0);
    assert!(metrics.quality_metrics.mass_error >= 0.0);
    assert!(metrics.quality_metrics.rel_l2_error > 0.0);
    assert!(metrics.quality_metrics.rel_l1_error > 0.0);
    assert!(metrics.quality_metrics.rel_linf_error > 0.0);
    assert!(metrics.quality_metrics.max_log_error >= 0.0);
    assert!(metrics.quality_metrics.symmetry_error >= 0.0);
    assert!(metrics.quality_metrics.quantile_95_error >= 0.0);

    // Test peak metrics
    assert!(metrics.peak_metrics.value_error >= 0.0);
    assert!(metrics.peak_metrics.location_error >= 0.0);

    // Test tolerance fractions
    assert!(metrics.abs_tolerance_fraction(0) <= 1.0);
    assert!(metrics.rel_tolerance_fraction(0) <= 1.0);

    // Test max error location
    assert!(metrics.max_error_location.idx < metrics.total_points);
    assert!(metrics.max_error_location.x_value >= x[0] && metrics.max_error_location.x_value <= x[x.len() - 1]);
  }

  #[test]
  fn test_new_advanced_metrics() {
    // Test data with known characteristics for specific metric validation
    let x = array![-2.0, -1.0, 0.0, 1.0, 2.0];
    let expected = array![0.1, 0.5, 1.0, 0.5, 0.1]; // Symmetric about x=0
    let actual = array![0.12, 0.48, 0.95, 0.52, 0.11]; // Slight asymmetry and errors

    let metrics = DomainAgreementMetrics::new(&x, &actual, &expected).unwrap();

    // Test L1 norm error
    assert!(metrics.quality_metrics.rel_l1_error > 0.0);
    assert!(metrics.quality_metrics.rel_l1_error < 1.0);

    // Test L∞ norm error
    assert!(metrics.quality_metrics.rel_linf_error > 0.0);
    assert!(metrics.quality_metrics.rel_linf_error < 1.0);

    // Test log error (should be computed for positive values)
    assert!(metrics.quality_metrics.max_log_error >= 0.0);

    // Test symmetry error (should detect asymmetry in actual vs expected data)
    assert!(metrics.quality_metrics.symmetry_error >= 0.0);
    assert!(metrics.quality_metrics.symmetry_error < 1.0);

    // Test quantile error
    assert!(metrics.quality_metrics.quantile_95_error >= 0.0);
    assert!(metrics.quality_metrics.quantile_95_error >= metrics.abs_error_stats.mean);
  }

  #[test]
  fn test_symmetry_specific() {
    // Test symmetric grid with perfect symmetry
    let x = array![-1.0, 0.0, 1.0];
    let actual = array![0.5, 1.0, 0.5]; // Perfectly symmetric
    let expected = array![0.5, 1.0, 0.5];

    let metrics = DomainAgreementMetrics::new(&x, &actual, &expected).unwrap();
    assert_ulps_eq!(metrics.quality_metrics.symmetry_error, 0.0);

    // Test asymmetric data
    let actual_asym = array![0.4, 1.0, 0.6]; // Asymmetric
    let metrics_asym = DomainAgreementMetrics::new(&x, &actual_asym, &expected).unwrap();
    assert!(metrics_asym.quality_metrics.symmetry_error > 0.0);
  }
}
