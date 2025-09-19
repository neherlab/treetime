use crate::make_error;
use crate::utils::float_fmt::float_to_digits;
use itertools::{Itertools, izip};
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt;

/// Comprehensive domain-wide agreement metrics between actual and expected solutions
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DomainAgreementMetrics {
  pub total_points: usize,
  pub abs_error_stats: AbsoluteErrorStats,
  pub rel_error_stats: RelativeErrorStats,
  pub quality_metrics: QualityMetrics,
  pub tolerance_counts: ToleranceCounts,
  pub max_error_location: MaxErrorLocation,
  pub thresholds: ToleranceThresholds,
}

impl DomainAgreementMetrics {
  /// Creates new domain agreement metrics from actual and expected values
  pub fn new(x: &Array1<f64>, actual: &Array1<f64>, expected: &Array1<f64>) -> eyre::Result<Self> {
    Self::new_with_thresholds(x, actual, expected, ToleranceThresholds::default())
  }

  /// Creates new domain agreement metrics with custom tolerance thresholds
  pub fn new_with_thresholds(
    x: &Array1<f64>,
    actual: &Array1<f64>,
    expected: &Array1<f64>,
    thresholds: ToleranceThresholds,
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
    let tolerance_counts = compute_tolerance_counts(actual, expected, &thresholds);
    let max_error_location = find_max_error_location(x, actual, expected);

    let quality_metrics = QualityMetrics {
      rmse,
      r_squared,
      correlation,
    };

    Ok(Self {
      total_points,
      abs_error_stats,
      rel_error_stats,
      quality_metrics,
      tolerance_counts,
      max_error_location,
      thresholds,
    })
  }

  /// Provides overall assessment based on R² value
  pub fn overall_assessment(&self) -> AgreementAssessment {
    let r2 = self.quality_metrics.r_squared;
    let thresholds = &self.thresholds.r2_thresholds;

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

  /// Returns percentage of points within absolute tolerance at given level (0-2)
  pub fn abs_tolerance_percentage(&self, level: usize) -> f64 {
    if level >= 3 {
      return 0.0;
    }
    100.0 * self.tolerance_counts.within_abs_tolerances[level] as f64 / self.total_points as f64
  }

  /// Returns percentage of points within relative tolerance at given level (0-2)
  pub fn rel_tolerance_percentage(&self, level: usize) -> f64 {
    if level >= 3 {
      return 0.0;
    }
    100.0 * self.tolerance_counts.within_rel_tolerances[level] as f64 / self.total_points as f64
  }

  /// Returns absolute tolerance value at given level (0-2)
  pub fn abs_tolerance(&self, level: usize) -> f64 {
    self.thresholds.abs_tolerances[level]
  }

  /// Returns relative tolerance value at given level (0-2)
  pub fn rel_tolerance(&self, level: usize) -> f64 {
    self.thresholds.rel_tolerances[level]
  }

  /// Returns count of points within absolute tolerance at given level (0-2)
  pub fn abs_tolerance_count(&self, level: usize) -> usize {
    self.tolerance_counts.within_abs_tolerances[level]
  }

  /// Returns count of points within relative tolerance at given level (0-2)
  pub fn rel_tolerance_count(&self, level: usize) -> usize {
    self.tolerance_counts.within_rel_tolerances[level]
  }
}

/// Assessment levels for domain agreement quality
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

/// Absolute error statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AbsoluteErrorStats {
  pub mean: f64,
  pub max: f64,
  pub std: f64,
}

/// Relative error statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RelativeErrorStats {
  pub mean: f64,
  pub max: f64,
  pub mape: f64,
}

/// Quality metrics for agreement assessment
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityMetrics {
  pub rmse: f64,
  pub r_squared: f64,
  pub correlation: f64,
}

/// Tolerance counts for different threshold levels
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ToleranceCounts {
  pub within_abs_tolerances: [usize; 3],
  pub within_rel_tolerances: [usize; 3],
}

/// Maximum error location information
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MaxErrorLocation {
  pub idx: usize,
  pub x_value: f64,
}

/// Tolerance thresholds for agreement assessment
#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
pub struct ToleranceThresholds {
  #[default(_code = "[1e-6, 1e-9, 1e-12]")]
  pub abs_tolerances: [f64; 3],
  #[default(_code = "[0.01, 0.001, 0.0001]")] // 1%, 0.1%, 0.01%
  pub rel_tolerances: [f64; 3],
  #[default(_code = "[0.999999, 0.9999, 0.99]")] // [excellent, very_good, good]
  pub r2_thresholds: [f64; 3],
}

/// Trait for displaying domain agreement metrics
pub trait DomainAgreementDisplay {
  fn display_metrics(&self) -> String;
  fn display_summary(&self) -> String;
}

impl DomainAgreementDisplay for DomainAgreementMetrics {
  fn display_metrics(&self) -> String {
    let Self {
      total_points,
      abs_error_stats,
      rel_error_stats,
      quality_metrics,
      tolerance_counts,
      max_error_location,
      thresholds,
    } = self;

    let abs_tolerance_lines = (0..3)
      .map(|i| {
        format!(
          "  < {:.0e}: {:3}/{total_points} ({:5.1}%)",
          thresholds.abs_tolerances[i],
          tolerance_counts.within_abs_tolerances[i],
          self.abs_tolerance_percentage(i)
        )
      })
      .join("\n");

    let rel_tolerance_lines = (0..3)
      .map(|i| {
        let percentage = thresholds.rel_tolerances[i] * 100.0;
        let formatted_percentage = float_to_digits(percentage, Some(3), None);
        format!(
          "  < {:>6}%: {:3}/{total_points} ({:5.1}%)",
          formatted_percentage,
          tolerance_counts.within_rel_tolerances[i],
          self.rel_tolerance_percentage(i)
        )
      })
      .join("\n");

    let overall_assessment = self.overall_assessment();

    let mean_abs_error = abs_error_stats.mean;
    let max_abs_error = abs_error_stats.max;
    let std_abs_error = abs_error_stats.std;
    let rmse = quality_metrics.rmse;
    let max_error_x_value = max_error_location.x_value;
    let mean_rel_error = rel_error_stats.mean;
    let mean_rel_error_pct = rel_error_stats.mean * 100.0;
    let max_rel_error = rel_error_stats.max;
    let max_rel_error_pct = rel_error_stats.max * 100.0;
    let mape = rel_error_stats.mape;
    let r_squared = quality_metrics.r_squared;
    let correlation = quality_metrics.correlation;

    format!(
      r#"=== Domain-Wide Agreement Metrics ===
Total evaluation points: {total_points}

Absolute Error Statistics:
  Mean absolute error:    {mean_abs_error:.6e}
  Max absolute error:     {max_abs_error:.6e}
  Std absolute error:     {std_abs_error:.6e}
  RMSE:                   {rmse:.6e}
  Max error at x:         {max_error_x_value:.3}

Relative Error Statistics:
  Mean relative error:    {mean_rel_error:.6e} ({mean_rel_error_pct:.4}%)
  Max relative error:     {max_rel_error:.6e} ({max_rel_error_pct:.4}%)
  MAPE:                   {mape:.6}%

Agreement Quality:
  R² (coefficient of determination): {r_squared:.12}
  Correlation coefficient:           {correlation:.12}

Points within absolute tolerance:
{abs_tolerance_lines}

Points within relative tolerance:
{rel_tolerance_lines}

{overall_assessment}"#
    )
  }

  fn display_summary(&self) -> String {
    format!(
      "R²={:.6}, RMSE={:.2e}, Max err={:.2e}, {}",
      self.quality_metrics.r_squared,
      self.quality_metrics.rmse,
      self.abs_error_stats.max,
      self.overall_assessment()
    )
  }
}

impl fmt::Display for DomainAgreementMetrics {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    write!(f, "{}", self.display_metrics())
  }
}

/// Compute absolute error statistics (mean, max, standard deviation)
fn compute_absolute_error_statistics(actual: &Array1<f64>, expected: &Array1<f64>) -> AbsoluteErrorStats {
  let abs_errors = (actual - expected).mapv(|x| x.abs());

  let mean = abs_errors.mean().unwrap_or(0.0);
  let max = abs_errors.fold(0.0, |acc, &x| acc.max(x));

  let variance = abs_errors.mapv(|x| (x - mean).powi(2)).mean().unwrap_or(0.0);
  let std = variance.sqrt();

  AbsoluteErrorStats { mean, max, std }
}

/// Compute relative error statistics with protection against division by zero
fn compute_relative_error_statistics(actual: &Array1<f64>, expected: &Array1<f64>) -> eyre::Result<RelativeErrorStats> {
  let mut rel_errors = Vec::new();
  let mut abs_percentage_errors = Vec::new();

  for (&a, &e) in izip!(actual, expected) {
    if e.abs() < f64::EPSILON {
      return make_error!("Expected value too close to zero: {:.2e}", e);
    }

    let rel_error = (a - e) / e;
    let abs_percentage_error = rel_error.abs() * 100.0;

    rel_errors.push(rel_error);
    abs_percentage_errors.push(abs_percentage_error);
  }

  let mean = rel_errors.iter().sum::<f64>() / rel_errors.len() as f64;
  let max = rel_errors.iter().fold(0.0, |acc, &x| acc.max(x.abs()));
  let mape = abs_percentage_errors.iter().sum::<f64>() / abs_percentage_errors.len() as f64;

  Ok(RelativeErrorStats { mean, max, mape })
}

/// Compute Root Mean Square Error
fn compute_rmse(actual: &Array1<f64>, expected: &Array1<f64>) -> f64 {
  let squared_errors = (actual - expected).mapv(|x| x.powi(2));
  let mse = squared_errors.mean().unwrap_or(0.0);
  mse.sqrt()
}

/// Compute coefficient of determination (R²)
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
fn compute_correlation(actual: &Array1<f64>, expected: &Array1<f64>) -> eyre::Result<f64> {
  let n = actual.len() as f64;
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

/// Find location of maximum absolute error
fn find_max_error_location(x: &Array1<f64>, actual: &Array1<f64>, expected: &Array1<f64>) -> MaxErrorLocation {
  let abs_errors = (actual - expected).mapv(|x| x.abs());

  let (max_idx, _) = abs_errors
    .indexed_iter()
    .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
    .unwrap_or((0, &0.0));

  MaxErrorLocation {
    idx: max_idx,
    x_value: x[max_idx],
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use ndarray::array;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_perfect_agreement() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let values = array![1.0, 2.0, 3.0, 4.0, 5.0];
    let expected = values.clone();

    let metrics = DomainAgreementMetrics::new(&x, &values, &expected).unwrap();

    assert_eq!(metrics.abs_error_stats.mean, 0.0);
    assert_eq!(metrics.abs_error_stats.max, 0.0);
    assert_eq!(metrics.abs_error_stats.std, 0.0);
    assert_eq!(metrics.rel_error_stats.mean, 0.0);
    assert_eq!(metrics.rel_error_stats.max, 0.0);
    assert_eq!(metrics.rel_error_stats.mape, 0.0);
    assert_eq!(metrics.quality_metrics.rmse, 0.0);
    assert_eq!(metrics.quality_metrics.r_squared, 1.0);
    assert_eq!(metrics.quality_metrics.correlation, 1.0);

    for i in 0..3 {
      assert_eq!(metrics.tolerance_counts.within_abs_tolerances[i], 5);
      assert_eq!(metrics.tolerance_counts.within_rel_tolerances[i], 5);
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
    assert!(matches!(metrics.overall_assessment(), AgreementAssessment::Excellent));
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

    let metrics = DomainAgreementMetrics::new_with_thresholds(&x, &actual, &expected, custom_thresholds).unwrap();

    assert_eq!(metrics.tolerance_counts.within_abs_tolerances[0], 3);
    assert_eq!(metrics.tolerance_counts.within_abs_tolerances[1], 0);
  }

  #[test]
  fn test_error_cases() {
    let x = array![0.0, 1.0];
    let actual = array![1.0, 2.0];
    let expected_wrong_size = array![1.0];

    assert!(DomainAgreementMetrics::new(&x, &actual, &expected_wrong_size).is_err());

    let empty_x = array![];
    let empty_actual = array![];
    let empty_expected = array![];

    assert!(DomainAgreementMetrics::new(&empty_x, &empty_actual, &empty_expected).is_err());
  }

  #[test]
  fn test_tolerance_percentages() {
    let x = array![0.0, 1.0, 2.0, 3.0];
    let actual = array![1.0, 2.0, 3.0, 4.0];
    let expected = array![1.0, 2.0, 3.0, 4.0];

    let metrics = DomainAgreementMetrics::new(&x, &actual, &expected).unwrap();

    assert_eq!(metrics.abs_tolerance_percentage(0), 100.0);
    assert_eq!(metrics.rel_tolerance_percentage(0), 100.0);
    assert_eq!(metrics.abs_tolerance_percentage(3), 0.0); // Invalid level
  }

  #[test]
  fn test_accessor_methods() {
    let x = array![0.0, 1.0, 2.0];
    let actual = array![1.0, 2.0, 3.0];
    let expected = array![1.0, 2.0, 3.0];

    let metrics = DomainAgreementMetrics::new(&x, &actual, &expected).unwrap();

    assert_eq!(metrics.abs_tolerance(0), 1e-6);
    assert_eq!(metrics.rel_tolerance(0), 0.01);
    assert_eq!(metrics.abs_tolerance_count(0), 3);
    assert_eq!(metrics.rel_tolerance_count(0), 3);
  }

  #[test]
  fn test_display_trait() {
    let x = array![0.0, 1.0, 2.0];
    let actual = array![1.0, 2.0, 3.0];
    let expected = array![1.0, 2.0, 3.0];

    let metrics = DomainAgreementMetrics::new(&x, &actual, &expected).unwrap();

    let display_output = metrics.display_metrics();
    assert!(display_output.contains("Domain-Wide Agreement Metrics"));
    assert!(display_output.contains("R²"));

    let summary = metrics.display_summary();
    assert!(summary.contains("R²="));
    assert!(summary.contains("RMSE="));
  }

  #[test]
  fn test_agreement_assessment_display() {
    assert_eq!(
      format!("{}", AgreementAssessment::Excellent),
      "🟢 EXCELLENT: Near-perfect agreement"
    );
    assert_eq!(
      format!("{}", AgreementAssessment::VeryGood),
      "🟡 VERY GOOD: High agreement"
    );
    assert_eq!(
      format!("{}", AgreementAssessment::Good),
      "🟠 GOOD: Reasonable agreement"
    );
    assert_eq!(format!("{}", AgreementAssessment::Poor), "🔴 POOR: Low agreement");
  }
}
