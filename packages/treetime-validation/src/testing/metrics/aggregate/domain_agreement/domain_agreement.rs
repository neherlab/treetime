use crate::testing::metrics::aggregate::domain_agreement::error_stats::{
  AbsoluteErrorStats, RelativeErrorStats, compute_absolute_error_statistics, compute_relative_error_statistics,
};
use crate::testing::metrics::aggregate::domain_agreement::peak_metrics::{PeakMetrics, compute_peak_metrics};
use crate::testing::metrics::aggregate::domain_agreement::quality_metrics::{
  QualityMetrics, compute_correlation, compute_mass_error, compute_max_log_error, compute_quantile_error,
  compute_r_squared, compute_relative_l1_norm_error, compute_relative_l2_norm_error, compute_relative_linf_norm_error,
  compute_rmse, compute_symmetry_error,
};
use crate::testing::metrics::aggregate::domain_agreement::tolerance::{
  MaxErrorLocation, compute_tolerance_counts, find_max_error_location,
};
use crate::testing::metrics::config::ToleranceThresholds;
use itertools::Itertools;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::fmt;
use treetime_utils::fmt::float::float_to_digits;
use treetime_utils::make_error;

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
    let rel_error_stats = compute_relative_error_statistics(actual, expected);
    let rmse = compute_rmse(actual, expected);
    let r_squared = compute_r_squared(actual, expected);
    let correlation = compute_correlation(actual, expected);
    let mass_error = compute_mass_error(x, actual, expected);
    let rel_l2_error = compute_relative_l2_norm_error(actual, expected);
    let rel_l1_error = compute_relative_l1_norm_error(actual, expected);
    let rel_linf_error = compute_relative_linf_norm_error(actual, expected);
    let max_log_error = compute_max_log_error(actual, expected, 1e-300);
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
        idx: max_error_idx,
        x_value: max_error_x_value,
      },
      abs_tolerance_fractions,
      rel_tolerance_fractions,
      overall_assessment,
      ..
    } = self;

    let abs_tolerances = ToleranceThresholds::default().abs_tolerances;
    let rel_tolerances = ToleranceThresholds::default().rel_tolerances;

    let mean_rel_error_pct = mean_rel_error * 100.0;
    let max_rel_error_pct = max_rel_error * 100.0;
    let median_rel_error_pct = median_rel_error * 100.0;
    let peak_value_error_pct = peak_value_error * 100.0;

    let abs_tolerance_lines = (0..3)
      .map(|i| {
        let count = (abs_tolerance_fractions[i] * *total_points as f64).round() as usize;
        let formatted_threshold = float_to_digits(abs_tolerances[i], Some(1), None);
        let fraction_percentage = abs_tolerance_fractions[i] * 100.0;
        format!("  < {formatted_threshold:>9}: {count:3}/{total_points} ({fraction_percentage:5.1}%)")
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
