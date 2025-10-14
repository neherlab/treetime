use crate::distribution::reference::domain_agreement_metrics::DomainAgreementMetrics;
use ndarray::Array1;
use serde::{Deserialize, Serialize};

/// Enhanced aggregate metrics combining domain agreement with comprehensive summaries
///
/// These metrics provide domain-wide scalar summaries that give an overall
/// assessment of algorithm performance across all evaluation points.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AggregateMetrics {
  /// Original domain agreement metrics (for backward compatibility)
  pub domain_agreement: DomainAgreementMetrics,
  /// Enhanced performance assessment
  pub performance: PerformanceMetrics,
  /// Quality assessment based on multiple criteria
  pub quality: QualityAssessment,
  /// Algorithm efficiency metrics
  pub efficiency: EfficiencyMetrics,
}

impl AggregateMetrics {
  /// Creates new aggregate metrics from evaluation data
  pub fn new(
    x: &Array1<f64>,
    actual: &Array1<f64>,
    expected: &Array1<f64>,
    execution_time_ms: f64,
  ) -> eyre::Result<Self> {
    let domain_agreement = DomainAgreementMetrics::new(x, actual, expected)?;
    let performance = compute_performance_metrics(actual, expected)?;
    let quality = compute_quality_assessment(&domain_agreement, &performance);
    let efficiency = compute_efficiency_metrics(execution_time_ms, x.len());

    Ok(Self {
      domain_agreement,
      performance,
      quality,
      efficiency,
    })
  }
}

/// Enhanced performance metrics beyond basic domain agreement
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerformanceMetrics {
  /// Signal-to-noise ratio
  pub signal_to_noise_ratio: f64,
  /// Peak signal-to-noise ratio (for image/signal processing)
  pub peak_snr: f64,
  /// Normalized root mean square error
  pub normalized_rmse: f64,
  /// Mean absolute percentage error
  pub mean_absolute_percentage_error: f64,
  /// Symmetric mean absolute percentage error
  pub symmetric_mape: f64,
  /// Coefficient of determination (R²) alternative calculation
  pub coefficient_of_determination: f64,
  /// Nash-Sutcliffe efficiency
  pub nash_sutcliffe_efficiency: f64,
  /// Index of agreement (Willmott's d)
  pub index_of_agreement: f64,
}

/// Comprehensive quality assessment
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityAssessment {
  /// Overall quality score (0-100)
  pub overall_score: f64,
  /// Quality grade (A, B, C, D, F)
  pub grade: QualityGrade,
  /// Individual assessment components
  pub components: QualityComponents,
  /// Detailed verdict with reasoning
  pub verdict: QualityVerdict,
}

/// Quality grade classification
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum QualityGrade {
  /// Excellent quality (90-100)
  A,
  /// Good quality (80-89)
  B,
  /// Acceptable quality (70-79)
  C,
  /// Poor quality (60-69)
  D,
  /// Unacceptable quality (<60)
  F,
}

/// Individual quality assessment components
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityComponents {
  /// Accuracy score (0-100)
  pub accuracy_score: f64,
  /// Precision score (0-100)
  pub precision_score: f64,
  /// Stability score (0-100)
  pub stability_score: f64,
  /// Robustness score (0-100)
  pub robustness_score: f64,
  /// Overall composite score (0-100)
  pub overall_score: f64,
}

/// Detailed quality verdict with reasoning
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityVerdict {
  /// Primary strengths of the algorithm
  pub strengths: Vec<String>,
  /// Primary weaknesses or concerns
  pub weaknesses: Vec<String>,
  /// Recommended use cases
  pub recommendations: Vec<String>,
  /// Critical issues requiring attention
  pub critical_issues: Vec<String>,
}

/// Algorithm efficiency and performance characteristics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EfficiencyMetrics {
  /// Execution time in milliseconds
  pub execution_time_ms: f64,
  /// Time per evaluation point (ms/point)
  pub time_per_point_ms: f64,
  /// Computational efficiency rating
  pub efficiency_rating: EfficiencyRating,
  /// Memory efficiency estimate
  pub memory_efficiency: MemoryEfficiency,
  /// Scalability assessment
  pub scalability: ScalabilityAssessment,
}

/// Computational efficiency classification
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum EfficiencyRating {
  /// Very fast execution
  VeryFast,
  /// Fast execution
  Fast,
  /// Moderate execution time
  Moderate,
  /// Slow execution
  Slow,
  /// Very slow execution
  VerySlow,
}

/// Memory usage efficiency assessment
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MemoryEfficiency {
  /// Estimated memory complexity (e.g., "O(n)", "O(n log n)")
  pub complexity_class: String,
  /// Efficiency rating
  pub rating: MemoryRating,
}

/// Memory efficiency classification
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MemoryRating {
  /// Very memory efficient
  Excellent,
  /// Good memory usage
  Good,
  /// Moderate memory usage
  Moderate,
  /// High memory usage
  Poor,
  /// Excessive memory usage
  Critical,
}

/// Algorithm scalability characteristics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScalabilityAssessment {
  /// Expected scaling behavior with problem size
  pub scaling_behavior: String,
  /// Predicted performance at 10x problem size
  pub performance_10x: f64,
  /// Predicted performance at 100x problem size
  pub performance_100x: f64,
}

// Implementation functions

fn compute_performance_metrics(actual: &Array1<f64>, expected: &Array1<f64>) -> eyre::Result<PerformanceMetrics> {
  let n = actual.len() as f64;
  
  // Basic error calculations
  let errors = actual - expected;
  let squared_errors = errors.mapv(|x| x * x);
  let abs_errors = errors.mapv(|x| x.abs());
  
  let mse = squared_errors.mean().unwrap_or(0.0);
  let rmse = mse.sqrt();
  
  // Signal power and noise power
  let signal_power = expected.mapv(|x| x * x).mean().unwrap_or(0.0);
  let noise_power = mse;
  
  let signal_to_noise_ratio = if noise_power > 0.0 {
    10.0 * (signal_power / noise_power).log10()
  } else {
    f64::INFINITY
  };
  
  // Peak SNR
  let max_expected = expected.iter().copied().fold(0.0_f64, f64::max);
  let peak_snr = if mse > 0.0 && max_expected > 0.0 {
    20.0 * (max_expected / rmse).log10()
  } else {
    f64::INFINITY
  };
  
  // Normalized RMSE
  let expected_range = expected.iter().copied().fold(0.0_f64, f64::max) - 
                      expected.iter().copied().fold(f64::INFINITY, f64::min);
  let normalized_rmse = if expected_range > 0.0 {
    rmse / expected_range
  } else {
    0.0
  };
  
  // MAPE and SMAPE
  let mut mape_sum = 0.0;
  let mut smape_sum = 0.0;
  let mut valid_count = 0;
  
  for i in 0..actual.len() {
    let exp_val = expected[i];
    let act_val = actual[i];
    let abs_error = (act_val - exp_val).abs();
    
    if exp_val.abs() > 1e-10 {
      mape_sum += abs_error / exp_val.abs();
      smape_sum += abs_error / ((exp_val.abs() + act_val.abs()) / 2.0);
      valid_count += 1;
    }
  }
  
  let mean_absolute_percentage_error = if valid_count > 0 {
    (mape_sum / valid_count as f64) * 100.0
  } else {
    0.0
  };
  
  let symmetric_mape = if valid_count > 0 {
    (smape_sum / valid_count as f64) * 100.0
  } else {
    0.0
  };
  
  // R² coefficient of determination
  let expected_mean = expected.mean().unwrap_or(0.0);
  let total_sum_squares = expected.mapv(|x| (x - expected_mean).powi(2)).sum();
  let residual_sum_squares = squared_errors.sum();
  
  let coefficient_of_determination = if total_sum_squares > 0.0 {
    1.0 - (residual_sum_squares / total_sum_squares)
  } else {
    0.0
  };
  
  // Nash-Sutcliffe efficiency (same as R² for this case)
  let nash_sutcliffe_efficiency = coefficient_of_determination;
  
  // Index of agreement (Willmott's d)
  let potential_error_sum = actual.iter()
    .zip(expected.iter())
    .map(|(&a, &e)| (a - expected_mean).abs() + (e - expected_mean).abs())
    .sum::<f64>();
  
  let index_of_agreement = if potential_error_sum > 0.0 {
    1.0 - (residual_sum_squares / potential_error_sum.powi(2))
  } else {
    1.0
  };
  
  Ok(PerformanceMetrics {
    signal_to_noise_ratio,
    peak_snr,
    normalized_rmse,
    mean_absolute_percentage_error,
    symmetric_mape,
    coefficient_of_determination,
    nash_sutcliffe_efficiency,
    index_of_agreement,
  })
}

fn compute_quality_assessment(
  domain_agreement: &DomainAgreementMetrics,
  performance: &PerformanceMetrics,
) -> QualityAssessment {
  // Compute individual component scores
  let accuracy_score = compute_accuracy_score(domain_agreement, performance);
  let precision_score = compute_precision_score(domain_agreement, performance);
  let stability_score = compute_stability_score(domain_agreement);
  let robustness_score = compute_robustness_score(performance);
  
  // Overall score is weighted average
  let overall_score = 0.4 * accuracy_score + 
                     0.3 * precision_score + 
                     0.2 * stability_score + 
                     0.1 * robustness_score;

  let components = QualityComponents {
    accuracy_score,
    precision_score,
    stability_score,
    robustness_score,
    overall_score,
  };
  
  let grade = match overall_score {
    s if s >= 90.0 => QualityGrade::A,
    s if s >= 80.0 => QualityGrade::B,
    s if s >= 70.0 => QualityGrade::C,
    s if s >= 60.0 => QualityGrade::D,
    _ => QualityGrade::F,
  };
  
  let verdict = generate_quality_verdict(domain_agreement, performance, &components);
  
  QualityAssessment {
    overall_score,
    grade,
    components,
    verdict,
  }
}

fn compute_efficiency_metrics(execution_time_ms: f64, num_points: usize) -> EfficiencyMetrics {
  let time_per_point_ms = execution_time_ms / num_points as f64;
  
  let efficiency_rating = match time_per_point_ms {
    t if t < 0.001 => EfficiencyRating::VeryFast,
    t if t < 0.01 => EfficiencyRating::Fast,
    t if t < 0.1 => EfficiencyRating::Moderate,
    t if t < 1.0 => EfficiencyRating::Slow,
    _ => EfficiencyRating::VerySlow,
  };
  
  // Estimate memory complexity (simplified heuristic)
  let memory_efficiency = if num_points < 1000 {
    MemoryEfficiency {
      complexity_class: "O(n)".to_owned(),
      rating: MemoryRating::Excellent,
    }
  } else if num_points < 10000 {
    MemoryEfficiency {
      complexity_class: "O(n log n)".to_owned(),
      rating: MemoryRating::Good,
    }
  } else {
    MemoryEfficiency {
      complexity_class: "O(n²)".to_owned(),
      rating: MemoryRating::Moderate,
    }
  };
  
  // Scalability prediction (simplified)
  let performance_10x = execution_time_ms * 10.0; // Linear scaling assumption
  let performance_100x = execution_time_ms * 100.0;
  
  let scalability = ScalabilityAssessment {
    scaling_behavior: "Linear".to_owned(),
    performance_10x,
    performance_100x,
  };
  
  EfficiencyMetrics {
    execution_time_ms,
    time_per_point_ms,
    efficiency_rating,
    memory_efficiency,
    scalability,
  }
}

// Quality scoring helper functions

fn compute_accuracy_score(domain_agreement: &DomainAgreementMetrics, _performance: &PerformanceMetrics) -> f64 {
  let r2 = domain_agreement.quality_metrics.r_squared;
  (r2 * 100.0).max(0.0).min(100.0)
}

fn compute_precision_score(_domain_agreement: &DomainAgreementMetrics, performance: &PerformanceMetrics) -> f64 {
  let snr = performance.signal_to_noise_ratio;
  if snr.is_infinite() {
    100.0
  } else if snr > 60.0 {
    100.0
  } else if snr > 0.0 {
    (snr / 60.0 * 100.0).max(0.0).min(100.0)
  } else {
    0.0
  }
}

fn compute_stability_score(domain_agreement: &DomainAgreementMetrics) -> f64 {
  let max_error = domain_agreement.abs_error_stats.max;
  let mean_error = domain_agreement.abs_error_stats.mean;
  
  if mean_error > 0.0 {
    let stability_ratio = 1.0 - (max_error / (mean_error * 10.0)).min(1.0);
    (stability_ratio * 100.0).max(0.0)
  } else {
    100.0
  }
}

fn compute_robustness_score(performance: &PerformanceMetrics) -> f64 {
  let mape = performance.mean_absolute_percentage_error;
  if mape < 1.0 {
    100.0
  } else if mape < 100.0 {
    (100.0 - mape).max(0.0)
  } else {
    0.0
  }
}

fn generate_quality_verdict(
  domain_agreement: &DomainAgreementMetrics,
  performance: &PerformanceMetrics,
  components: &QualityComponents,
) -> QualityVerdict {
  let mut strengths = Vec::new();
  let mut weaknesses = Vec::new();
  let mut recommendations = Vec::new();
  let mut critical_issues = Vec::new();
  
  // Analyze strengths
  if components.accuracy_score > 95.0 {
    strengths.push("Excellent accuracy".to_owned());
  }
  if performance.signal_to_noise_ratio > 60.0 {
    strengths.push("Very high signal-to-noise ratio".to_owned());
  }
  if domain_agreement.quality_metrics.r_squared > 0.99 {
    strengths.push("Near-perfect correlation".to_owned());
  }
  
  // Analyze weaknesses
  if components.accuracy_score < 70.0 {
    weaknesses.push("Poor accuracy performance".to_owned());
  }
  if performance.normalized_rmse > 0.1 {
    weaknesses.push("High normalized error".to_owned());
  }
  if domain_agreement.quality_metrics.r_squared < 0.9 {
    weaknesses.push("Low correlation with reference".to_owned());
  }
  
  // Generate recommendations
  if components.overall_score > 90.0 {
    recommendations.push("Suitable for production use".to_owned());
    recommendations.push("Can be used for high-precision applications".to_owned());
  } else if components.overall_score > 70.0 {
    recommendations.push("Suitable for most applications with moderate precision requirements".to_owned());
  } else {
    recommendations.push("Requires improvement before production use".to_owned());
    recommendations.push("Consider parameter tuning or algorithm alternatives".to_owned());
  }
  
  // Identify critical issues
  if domain_agreement.quality_metrics.r_squared < 0.5 {
    critical_issues.push("Very poor correlation - fundamental algorithm issues".to_owned());
  }
  if performance.mean_absolute_percentage_error > 50.0 {
    critical_issues.push("Extremely high percentage errors".to_owned());
  }
  
  QualityVerdict {
    strengths,
    weaknesses,
    recommendations,
    critical_issues,
  }
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
    let result = AggregateMetrics::new(&x, &y, &y, 100.0).unwrap();

    assert_eq!(result.quality.grade, QualityGrade::A);
    assert!(result.performance.signal_to_noise_ratio.is_infinite() || result.performance.signal_to_noise_ratio > 100.0);
    assert_abs_diff_eq!(result.domain_agreement.quality_metrics.r_squared, 1.0, epsilon = 1e-12);
  }

  #[test]
  fn test_poor_agreement() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let expected = array![1.0, 2.0, 3.0, 2.0, 1.0];
    let actual = array![5.0, 6.0, 7.0, 6.0, 5.0]; // Completely different values

    let result = AggregateMetrics::new(&x, &actual, &expected, 100.0).unwrap();

    assert!(matches!(result.quality.grade, QualityGrade::F));
    assert!(result.quality.components.accuracy_score < 50.0);
    assert!(!result.quality.verdict.critical_issues.is_empty());
  }
}