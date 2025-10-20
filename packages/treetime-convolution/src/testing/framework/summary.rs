use crate::testing::algorithms::ConvolutionAlgorithm;
use crate::testing::framework::results::{TestFailure, TestResult};
use crate::testing::framework::test_case::TestCase;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestSummary {
  pub function_type: String,
  pub total_tests: usize,
  pub total_successes: usize,
  pub total_failures: usize,
  pub total_algorithms: usize,
  pub execution_time_total_ms: f64,
  pub algorithm_summaries: Vec<AlgorithmSummary>,
  pub overall_assessment: String,
}

/// Summary for a specific algorithm across all test cases
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlgorithmSummary {
  pub algorithm: ConvolutionAlgorithm,
  pub test_cases_count: usize,
  pub execution_time_total_ms: f64,
  pub execution_time_avg_ms: f64,
  pub r2_min: f64,
  pub r2_max: f64,
  pub r2_mean: f64,
  pub max_abs_error_overall: f64,
  pub max_rel_error_overall: f64,
  pub passed_tests: usize,
  pub failed_tests: usize,
  pub error_failures: usize,
  pub success_rate: f64,
}

impl AlgorithmSummary {
  /// Create a new AlgorithmSummary from test results
  pub fn new<T: TestCase>(
    algorithm: ConvolutionAlgorithm,
    successes: &[&&TestResult<T>],
    failures: &[&&TestFailure<T>],
  ) -> Self {
    let total_runs = successes.len() + failures.len();

    let r2_values: Vec<f64> = successes
      .iter()
      .map(|result| result.metrics.aggregate.domain_agreement.quality_metrics.r_squared)
      .collect_vec();
    let execution_time_total = successes.iter().map(|result| result.execution_time_ms).sum::<f64>()
      + failures.iter().map(|failure| failure.execution_time_ms).sum::<f64>();
    let execution_time_avg = execution_time_total / total_runs as f64;

    let r2_min = r2_values
      .iter()
      .map(|&value| OrderedFloat(value))
      .min()
      .map_or(f64::NAN, |value| value.0);
    let r2_max = r2_values
      .iter()
      .map(|&value| OrderedFloat(value))
      .max()
      .map_or(f64::NAN, |value| value.0);
    let r2_mean = if r2_values.is_empty() {
      f64::NAN
    } else {
      r2_values.iter().sum::<f64>() / r2_values.len() as f64
    };

    let max_abs_error_overall = successes
      .iter()
      .map(|result| OrderedFloat(result.metrics.aggregate.domain_agreement.abs_error_stats.max))
      .max()
      .map_or(f64::NAN, |value| value.0);

    let max_rel_error_overall = successes
      .iter()
      .map(|result| OrderedFloat(result.metrics.aggregate.domain_agreement.rel_error_stats.max))
      .max()
      .map_or(f64::NAN, |value| value.0);

    let metric_failures = successes
      .iter()
      .filter(|result| result.metrics.aggregate.domain_agreement.quality_metrics.r_squared <= 0.95)
      .count();
    let passed_tests = successes.len() - metric_failures;
    let failed_tests = metric_failures + failures.len();
    let success_rate = passed_tests as f64 / total_runs as f64;

    Self {
      algorithm,
      test_cases_count: total_runs,
      execution_time_total_ms: execution_time_total,
      execution_time_avg_ms: execution_time_avg,
      r2_min,
      r2_max,
      r2_mean,
      max_abs_error_overall,
      max_rel_error_overall,
      passed_tests,
      failed_tests,
      error_failures: failures.len(),
      success_rate,
    }
  }
}
