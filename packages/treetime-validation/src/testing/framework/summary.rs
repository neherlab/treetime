use crate::testing::framework::results::{TestFailure, TestResult};
use crate::testing::framework::test_case::TestCase;
use bon::bon;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};

const R2_PASS_THRESHOLD: f64 = 0.95;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestSummary {
  pub(crate) test_suite_name: String,
  pub(crate) total_tests: usize,
  pub(crate) total_successes: usize,
  pub(crate) total_failures: usize,
  pub(crate) total_algorithms: usize,
  pub(crate) execution_time_total_ms: f64,
  pub(crate) algorithm_summaries: Vec<AlgorithmSummary>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlgorithmSummary {
  pub(crate) algorithm_name: String,
  pub(crate) test_cases_count: usize,
  pub(crate) execution_time_total_ms: f64,
  pub(crate) execution_time_avg_ms: f64,
  pub(crate) r2_min: f64,
  pub(crate) r2_max: f64,
  pub(crate) r2_mean: f64,
  pub(crate) max_abs_error_overall: f64,
  pub(crate) max_rel_error_overall: f64,
  pub(crate) passed_tests: usize,
  pub(crate) failed_tests: usize,
  pub(crate) error_failures: usize,
  pub(crate) success_rate: f64,
}

#[bon]
impl AlgorithmSummary {
  #[allow(
    clippy::as_conversions,
    reason = "count/index numeric cast is exact for the domain range"
  )]
  #[builder]
  pub(crate) fn new<T: TestCase>(
    algorithm_name: &str,
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
      .filter(|result| result.metrics.aggregate.domain_agreement.quality_metrics.r_squared <= R2_PASS_THRESHOLD)
      .count();
    let passed_tests = successes.len() - metric_failures;
    let failed_tests = metric_failures + failures.len();
    let success_rate = passed_tests as f64 / total_runs as f64;

    Self {
      algorithm_name: algorithm_name.to_owned(),
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
