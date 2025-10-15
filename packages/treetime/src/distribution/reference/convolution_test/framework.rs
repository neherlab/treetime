use crate::distribution::reference::convolution_test::metrics::metrics::ConvolutionMetrics;
use crate::io::json::{JsonPretty, json_write_file, json_write_str};
use crate::io::serde::{array1_as_vec, array1_from_vec};
use crate::utils::float_fmt::float_to_significant_digits;
use crate::utils::iterator::mean_by_key::MeanByKey;
use eyre::Report;
use itertools::Itertools;
use ndarray::Array1;
use ordered_float::OrderedFloat;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::fs;
use std::marker::PhantomData;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use super::algorithms::ConvolutionAlgorithm;

/// Generic trait for test cases that can be used with any function type
pub trait TestCase: Clone + Send + Sync + Serialize {
  /// Get the unique name of this test case
  fn name(&self) -> &str;

  /// Get the description of what this test case tests
  fn description(&self) -> &str;

  /// Get the type of stress this test case applies
  fn stress_type(&self) -> &str;

  /// Get any analytical cautions for this test case
  fn analytical_caution(&self) -> &str;

  /// Get the grid spacing (dx) for this test case
  fn dx(&self) -> f64;
}

/// Generic trait for function-specific test frameworks
pub trait ConvolutionTestRunner<T: TestCase>: Send + Sync {
  /// Run a single test case with the specified algorithm
  fn run_test(&self, test_case: &T, algorithm: ConvolutionAlgorithm) -> Result<TestResult<T>, Report>;

  /// Get all test cases for this function type
  fn test_cases(&self) -> &[T];

  /// Get the function type name (e.g., "gaussian", "exponential")
  fn function_type(&self) -> &str;
}

/// Results for a single test case and algorithm combination.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestResult<T: TestCase> {
  /// Convolution algorithm that was tested.
  pub algorithm: ConvolutionAlgorithm,
  /// Full test case configuration and parameters.
  pub test_case: T,
  /// Wall-clock time taken to execute the test, measured in milliseconds.
  pub execution_time_ms: f64,

  // Input functions
  /// Domain coordinates of the first input function f(x); sample points where `f_y_values` are defined.
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub f_x_values: Array1<f64>,
  /// Function values of f(x) at coordinates in `f_x_values`; the first input waveform to the convolution.
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub f_y_values: Array1<f64>,
  /// Domain coordinates of the second input function g(x); sample points where `g_y_values` are defined.
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub g_x_values: Array1<f64>,
  /// Function values of g(x) at coordinates in `g_x_values`; the second input waveform to the convolution.
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub g_y_values: Array1<f64>,

  // Convolution results
  /// Domain coordinates where the convolution result is sampled; the shared grid used to compare
  /// `actual_values` with `expected_values`.
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub evaluation_grid: Array1<f64>,
  /// Convolution values produced by the tested algorithm, evaluated at points in `evaluation_grid`;
  /// the numerical result under evaluation.
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub actual_values: Array1<f64>,
  /// Analytical or construction-ground-truth convolution values, evaluated at points in `evaluation_grid`;
  /// the reference curve used to compute `metrics`.
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub expected_values: Array1<f64>,

  /// Comprehensive convolution algorithm metrics
  pub metrics: ConvolutionMetrics,
}

/// Results for a test that failed to execute.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestFailure<T: TestCase> {
  /// Convolution algorithm that was tested.
  pub algorithm: ConvolutionAlgorithm,
  /// Full test case configuration and parameters.
  pub test_case: T,
  /// Error message describing the failure.
  pub error: String,
  /// Wall-clock time elapsed before failure, measured in milliseconds.
  pub execution_time_ms: f64,
}

#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum TestRunOutcome<T: TestCase> {
  Success(TestResult<T>),
  Failure(TestFailure<T>),
}

/// Aggregate summary across all tests
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

/// Generic test framework that works with any function type
pub struct GenericConvolutionTestFramework<T: TestCase, R: ConvolutionTestRunner<T>> {
  pub runner: R,
  pub algorithms: Vec<ConvolutionAlgorithm>,
  pub output_dir: String,
  pub _phantom: PhantomData<T>,
}

impl<T: TestCase, R: ConvolutionTestRunner<T>> GenericConvolutionTestFramework<T, R> {
  /// Create new framework with specific runner
  pub fn new(runner: R, output_dir: String) -> Self {
    Self {
      runner,
      algorithms: ConvolutionAlgorithm::all(),
      output_dir,
      _phantom: PhantomData,
    }
  }

  /// Set algorithms to test
  pub fn set_algorithms(&mut self, algorithms: Vec<ConvolutionAlgorithm>) {
    self.algorithms = algorithms;
  }

  /// Run all test cases for all algorithms
  pub fn run_all_tests(&self) -> Result<Vec<TestRunOutcome<T>>, Report> {
    let total_tests = self.runner.test_cases().len() * self.algorithms.len();
    let completed = AtomicUsize::new(0);
    let counter_width = total_tests.to_string().len();

    println!(
      "=== {} Convolution Test Framework ===",
      self.runner.function_type().to_uppercase()
    );
    println!(
      "Running {} test cases with {} algorithms ({} total)\n",
      self.runner.test_cases().len(),
      self.algorithms.len(),
      total_tests
    );

    // Print table header
    println!(
      "| {:^12} | {:^2} | {:^12} | {:^30} | {:^8} | {:^8} |",
      "Progress", "S", "Algorithm", "Test Case", "Time, ms", "R²"
    );
    println!(
      "|{:-<14}|{:-^4}|{:-<14}|{:-<32}|{:->10}|{:->10}|",
      "", "", "", "", "", ""
    );

    // Generate all test-algorithm combinations
    let test_combinations: Vec<_> = self
      .runner
      .test_cases()
      .iter()
      .cartesian_product(&self.algorithms)
      .collect();

    // Run tests in parallel
    let outcomes: Vec<TestRunOutcome<T>> = test_combinations
      .into_par_iter()
      .map(|(test_case, &algorithm)| {
        let start_time = Instant::now();

        match self.runner.run_test(test_case, algorithm) {
          Ok(result) => {
            let completed_count = completed.fetch_add(1, Ordering::Relaxed) + 1;
            let elapsed_ms = result.execution_time_ms;
            let r_squared = result.metrics.aggregate.domain_agreement.quality_metrics.r_squared;
            let progress = format!("[{completed_count}/{total_tests}]");
            println!(
              "| {:<12} | {:^1} | {:<12} | {:<30} | {:>8.1} | {:>8.6} |",
              progress,
              "✅",
              format!("{algorithm}"),
              test_case.name(),
              elapsed_ms,
              r_squared
            );
            TestRunOutcome::Success(result)
          },
          Err(error) => {
            let completed_count = completed.fetch_add(1, Ordering::Relaxed) + 1;
            let elapsed_ms = start_time.elapsed().as_secs_f64() * 1000.0;
            let progress = format!("[{completed_count}/{total_tests}]");
            println!(
              "| {:<12} | {:^1} | {:<12} | {:<30} | {:>8.1} | {:>8} |",
              progress,
              "❌",
              format!("{algorithm}"),
              test_case.name(),
              elapsed_ms,
              "FAILED"
            );
            TestRunOutcome::Failure(TestFailure {
              algorithm,
              test_case: test_case.clone(),
              error: format!("{error}"),
              execution_time_ms: elapsed_ms,
            })
          },
        }
      })
      .collect();

    // Report any errors that occurred
    let failures: Vec<_> = outcomes
      .iter()
      .filter_map(|outcome| match outcome {
        TestRunOutcome::Failure(failure) => Some(failure),
        TestRunOutcome::Success(_) => None,
      })
      .collect();

    if !failures.is_empty() {
      println!("\n\n=== ERRORS ENCOUNTERED ===");
      for failure in failures {
        let test_case_json = json_write_str(&failure.test_case, JsonPretty(true))?;
        println!(
          "\n❌ ERROR: {} + {}: {}\n{test_case_json}",
          failure.test_case.name(),
          failure.algorithm,
          failure.error
        );
      }
      println!();
    }

    Ok(outcomes)
  }

  /// Generate comprehensive summary from test outcomes
  pub fn generate_summary(&self, outcomes: &[TestRunOutcome<T>]) -> TestSummary {
    let successes = outcomes
      .iter()
      .filter_map(|outcome| match outcome {
        TestRunOutcome::Success(result) => Some(result),
        TestRunOutcome::Failure(_) => None,
      })
      .collect_vec();
    let failures = outcomes
      .iter()
      .filter_map(|outcome| match outcome {
        TestRunOutcome::Success(_) => None,
        TestRunOutcome::Failure(failure) => Some(failure),
      })
      .collect_vec();
    let total_execution_time = successes.iter().map(|result| result.execution_time_ms).sum::<f64>()
      + failures.iter().map(|failure| failure.execution_time_ms).sum::<f64>();

    let mut algorithm_summaries = Vec::new();

    for &algorithm in &self.algorithms {
      let algo_successes = successes
        .iter()
        .filter(|result| result.algorithm == algorithm)
        .collect_vec();
      let algo_failures = failures
        .iter()
        .filter(|failure| failure.algorithm == algorithm)
        .collect_vec();
      let total_runs = algo_successes.len() + algo_failures.len();

      if total_runs == 0 {
        continue;
      }

      let r2_values: Vec<f64> = algo_successes
        .iter()
        .map(|result| result.metrics.aggregate.domain_agreement.quality_metrics.r_squared)
        .collect_vec();
      let execution_time_total = algo_successes
        .iter()
        .map(|result| result.execution_time_ms)
        .sum::<f64>()
        + algo_failures
          .iter()
          .map(|failure| failure.execution_time_ms)
          .sum::<f64>();
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

      let max_abs_error_overall = algo_successes
        .iter()
        .map(|result| OrderedFloat(result.metrics.aggregate.domain_agreement.abs_error_stats.max))
        .max()
        .map_or(f64::NAN, |value| value.0);

      let max_rel_error_overall = algo_successes
        .iter()
        .map(|result| OrderedFloat(result.metrics.aggregate.domain_agreement.rel_error_stats.max))
        .max()
        .map_or(f64::NAN, |value| value.0);

      let metric_failures = algo_successes
        .iter()
        .filter(|result| result.metrics.aggregate.domain_agreement.quality_metrics.r_squared <= 0.95)
        .count();
      let passed_tests = algo_successes.len() - metric_failures;
      let failed_tests = metric_failures + algo_failures.len();
      let success_rate = passed_tests as f64 / total_runs as f64;

      algorithm_summaries.push(AlgorithmSummary {
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
        error_failures: algo_failures.len(),
        success_rate,
      });
    }

    let overall_assessment = if algorithm_summaries.iter().all(|summary| summary.success_rate > 0.9) {
      "Excellent - All algorithms perform well".to_owned()
    } else if algorithm_summaries.iter().any(|summary| summary.success_rate > 0.8) {
      "Good - Some algorithms perform well".to_owned()
    } else {
      "Poor - Significant issues detected".to_owned()
    };

    TestSummary {
      function_type: self.runner.function_type().to_owned(),
      total_tests: outcomes.len(),
      total_successes: successes.len(),
      total_failures: failures.len(),
      total_algorithms: self.algorithms.len(),
      execution_time_total_ms: total_execution_time,
      algorithm_summaries,
      overall_assessment,
    }
  }

  /// Print comprehensive summary to console
  pub fn print_summary(&self, summary: &TestSummary, outcomes: &[TestRunOutcome<T>]) {
    println!(
      "=== {} CONVOLUTION TEST SUMMARY ===\n",
      summary.function_type.to_uppercase()
    );

    println!("Overall Statistics:");
    let function_type = &summary.function_type;
    let total_tests = summary.total_tests;
    let total_successes = summary.total_successes;
    let total_failures = summary.total_failures;
    let total_algorithms = summary.total_algorithms;
    let execution_time_total_ms = summary.execution_time_total_ms;
    let overall_assessment = &summary.overall_assessment;
    println!("  Function type: {function_type}");
    println!("  Total tests run: {total_tests}");
    println!("  Successful runs: {total_successes}");
    println!("  Failed runs: {total_failures}");
    println!("  Total algorithms: {total_algorithms}");
    println!("  Total execution time: {execution_time_total_ms:.1}ms");
    println!("  Assessment: {overall_assessment}\n");

    self.print_unified_metrics_table(summary, outcomes);
  }

  fn print_unified_metrics_table(&self, summary: &TestSummary, outcomes: &[TestRunOutcome<T>]) {
    let successes: Vec<_> = outcomes
      .iter()
      .filter_map(|outcome| match outcome {
        TestRunOutcome::Success(result) => Some(result),
        TestRunOutcome::Failure(_) => None,
      })
      .collect();

    if successes.is_empty() {
      return;
    }

    println!("Comprehensive Metrics by Algorithm:");

    let grouped_by_algorithm = successes.into_iter().into_group_map_by(|result| result.algorithm);

    let algorithms: Vec<_> = grouped_by_algorithm.keys().sorted().copied().collect();

    let mut all_metrics = std::collections::HashMap::new();

    for (algorithm, algorithm_results) in &grouped_by_algorithm {
      let algo_summary = summary
        .algorithm_summaries
        .iter()
        .find(|s| s.algorithm == *algorithm)
        .unwrap();

      let correlation_mean: f64 = algorithm_results
        .iter()
        .mean_by_key(|r| r.metrics.aggregate.domain_agreement.quality_metrics.correlation);

      let rmse_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.aggregate.domain_agreement.quality_metrics.rmse))
        .map_or(0.0, |r| r.metrics.aggregate.domain_agreement.quality_metrics.rmse);

      let mass_error_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.aggregate.domain_agreement.quality_metrics.mass_error.abs()))
        .map_or(0.0, |r| {
          r.metrics.aggregate.domain_agreement.quality_metrics.mass_error.abs()
        });

      let rel_l2_error_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.aggregate.domain_agreement.quality_metrics.rel_l2_error))
        .map_or(0.0, |r| {
          r.metrics.aggregate.domain_agreement.quality_metrics.rel_l2_error
        });

      let agg_abs_mean: f64 = algorithm_results
        .iter()
        .mean_by_key(|r| r.metrics.aggregate.domain_agreement.abs_error_stats.mean);

      let agg_rel_mean: f64 = algorithm_results
        .iter()
        .mean_by_key(|r| r.metrics.aggregate.domain_agreement.rel_error_stats.mean);

      let pw_abs_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.pointwise.errors.summary.abs_max))
        .map_or(0.0, |r| r.metrics.pointwise.errors.summary.abs_max);

      let pw_abs_mean = algorithm_results
        .iter()
        .map(|r| r.metrics.pointwise.errors.summary.abs_mean)
        .sum::<f64>()
        / algorithm_results.len() as f64;

      let pw_abs_std = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.pointwise.errors.summary.abs_std))
        .map_or(0.0, |r| r.metrics.pointwise.errors.summary.abs_std);

      let pw_rel_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.pointwise.errors.summary.rel_max))
        .map_or(0.0, |r| r.metrics.pointwise.errors.summary.rel_max);

      let pw_rel_mean = algorithm_results
        .iter()
        .map(|r| r.metrics.pointwise.errors.summary.rel_mean)
        .sum::<f64>()
        / algorithm_results.len() as f64;

      let pw_rel_median = algorithm_results
        .iter()
        .map(|r| r.metrics.pointwise.errors.summary.rel_median)
        .sum::<f64>()
        / algorithm_results.len() as f64;

      let pw_signed_bias_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.pointwise.errors.summary.signed_bias.abs()))
        .map_or(0.0, |r| r.metrics.pointwise.errors.summary.signed_bias.abs());

      let pw_log_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.pointwise.errors.summary.log_max))
        .map_or(0.0, |r| r.metrics.pointwise.errors.summary.log_max);

      let d1_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.pointwise.structural.summary.d1_max))
        .map_or(0.0, |r| r.metrics.pointwise.structural.summary.d1_max);

      let d1_mean = algorithm_results
        .iter()
        .map(|r| r.metrics.pointwise.structural.summary.d1_mean)
        .sum::<f64>()
        / algorithm_results.len() as f64;

      let d2_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.pointwise.structural.summary.d2_max))
        .map_or(0.0, |r| r.metrics.pointwise.structural.summary.d2_max);

      let d2_mean = algorithm_results
        .iter()
        .map(|r| r.metrics.pointwise.structural.summary.d2_mean)
        .sum::<f64>()
        / algorithm_results.len() as f64;

      let symmetry_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.pointwise.structural.summary.symmetry_max))
        .map_or(0.0, |r| r.metrics.pointwise.structural.summary.symmetry_max);

      let monotonicity_violations_max = algorithm_results
        .iter()
        .map(|r| r.metrics.pointwise.structural.summary.monotonicity_violation_count)
        .max()
        .unwrap_or(0);

      let peak_mean = algorithm_results
        .iter()
        .map(|r| r.metrics.spatial.regional.summary.peak_region.mean_error)
        .sum::<f64>()
        / algorithm_results.len() as f64;

      let peak_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.spatial.regional.summary.peak_region.max_error))
        .map_or(0.0, |r| r.metrics.spatial.regional.summary.peak_region.max_error);

      let tail_mean = algorithm_results
        .iter()
        .map(|r| r.metrics.spatial.regional.summary.tail_region.mean_error)
        .sum::<f64>()
        / algorithm_results.len() as f64;

      let tail_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.spatial.regional.summary.tail_region.max_error))
        .map_or(0.0, |r| r.metrics.spatial.regional.summary.tail_region.max_error);

      let cumulative_final = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.spatial.cumulative.summary.final_value.abs()))
        .map_or(0.0, |r| r.metrics.spatial.cumulative.summary.final_value.abs());

      let cumulative_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.spatial.cumulative.summary.max_abs))
        .map_or(0.0, |r| r.metrics.spatial.cumulative.summary.max_abs);

      let sliding_rms_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.spatial.windowed.summary.sliding_rms_max))
        .map_or(0.0, |r| r.metrics.spatial.windowed.summary.sliding_rms_max);

      let sliding_max_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.spatial.windowed.summary.sliding_max_max))
        .map_or(0.0, |r| r.metrics.spatial.windowed.summary.sliding_max_max);

      let tolerance_strict = algorithm_results
        .iter()
        .min_by_key(|r| OrderedFloat(r.metrics.pointwise.tolerance.summary.pass_fractions[0] * 100.0))
        .map_or(0.0, |r| r.metrics.pointwise.tolerance.summary.pass_fractions[0] * 100.0);

      let tolerance_moderate = algorithm_results
        .iter()
        .min_by_key(|r| OrderedFloat(r.metrics.pointwise.tolerance.summary.pass_fractions[1] * 100.0))
        .map_or(0.0, |r| r.metrics.pointwise.tolerance.summary.pass_fractions[1] * 100.0);

      let tolerance_loose = algorithm_results
        .iter()
        .min_by_key(|r| OrderedFloat(r.metrics.pointwise.tolerance.summary.pass_fractions[2] * 100.0))
        .map_or(0.0, |r| r.metrics.pointwise.tolerance.summary.pass_fractions[2] * 100.0);

      let support_mismatch_max = algorithm_results
        .iter()
        .map(|r| r.metrics.pointwise.tolerance.summary.support_mismatch_count)
        .max()
        .unwrap_or(0);

      let metrics = vec![
        ("test_count", format!("{}", algo_summary.test_cases_count)),
        ("exec_time", format!("{:.1}", algo_summary.execution_time_total_ms)),
        ("success_rate", format!("{:.1}%", algo_summary.success_rate * 100.0)),
        ("error_count", format!("{}", algo_summary.error_failures)),
        ("r2_min", format!("{:.4}", algo_summary.r2_min)),
        ("r2_max", format!("{:.4}", algo_summary.r2_max)),
        ("r2_mean", format!("{:.4}", algo_summary.r2_mean)),
        ("correlation_mean", format!("{correlation_mean:.4}")),
        ("rmse_max", float_to_significant_digits(rmse_max, 3)),
        ("mass_error_max", format!("{mass_error_max:.3e}")),
        ("rel_l2_error_max", format!("{rel_l2_error_max:.3e}")),
        (
          "agg_abs_max",
          float_to_significant_digits(algo_summary.max_abs_error_overall, 3),
        ),
        ("agg_abs_mean", float_to_significant_digits(agg_abs_mean, 3)),
        ("agg_rel_max", format!("{:.3e}", algo_summary.max_rel_error_overall)),
        ("agg_rel_mean", format!("{agg_rel_mean:.3e}")),
        ("pw_abs_max", float_to_significant_digits(pw_abs_max, 3)),
        ("pw_abs_mean", float_to_significant_digits(pw_abs_mean, 3)),
        ("pw_abs_std", float_to_significant_digits(pw_abs_std, 3)),
        ("pw_rel_max", format!("{pw_rel_max:.3e}")),
        ("pw_rel_mean", format!("{pw_rel_mean:.3e}")),
        ("pw_rel_median", format!("{pw_rel_median:.3e}")),
        ("pw_signed_bias_max", float_to_significant_digits(pw_signed_bias_max, 3)),
        ("pw_log_max", float_to_significant_digits(pw_log_max, 3)),
        ("d1_max", float_to_significant_digits(d1_max, 3)),
        ("d1_mean", float_to_significant_digits(d1_mean, 3)),
        ("d2_max", float_to_significant_digits(d2_max, 3)),
        ("d2_mean", float_to_significant_digits(d2_mean, 3)),
        ("symmetry_max", float_to_significant_digits(symmetry_max, 3)),
        ("monotonicity_violations", format!("{monotonicity_violations_max}")),
        ("peak_mean", format!("{peak_mean:.3e}")),
        ("peak_max", format!("{peak_max:.3e}")),
        ("tail_mean", format!("{tail_mean:.3e}")),
        ("tail_max", format!("{tail_max:.3e}")),
        ("cumulative_final", float_to_significant_digits(cumulative_final, 3)),
        ("cumulative_max", float_to_significant_digits(cumulative_max, 3)),
        ("sliding_rms_max", float_to_significant_digits(sliding_rms_max, 3)),
        ("sliding_max_max", float_to_significant_digits(sliding_max_max, 3)),
        ("tolerance_strict", format!("{tolerance_strict:.1}%")),
        ("tolerance_moderate", format!("{tolerance_moderate:.1}%")),
        ("tolerance_loose", format!("{tolerance_loose:.1}%")),
        ("support_mismatch", format!("{support_mismatch_max}")),
      ];

      all_metrics.insert(
        *algorithm,
        metrics.into_iter().collect::<std::collections::HashMap<_, _>>(),
      );
    }

    let metric_col_width = "Moderate tolerance pass (min%)".len();
    let mut algo_col_widths = std::collections::HashMap::new();

    for algo in &algorithms {
      let name_w = algo.to_string().len();
      let max_value_w = all_metrics[algo].values().map(|s| s.len()).max().unwrap_or(0);
      algo_col_widths.insert(*algo, name_w.max(max_value_w));
    }

    print!("| {:^width$} |", "Metric", width = metric_col_width);
    for algo in &algorithms {
      print!(" {:^width$} |", algo.to_string(), width = algo_col_widths[algo]);
    }
    println!();

    print!("|{:-^width$}|", "", width = metric_col_width + 2);
    for algo in &algorithms {
      print!("{:-^width$}|", "", width = algo_col_widths[algo] + 2);
    }
    println!();

    let categories = vec![
      (
        "PERFORMANCE",
        vec![
          ("Test count", "test_count"),
          ("Execution time (ms)", "exec_time"),
          ("Success rate", "success_rate"),
          ("Error count", "error_count"),
        ],
      ),
      (
        "QUALITY METRICS",
        vec![
          ("R² min", "r2_min"),
          ("R² max", "r2_max"),
          ("R² mean", "r2_mean"),
          ("Correlation (mean)", "correlation_mean"),
          ("RMSE (max)", "rmse_max"),
          ("Mass conservation error (max)", "mass_error_max"),
          ("Relative L2 error (max)", "rel_l2_error_max"),
        ],
      ),
      (
        "AGGREGATE ERRORS",
        vec![
          ("Absolute max", "agg_abs_max"),
          ("Absolute mean", "agg_abs_mean"),
          ("Relative max", "agg_rel_max"),
          ("Relative mean", "agg_rel_mean"),
        ],
      ),
      (
        "POINTWISE ERRORS",
        vec![
          ("Absolute max", "pw_abs_max"),
          ("Absolute mean", "pw_abs_mean"),
          ("Absolute std dev", "pw_abs_std"),
          ("Relative max", "pw_rel_max"),
          ("Relative mean", "pw_rel_mean"),
          ("Relative median", "pw_rel_median"),
          ("Signed bias (max)", "pw_signed_bias_max"),
          ("Logarithmic max", "pw_log_max"),
        ],
      ),
      (
        "STRUCTURAL",
        vec![
          ("1st derivative error (max)", "d1_max"),
          ("1st derivative error (mean)", "d1_mean"),
          ("2nd derivative error (max)", "d2_max"),
          ("2nd derivative error (mean)", "d2_mean"),
          ("Symmetry residual (max)", "symmetry_max"),
          ("Monotonicity violations (max)", "monotonicity_violations"),
        ],
      ),
      (
        "REGIONAL",
        vec![
          ("Peak region error (mean)", "peak_mean"),
          ("Peak region error (max)", "peak_max"),
          ("Tail region error (mean)", "tail_mean"),
          ("Tail region error (max)", "tail_max"),
        ],
      ),
      (
        "CUMULATIVE & WINDOWED",
        vec![
          ("Cumulative error (final)", "cumulative_final"),
          ("Cumulative error (max abs)", "cumulative_max"),
          ("Sliding RMS error (max)", "sliding_rms_max"),
          ("Sliding max error (max)", "sliding_max_max"),
        ],
      ),
      (
        "TOLERANCE COMPLIANCE",
        vec![
          ("Strict tolerance pass (min%)", "tolerance_strict"),
          ("Moderate tolerance pass (min%)", "tolerance_moderate"),
          ("Loose tolerance pass (min%)", "tolerance_loose"),
          ("Support mismatches (max)", "support_mismatch"),
        ],
      ),
    ];

    for (category, metrics) in categories {
      print!("| {:width$} |", "", width = metric_col_width);
      for algo in &algorithms {
        print!(" {:width$} |", "", width = algo_col_widths[algo]);
      }
      println!();

      print!("| {category:<metric_col_width$} |");
      for algo in &algorithms {
        print!(" {:width$} |", "", width = algo_col_widths[algo]);
      }
      println!();

      for (metric_name, metric_key) in metrics {
        print!("| {metric_name:<metric_col_width$} |");
        for algo in &algorithms {
          let value = &all_metrics[algo][metric_key];
          print!(" {:>width$} |", value, width = algo_col_widths[algo]);
        }
        println!();
      }
    }

    println!();
  }
}

#[derive(Serialize)]
struct FullResults<'a, T: TestCase> {
  summary: &'a TestSummary,
  outcomes: &'a [TestRunOutcome<T>],
}

impl<T: TestCase, R: ConvolutionTestRunner<T>> GenericConvolutionTestFramework<T, R> {
  /// Save results to JSON file
  pub fn save_results_json(&self, outcomes: &[TestRunOutcome<T>], summary: &TestSummary) -> Result<(), Report> {
    fs::create_dir_all(&self.output_dir)?;

    let full_results = FullResults { summary, outcomes };
    let output_dir = &self.output_dir;
    let json_path = format!("{output_dir}/convolution_test_results.json");
    json_write_file(&json_path, &full_results, JsonPretty(true))?;
    println!("Saved detailed JSON results to: {json_path}");
    Ok(())
  }
}
