use crate::distribution::reference::domain_agreement_metrics::DomainAgreementMetrics;
use crate::distribution::reference::pointwise_metrics::PointwiseMetrics;
use crate::io::json::{JsonPretty, json_write_file, json_write_str};
use crate::io::serde::{array1_as_vec, array1_from_vec};
use crate::utils::float_fmt::float_to_significant_digits;
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

  /// Accuracy metrics comparing `actual_values` vs `expected_values`.
  pub metrics: DomainAgreementMetrics,

  /// Detailed pointwise and spatial accuracy metrics.
  pub pointwise_metrics: PointwiseMetrics,
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
            let r_squared = result.metrics.quality_metrics.r_squared;
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
        .map(|result| result.metrics.quality_metrics.r_squared)
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
        .map(|result| OrderedFloat(result.metrics.abs_error_stats.max))
        .max()
        .map_or(f64::NAN, |value| value.0);

      let max_rel_error_overall = algo_successes
        .iter()
        .map(|result| OrderedFloat(result.metrics.rel_error_stats.max))
        .max()
        .map_or(f64::NAN, |value| value.0);

      let metric_failures = algo_successes
        .iter()
        .filter(|result| result.metrics.quality_metrics.r_squared <= 0.95)
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

    println!("Algorithm Performance:");
    println!(
      "{:>10} {:>8} {:>10} {:>8} {:>8} {:>8} {:>12} {:>12} {:>8} {:>8}",
      "Algorithm", "Tests", "Time(ms)", "R²Min", "R²Max", "R²Mean", "MaxAbsErr", "MaxRelErr%", "Errors", "Success%"
    );
    println!("{}", "-".repeat(110));

    for algo_summary in &summary.algorithm_summaries {
      println!(
        "{:>10} {:>8} {:>10.1} {:>8.4} {:>8.4} {:>8.4} {:>12} {:>12.1} {:>8} {:>8.1}",
        algo_summary.algorithm,
        algo_summary.test_cases_count,
        algo_summary.execution_time_total_ms,
        algo_summary.r2_min,
        algo_summary.r2_max,
        algo_summary.r2_mean,
        float_to_significant_digits(algo_summary.max_abs_error_overall, 3),
        algo_summary.max_rel_error_overall * 100.0,
        algo_summary.error_failures,
        algo_summary.success_rate * 100.0,
      );
    }
    println!();

    self.print_pointwise_summary(outcomes);
  }

  fn print_pointwise_summary(&self, outcomes: &[TestRunOutcome<T>]) {
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

    println!("Pointwise Metrics Across All Tests:");

    let pw_abs_max_overall = successes
      .iter()
      .map(|result| OrderedFloat(result.pointwise_metrics.pointwise_errors.summary.abs_max))
      .max()
      .map_or(0.0, |value| value.0);

    let pw_rel_max_overall = successes
      .iter()
      .map(|result| OrderedFloat(result.pointwise_metrics.pointwise_errors.summary.rel_max))
      .max()
      .map_or(0.0, |value| value.0);

    let pw_d1_max_overall = successes
      .iter()
      .map(|result| OrderedFloat(result.pointwise_metrics.structural_errors.summary.d1_max))
      .max()
      .map_or(0.0, |value| value.0);

    let pw_d2_max_overall = successes
      .iter()
      .map(|result| OrderedFloat(result.pointwise_metrics.structural_errors.summary.d2_max))
      .max()
      .map_or(0.0, |value| value.0);

    let pw_cumulative_max_overall = successes
      .iter()
      .map(|result| OrderedFloat(result.pointwise_metrics.spatial_metrics.summary.cumulative_max_abs))
      .max()
      .map_or(0.0, |value| value.0);

    let pw_sliding_rms_max_overall = successes
      .iter()
      .map(|result| OrderedFloat(result.pointwise_metrics.spatial_metrics.summary.sliding_rms_max))
      .max()
      .map_or(0.0, |value| value.0);

    let pw_tolerance_strict_min_pct = successes
      .iter()
      .map(|result| {
        OrderedFloat(
          result
            .pointwise_metrics
            .error_distributions
            .summary
            .tolerance_pass_fractions[0]
            * 100.0,
        )
      })
      .min()
      .map_or(0.0, |value| value.0);

    println!("  Pointwise Errors:");
    println!(
      "    Max absolute error across all tests: {}",
      float_to_significant_digits(pw_abs_max_overall, 3)
    );
    println!("    Max relative error across all tests: {pw_rel_max_overall:.3e}");

    println!("  Structural Errors:");
    println!(
      "    Max 1st derivative error: {}",
      float_to_significant_digits(pw_d1_max_overall, 3)
    );
    println!(
      "    Max 2nd derivative error: {}",
      float_to_significant_digits(pw_d2_max_overall, 3)
    );

    println!("  Spatial Metrics:");
    println!(
      "    Max cumulative error: {}",
      float_to_significant_digits(pw_cumulative_max_overall, 3)
    );
    println!(
      "    Max sliding RMS error: {}",
      float_to_significant_digits(pw_sliding_rms_max_overall, 3)
    );

    println!("  Error Distributions:");
    println!("    Min strict tolerance pass rate: {pw_tolerance_strict_min_pct:.1}%");

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
