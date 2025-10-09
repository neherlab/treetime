use crate::distribution::reference::domain_agreement_metrics::DomainAgreementMetrics;
use crate::io::json::{JsonPretty, json_write_file};
use crate::utils::float_fmt::float_to_significant_digits;
use eyre::Report;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};
use std::fs;
use std::marker::PhantomData;
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
pub trait ConvolutionTestRunner<T: TestCase> {
  /// Run a single test case with the specified algorithm
  fn run_test(&self, test_case: &T, algorithm: ConvolutionAlgorithm) -> Result<TestResult<T>, Report>;

  /// Get all test cases for this function type
  fn test_cases(&self) -> &[T];

  /// Get the function type name (e.g., "gaussian", "exponential")
  fn function_type(&self) -> &str;
}

/// Results for a single test case and algorithm combination
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestResult<T: TestCase> {
  pub test_case_name: String,
  pub algorithm: ConvolutionAlgorithm,
  pub test_case: T,
  pub metrics: DomainAgreementMetrics,
  pub execution_time_ms: f64,
  pub grid_points: usize,
  pub peak_error_location: f64,
  pub peak_error_value: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestFailure<T: TestCase> {
  pub test_case_name: String,
  pub algorithm: ConvolutionAlgorithm,
  pub test_case: T,
  pub error: String,
  pub execution_time_ms: f64,
}

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
    let mut outcomes = Vec::new();

    println!(
      "=== {} Convolution Test Framework ===",
      self.runner.function_type().to_uppercase()
    );
    println!(
      "Running {} test cases with {} algorithms ({} total tests)\n",
      self.runner.test_cases().len(),
      self.algorithms.len(),
      self.runner.test_cases().len() * self.algorithms.len()
    );

    for test_case in self.runner.test_cases() {
      let name = test_case.name();
      let description = test_case.description();
      println!("Test Case: {name}");
      println!("  Description: {description}");

      for &algorithm in &self.algorithms {
        print!("  Running {algorithm} algorithm... ");
        let start_time = Instant::now();

        match self.runner.run_test(test_case, algorithm) {
          Ok(result) => {
            let elapsed_ms = result.execution_time_ms;
            let r_squared = result.metrics.quality_metrics.r_squared;
            println!("✓ ({elapsed_ms:.1}ms, R²={r_squared:.6})");
            outcomes.push(TestRunOutcome::Success(result));
          },
          Err(error) => {
            let elapsed_ms = start_time.elapsed().as_secs_f64() * 1000.0;
            println!("✗ Error: {error}");
            outcomes.push(TestRunOutcome::Failure(TestFailure {
              test_case_name: test_case.name().to_owned(),
              algorithm,
              test_case: test_case.clone(),
              error: format!("{error}"),
              execution_time_ms: elapsed_ms,
            }));
          },
        }
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
  pub fn print_summary(&self, summary: &TestSummary) {
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
