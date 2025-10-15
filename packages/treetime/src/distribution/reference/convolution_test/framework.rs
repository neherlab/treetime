use crate::distribution::reference::convolution_test::algorithm_summary::AlgorithmSummary;
use crate::distribution::reference::convolution_test::algorithms::ConvolutionAlgorithm;
use crate::distribution::reference::convolution_test::console::ConvolutionTestConsole;
use crate::distribution::reference::convolution_test::metrics::metrics::ConvolutionMetrics;
use crate::io::json::{JsonPretty, json_write_file};
use crate::io::serde::{array1_as_vec, array1_from_vec};
use eyre::Report;
use itertools::Itertools;
use ndarray::Array1;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::fs;
use std::marker::PhantomData;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

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

    ConvolutionTestConsole::print_header(
      self.runner.function_type(),
      self.runner.test_cases().len(),
      self.algorithms.len(),
    );

    ConvolutionTestConsole::print_progress_table_header();

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
            ConvolutionTestConsole::print_success_row(&result, completed_count, total_tests);
            TestRunOutcome::Success(result)
          },
          Err(error) => {
            let completed_count = completed.fetch_add(1, Ordering::Relaxed) + 1;
            let elapsed_ms = start_time.elapsed().as_secs_f64() * 1000.0;
            ConvolutionTestConsole::print_failure_row(test_case, algorithm, elapsed_ms, completed_count, total_tests);
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

    ConvolutionTestConsole::print_error_summary(&failures)?;

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

      algorithm_summaries.push(AlgorithmSummary::new(algorithm, &algo_successes, &algo_failures));
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
    ConvolutionTestConsole::print_summary(summary, outcomes);
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
