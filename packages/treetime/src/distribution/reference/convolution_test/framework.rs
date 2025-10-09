use crate::distribution::reference::domain_agreement_metrics::DomainAgreementMetrics;
use crate::io::json::{json_write_file, JsonPretty};
use crate::utils::float_fmt::float_to_significant_digits;
use eyre::Report;
use serde::{Deserialize, Serialize};

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

/// Aggregate summary across all tests
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestSummary {
  pub function_type: String,
  pub total_tests: usize,
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
  pub success_rate: f64,
}

/// Generic test framework that works with any function type
pub struct GenericConvolutionTestFramework<T: TestCase, R: ConvolutionTestRunner<T>> {
  pub runner: R,
  pub algorithms: Vec<ConvolutionAlgorithm>,
  pub output_dir: String,
  _phantom: std::marker::PhantomData<T>,
}

impl<T: TestCase, R: ConvolutionTestRunner<T>> GenericConvolutionTestFramework<T, R> {
  /// Create new framework with specific runner
  pub fn new(runner: R, output_dir: String) -> Self {
    Self {
      runner,
      algorithms: ConvolutionAlgorithm::all(),
      output_dir,
      _phantom: std::marker::PhantomData,
    }
  }

  /// Set algorithms to test
  pub fn set_algorithms(&mut self, algorithms: Vec<ConvolutionAlgorithm>) {
    self.algorithms = algorithms;
  }

  /// Run all test cases for all algorithms
  pub fn run_all_tests(&self) -> Result<Vec<TestResult<T>>, Report> {
    let mut results = Vec::new();

    println!("=== {} Convolution Test Framework ===", self.runner.function_type().to_uppercase());
    println!(
      "Running {} test cases with {} algorithms ({} total tests)\n",
      self.runner.test_cases().len(),
      self.algorithms.len(),
      self.runner.test_cases().len() * self.algorithms.len()
    );

    for test_case in self.runner.test_cases() {
      println!("Test Case: {}", test_case.name());
      println!("  Description: {}", test_case.description());

      for &algorithm in &self.algorithms {
        print!("  Running {} algorithm... ", algorithm);
        let start_time = std::time::Instant::now();

        match self.runner.run_test(test_case, algorithm) {
          Ok(result) => {
            let elapsed = start_time.elapsed().as_secs_f64() * 1000.0;
            println!("✓ ({:.1}ms, R²={:.6})", elapsed, result.metrics.quality_metrics.r_squared);
            results.push(result);
          }
          Err(e) => {
            println!("✗ Error: {}", e);
          }
        }
      }
      println!();
    }

    Ok(results)
  }

  /// Generate comprehensive summary from test results
  pub fn generate_summary(&self, results: &[TestResult<T>]) -> TestSummary {
    let total_execution_time: f64 = results.iter().map(|r| r.execution_time_ms).sum();

    let mut algorithm_summaries = Vec::new();

    for &algorithm in &self.algorithms {
      let algo_results: Vec<_> = results.iter().filter(|r| r.algorithm == algorithm).collect();

      if algo_results.is_empty() {
        continue;
      }

      let r2_values: Vec<f64> = algo_results.iter().map(|r| r.metrics.quality_metrics.r_squared).collect();
      let execution_times: Vec<f64> = algo_results.iter().map(|r| r.execution_time_ms).collect();

      let r2_min = r2_values.iter().fold(f64::INFINITY, |a, &b| a.min(b));
      let r2_max = r2_values.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
      let r2_mean = r2_values.iter().sum::<f64>() / r2_values.len() as f64;

      let execution_time_total = execution_times.iter().sum::<f64>();
      let execution_time_avg = execution_time_total / execution_times.len() as f64;

      let max_abs_error_overall = algo_results
        .iter()
        .map(|r| r.metrics.abs_error_stats.max)
        .fold(0.0f64, |a, b| a.max(b));

      let max_rel_error_overall = algo_results
        .iter()
        .map(|r| r.metrics.rel_error_stats.max)
        .fold(0.0f64, |a, b| a.max(b));

      // Simple pass/fail based on R² > 0.95
      let passed_tests = algo_results.iter().filter(|r| r.metrics.quality_metrics.r_squared > 0.95).count();
      let failed_tests = algo_results.len() - passed_tests;
      let success_rate = passed_tests as f64 / algo_results.len() as f64;

      algorithm_summaries.push(AlgorithmSummary {
        algorithm,
        test_cases_count: algo_results.len(),
        execution_time_total_ms: execution_time_total,
        execution_time_avg_ms: execution_time_avg,
        r2_min,
        r2_max,
        r2_mean,
        max_abs_error_overall,
        max_rel_error_overall,
        passed_tests,
        failed_tests,
        success_rate,
      });
    }

    let overall_assessment = if algorithm_summaries.iter().all(|s| s.success_rate > 0.9) {
      "Excellent - All algorithms perform well".to_owned()
    } else if algorithm_summaries.iter().any(|s| s.success_rate > 0.8) {
      "Good - Some algorithms perform well".to_owned()
    } else {
      "Poor - Significant issues detected".to_owned()
    };

    TestSummary {
      function_type: self.runner.function_type().to_owned(),
      total_tests: results.len(),
      total_algorithms: self.algorithms.len(),
      execution_time_total_ms: total_execution_time,
      algorithm_summaries,
      overall_assessment,
    }
  }

  /// Print comprehensive summary to console
  pub fn print_summary(&self, summary: &TestSummary) {
    println!("=== {} CONVOLUTION TEST SUMMARY ===\n", summary.function_type.to_uppercase());

    println!("Overall Statistics:");
    println!("  Function type: {}", summary.function_type);
    println!("  Total tests run: {}", summary.total_tests);
    println!("  Total algorithms: {}", summary.total_algorithms);
    println!("  Total execution time: {:.1}ms", summary.execution_time_total_ms);
    println!("  Assessment: {}\n", summary.overall_assessment);

    println!("Algorithm Performance:");
    println!("{:>10} {:>8} {:>10} {:>8} {:>8} {:>8} {:>12} {:>12} {:>8}",
      "Algorithm", "Tests", "Time(ms)", "R²Min", "R²Max", "R²Mean", "MaxAbsErr", "MaxRelErr%", "Success%");
    println!("{}", "-".repeat(100));

    for algo_summary in &summary.algorithm_summaries {
      println!("{:>10} {:>8} {:>10.1} {:>8.4} {:>8.4} {:>8.4} {:>12} {:>12.1} {:>8.1}",
        algo_summary.algorithm,
        algo_summary.test_cases_count,
        algo_summary.execution_time_total_ms,
        algo_summary.r2_min,
        algo_summary.r2_max,
        algo_summary.r2_mean,
        float_to_significant_digits(algo_summary.max_abs_error_overall, 3),
        algo_summary.max_rel_error_overall * 100.0,
        algo_summary.success_rate * 100.0,
      );
    }
    println!();
  }

  /// Save results to JSON file
  pub fn save_results_json(&self, results: &[TestResult<T>], summary: &TestSummary) -> Result<(), Report> {
    std::fs::create_dir_all(&self.output_dir)?;

    #[derive(Serialize)]
    struct FullResults<'a, T: TestCase> {
      summary: &'a TestSummary,
      results: &'a [TestResult<T>],
    }

    let full_results = FullResults { summary, results };
    let json_path = format!("{}/convolution_test_results.json", self.output_dir);
    json_write_file(&json_path, &full_results, JsonPretty(true))?;
    println!("Saved detailed JSON results to: {}", json_path);
    Ok(())
  }
}
