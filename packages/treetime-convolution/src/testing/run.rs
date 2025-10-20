use crate::algos::algos::Algo;
use crate::algos::algos::ConvolutionAlgorithm;
use crate::testing::console::console::ConvolutionTestConsole;
use crate::testing::framework::results::{TestFailure, TestResult, TestRunOutcome};
use crate::testing::framework::summary::AlgorithmSummary;
use crate::testing::framework::summary::TestSummary;
use crate::testing::framework::test_case::TestCase;
use crate::testing::functions::exponential::ExponentialConvInput;
use crate::testing::functions::functions::FunctionType;
use crate::testing::functions::gaussian::GaussianConvInput;
use crate::testing::metrics::metrics::ConvolutionMetrics;
use crate::testing::plots::plots::generate_plot_outputs;
use crate::testing::test_suites::TestSuite;
use clap::Parser;
use eyre::Report;
use itertools::Itertools;
use ndarray::Array1;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
use std::fs;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;
use treetime_io::json::{JsonPretty, json_write_file};
use treetime_utils::make_error;

#[derive(Parser, Clone, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
#[command(
  name = "convolution-test",
  about = "Comprehensive convolution accuracy and performance test framework for all function types",
  version
)]
pub struct Args {
  /// Function types to test
  #[arg(long, default_values_t = FunctionType::all())]
  functions: Vec<FunctionType>,

  /// Output directory for results files
  #[arg(long, default_value = "tmp/testing")]
  output_dir: String,

  /// Algorithms to test
  #[arg(long, default_values_t = ConvolutionAlgorithm::all())]
  algorithms: Vec<ConvolutionAlgorithm>,

  /// Run only specific test cases (comma-separated names, or "all" for all)
  #[arg(long, default_value = "all")]
  test_cases: String,

  /// Enable detailed progress output
  #[arg(long)]
  verbose: bool,

  /// List available test cases for the specified function types
  #[arg(long)]
  list_cases: bool,
}

#[derive(Serialize)]
struct FullResults<'a, T: TestCase> {
  summary: &'a TestSummary,
  outcomes: &'a [TestRunOutcome<T>],
}

pub fn run_convolution_tests() -> Result<(), Report> {
  let args = Args::parse();
  for function_type in &args.functions {
    match function_type {
      FunctionType::Gaussian => {
        run_convolution_tests_impl::<GaussianConvInput>(&args)?;
      },
      FunctionType::Exponential => {
        run_convolution_tests_impl::<ExponentialConvInput>(&args)?;
      },
    }
  }
  Ok(())
}

fn run_convolution_tests_impl<S>(args: &Args) -> Result<(), Report>
where
  S: TestSuite + Default,
{
  let suite = S::default();

  if args.list_cases {
    list_test_cases(&suite);
    return Ok(());
  }

  let output_dir = format!("{}/{}", args.output_dir, suite.function_type());
  let test_cases = filter_test_cases(&suite, Some(args.test_cases.as_str()))?;

  if args.verbose {
    println!("Test Configuration:");
    println!("  Function type: {}", suite.function_type());
    println!("  Algorithms: {:?}", args.algorithms);
    println!("  Test cases: {} selected", test_cases.len());
    println!("  Output directory: {output_dir}\n");
  }

  let outcomes = run_all_tests(&suite, &test_cases, &args.algorithms)?;

  generate_plot_outputs(&output_dir, &outcomes)?;

  let summary = generate_summary(&suite, &outcomes, &args.algorithms);

  ConvolutionTestConsole::print_summary(&summary, &outcomes);

  save_results_json(&output_dir, &outcomes, &summary)?;

  println!(
    "{} convolution test framework completed successfully!",
    suite.function_type()
  );
  println!("Check {output_dir} for detailed results.");

  Ok(())
}

fn list_test_cases<S: TestSuite>(suite: &S) {
  println!("Available {} test cases:", suite.function_type());
  for case in suite.create_test_cases() {
    println!("  - {} : {}", case.name(), case.description());
  }
}

fn filter_test_cases<S: TestSuite>(suite: &S, filter: Option<&str>) -> Result<Vec<S::TestCase>, Report> {
  let filter = filter.and_then(|value| {
    let trimmed = value.trim();
    if trimmed.is_empty() { None } else { Some(trimmed) }
  });

  let all_cases = suite.create_test_cases();
  match filter {
    None | Some("all") => Ok(all_cases),
    Some(filter_str) => {
      let requested_names: BTreeSet<&str> = filter_str.split(',').map(|name| name.trim()).collect();
      let filtered_cases: Vec<S::TestCase> = all_cases
        .iter()
        .filter(|case| requested_names.contains(case.name()))
        .cloned()
        .collect();

      if filtered_cases.is_empty() {
        let available_names = all_cases.iter().map(|case| case.name()).collect::<Vec<_>>().join(", ");
        return make_error!("No matching test cases found for: {filter_str}. Available test cases: {available_names}");
      }

      Ok(filtered_cases)
    },
  }
}

fn run_all_tests<S>(
  suite: &S,
  test_cases: &[S::TestCase],
  algorithms: &[ConvolutionAlgorithm],
) -> Result<Vec<TestRunOutcome<S::TestCase>>, Report>
where
  S: TestSuite,
{
  let total_tests = test_cases.len() * algorithms.len();
  let completed = AtomicUsize::new(0);

  ConvolutionTestConsole::print_header(suite.function_type(), test_cases.len(), algorithms.len());

  ConvolutionTestConsole::print_progress_table_header();

  let test_combinations: Vec<_> = test_cases.iter().cartesian_product(algorithms).collect();

  let outcomes = test_combinations
    .into_par_iter()
    .map(|(test_case, &algorithm)| execute_single_test(suite, test_case, algorithm, &completed, total_tests))
    .collect::<Vec<_>>();

  let failures = collect_failures(&outcomes);

  ConvolutionTestConsole::print_error_summary(&failures)?;

  Ok(outcomes)
}

fn execute_single_test<S>(
  suite: &S,
  test_case: &S::TestCase,
  algorithm: ConvolutionAlgorithm,
  completed: &AtomicUsize,
  total_tests: usize,
) -> TestRunOutcome<S::TestCase>
where
  S: TestSuite,
{
  let start_time = Instant::now();

  let algo = algorithm.instantiate();

  match run_test(suite, test_case, &*algo) {
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
        algorithm: algorithm.to_string(),
        test_case: test_case.clone(),
        error: format!("{error}"),
        execution_time_ms: elapsed_ms,
      })
    },
  }
}

fn run_test<S: TestSuite>(
  suite: &S,
  test_case: &S::TestCase,
  algo: &dyn Algo,
) -> Result<TestResult<S::TestCase>, Report> {
  let start_time = Instant::now();

  let (input_grid_min, input_grid_max) = test_case.input_grid_domain();
  let input_grid_n_points = test_case.input_grid_n_points();
  let input_grid = Array1::linspace(input_grid_min, input_grid_max, input_grid_n_points);

  let (output_grid_min, output_grid_max) = test_case.output_grid_domain();
  let output_grid_n_points = test_case.output_grid_n_points();
  let output_grid = Array1::linspace(output_grid_min, output_grid_max, output_grid_n_points);

  let f_values = suite.create_f(test_case, &input_grid)?;
  let g_values = suite.create_g(test_case, &input_grid)?;

  let actual_values = algo.convolve(&input_grid, &f_values, &g_values, &output_grid)?;
  let expected_values = suite.analytical_convolution(test_case, &output_grid)?;

  let execution_time = start_time.elapsed().as_secs_f64() * 1000.0;

  let metrics = ConvolutionMetrics::new(&output_grid, &actual_values, &expected_values, execution_time)?;

  Ok(TestResult {
    algorithm: algo.name().to_owned(),
    test_case: test_case.clone(),
    execution_time_ms: execution_time,
    f_x_values: input_grid.clone(),
    f_y_values: f_values,
    g_x_values: input_grid,
    g_y_values: g_values,
    evaluation_grid: output_grid,
    actual_values,
    expected_values,
    metrics,
  })
}

fn generate_summary<S>(
  suite: &S,
  outcomes: &[TestRunOutcome<S::TestCase>],
  algorithms: &[ConvolutionAlgorithm],
) -> TestSummary
where
  S: TestSuite,
{
  let successes = collect_successes(outcomes);
  let failures = collect_failures(outcomes);
  let total_execution_time = calculate_total_execution_time(outcomes);
  let algorithm_summaries = build_algorithm_summaries(algorithms, &successes, &failures);
  let overall_assessment = assess_overall_performance(&algorithm_summaries);

  TestSummary {
    function_type: suite.function_type().to_owned(),
    total_tests: outcomes.len(),
    total_successes: successes.len(),
    total_failures: failures.len(),
    total_algorithms: algorithms.len(),
    execution_time_total_ms: total_execution_time,
    algorithm_summaries,
    overall_assessment,
  }
}

fn collect_successes<T: TestCase>(outcomes: &[TestRunOutcome<T>]) -> Vec<&TestResult<T>> {
  outcomes
    .iter()
    .filter_map(|outcome| match outcome {
      TestRunOutcome::Success(result) => Some(result),
      TestRunOutcome::Failure(_) => None,
    })
    .collect()
}

fn collect_failures<T: TestCase>(outcomes: &[TestRunOutcome<T>]) -> Vec<&TestFailure<T>> {
  outcomes
    .iter()
    .filter_map(|outcome| match outcome {
      TestRunOutcome::Failure(failure) => Some(failure),
      TestRunOutcome::Success(_) => None,
    })
    .collect()
}

fn calculate_total_execution_time<T: TestCase>(outcomes: &[TestRunOutcome<T>]) -> f64 {
  outcomes
    .iter()
    .map(|outcome| match outcome {
      TestRunOutcome::Success(result) => result.execution_time_ms,
      TestRunOutcome::Failure(failure) => failure.execution_time_ms,
    })
    .sum()
}

fn build_algorithm_summaries<T: TestCase>(
  algorithms: &[ConvolutionAlgorithm],
  successes: &[&TestResult<T>],
  failures: &[&TestFailure<T>],
) -> Vec<AlgorithmSummary> {
  algorithms
    .iter()
    .filter_map(|&algorithm| {
      let algorithm_name = algorithm.to_string();
      let algo_successes = successes
        .iter()
        .filter(|result| result.algorithm == algorithm_name)
        .collect_vec();
      let algo_failures = failures
        .iter()
        .filter(|failure| failure.algorithm == algorithm_name)
        .collect_vec();
      let total_runs = algo_successes.len() + algo_failures.len();

      if total_runs == 0 {
        return None;
      }

      Some(AlgorithmSummary::new(algorithm, &algo_successes, &algo_failures))
    })
    .collect()
}

fn assess_overall_performance(algorithm_summaries: &[AlgorithmSummary]) -> String {
  if algorithm_summaries.iter().all(|summary| summary.success_rate > 0.9) {
    "Excellent - All algorithms perform well".to_owned()
  } else if algorithm_summaries.iter().any(|summary| summary.success_rate > 0.8) {
    "Good - Some algorithms perform well".to_owned()
  } else {
    "Poor - Significant issues detected".to_owned()
  }
}

fn save_results_json<T>(output_dir: &str, outcomes: &[TestRunOutcome<T>], summary: &TestSummary) -> Result<(), Report>
where
  T: Serialize + TestCase,
{
  fs::create_dir_all(output_dir)?;

  let full_results = FullResults { summary, outcomes };
  let json_path = format!("{output_dir}/convolution_test_results.json");
  json_write_file(&json_path, &full_results, JsonPretty(true))?;
  println!("Saved detailed JSON results to: {json_path}");
  Ok(())
}
