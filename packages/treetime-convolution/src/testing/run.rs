use crate::algos::algos::Algo;
use crate::algos::algos::ConvolutionAlgorithm;
use crate::algos::algos::MultiplicationAlgorithm;
use crate::algos::algos::MultiplyAlgo;
use crate::testing::console::console::ConvolutionTestConsole;
use crate::testing::framework::results::{TestFailure, TestResult, TestRunOutcome};
use crate::testing::framework::summary::AlgorithmSummary;
use crate::testing::framework::summary::TestSummary;
use crate::testing::framework::test_case::TestCase;
use crate::testing::framework::tsv_output::generate_tsv_outputs;
use crate::testing::metrics::metrics::ConvolutionMetrics;
use crate::testing::plots::plots::generate_plot_outputs;
use crate::testing::test_suites::test_suites::ChainMultiplicationTestSuite;
use crate::testing::test_suites::test_suites::MultiplicationTestSuite;
use crate::testing::test_suites::test_suites::TestSuite;
use crate::testing::test_suites::test_suites::TestSuiteName;
use clap::Parser;
use eyre::Report;
use itertools::Itertools;
use ndarray::Array1;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
use std::fmt::Display;
use std::fs;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;
use treetime_io::json::{JsonPretty, json_write_file};
use treetime_utils::make_error;

#[derive(Parser, Clone, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
#[command(
  name = "convolution-test",
  about = "Comprehensive convolution and multiplication accuracy test framework for all function types",
  version
)]
pub struct Args {
  /// Test suites to run
  #[arg(long, value_delimiter = ',', default_values_t = TestSuiteName::all())]
  pub test_suites: Vec<TestSuiteName>,

  /// Output directory for results files
  #[arg(long, default_value = "tmp/testing")]
  pub output_dir: String,

  /// Convolution algorithms to test
  #[arg(long, value_delimiter = ',', default_values_t = ConvolutionAlgorithm::all())]
  pub algorithms: Vec<ConvolutionAlgorithm>,

  /// Multiplication algorithms to test
  #[arg(long, value_delimiter = ',', default_values_t = MultiplicationAlgorithm::all())]
  pub mult_algorithms: Vec<MultiplicationAlgorithm>,

  /// Run only specific test cases (comma-separated names, or "all" for all)
  #[arg(long, default_value = "all")]
  pub test_cases: String,

  /// Slowness threshold: only run tests with slowness in [0.0, threshold]
  #[arg(long, default_value_t = 0.5)]
  pub slowness: f64,

  /// Enable detailed progress output
  #[arg(long)]
  pub verbose: bool,

  /// List available test cases for the specified function types
  #[arg(long)]
  pub list_cases: bool,
}

pub fn run_convolution_tests() -> Result<(), Report> {
  let mut args = Args::parse();
  args.test_suites = TestSuiteName::expand(&args.test_suites);
  args.algorithms = ConvolutionAlgorithm::expand(&args.algorithms);
  args.mult_algorithms = MultiplicationAlgorithm::expand(&args.mult_algorithms);
  for suite_name in &args.test_suites {
    suite_name.run_tests(&args)?;
  }
  Ok(())
}

pub fn run_convolution_tests_impl<S>(args: &Args) -> Result<(), Report>
where
  S: TestSuite + Default,
{
  let suite = S::default();

  if args.list_cases {
    list_test_cases(&suite);
    return Ok(());
  }

  let output_dir = format!("{}/{}", args.output_dir, suite.test_suite_name());
  let all_test_cases = suite.create_test_cases();
  let test_cases_after_name_filter = filter_test_cases(&suite, Some(args.test_cases.as_str()))?;
  let selected_test_cases = filter_by_slowness(&test_cases_after_name_filter, args.slowness);

  ConvolutionTestConsole::print_test_configuration(
    suite.test_suite_name(),
    &args.algorithms,
    all_test_cases.len(),
    selected_test_cases.len(),
    args.test_cases != "all",
    args.slowness,
    &output_dir,
  );

  let outcomes = run_all_tests(&suite, &all_test_cases, &selected_test_cases, &args.algorithms)?;

  generate_plot_outputs(&output_dir, &outcomes)?;

  generate_tsv_outputs(&output_dir, &outcomes)?;

  let summary = generate_summary(&suite, &outcomes, &args.algorithms);

  ConvolutionTestConsole::print_summary(&summary, &outcomes);

  save_results_json(&output_dir, &outcomes, &summary)?;

  println!(
    "{} convolution test framework completed successfully!",
    suite.test_suite_name()
  );
  println!("Check {output_dir} for detailed results.\n\n\n");

  Ok(())
}

fn list_test_cases<S: TestSuite>(suite: &S) {
  println!("Available {} test cases:", suite.test_suite_name());
  for case in suite.create_test_cases() {
    println!(
      "  - {} (slowness: {}) : {}",
      case.name(),
      case.slowness(),
      case.description()
    );
  }
}

fn filter_test_cases_by_name<T: TestCase + Clone>(
  all_cases: Vec<T>,
  filter: Option<&str>,
) -> Result<Vec<T>, Report> {
  let filter = filter.and_then(|value| {
    let trimmed = value.trim();
    if trimmed.is_empty() { None } else { Some(trimmed) }
  });

  match filter {
    None | Some("all") => Ok(all_cases),
    Some(filter_str) => {
      let requested_names: BTreeSet<&str> = filter_str.split(',').map(|name| name.trim()).collect();
      let filtered_cases: Vec<T> = all_cases
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

fn filter_test_cases<S: TestSuite>(suite: &S, filter: Option<&str>) -> Result<Vec<S::TestCase>, Report> {
  filter_test_cases_by_name(suite.create_test_cases(), filter)
}

fn filter_by_slowness<T: TestCase>(test_cases: &[T], threshold: f64) -> Vec<T> {
  test_cases
    .iter()
    .filter(|case| case.slowness() <= threshold)
    .cloned()
    .collect()
}

fn run_all_tests<S>(
  suite: &S,
  all_test_cases: &[S::TestCase],
  selected_test_cases: &[S::TestCase],
  algorithms: &[ConvolutionAlgorithm],
) -> Result<Vec<TestRunOutcome<S::TestCase>>, Report>
where
  S: TestSuite,
{
  let total_tests = selected_test_cases.len() * algorithms.len();
  let completed = AtomicUsize::new(0);

  ConvolutionTestConsole::print_header(suite.test_suite_name(), selected_test_cases.len(), algorithms.len());

  ConvolutionTestConsole::print_progress_table_header();

  let selected_names: BTreeSet<_> = selected_test_cases.iter().map(|tc| tc.name()).collect();

  let test_combinations: Vec<_> = all_test_cases
    .iter()
    .flat_map(|test_case| {
      let is_selected = selected_names.contains(test_case.name());
      algorithms
        .iter()
        .map(move |algorithm| (test_case, *algorithm, is_selected))
    })
    .collect();

  let outcomes = test_combinations
    .into_par_iter()
    .filter_map(|(test_case, algorithm, is_selected)| {
      if is_selected {
        Some(execute_single_test(
          suite,
          test_case,
          algorithm,
          &completed,
          total_tests,
        ))
      } else {
        ConvolutionTestConsole::print_skipped_row(test_case, algorithm);
        None
      }
    })
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
        error: format!("{error:?}"),
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

  let (evaluation_grid_min, evaluation_grid_max) = (input_grid_min * 2.0, input_grid_max * 2.0);
  let evaluation_grid_n_points = 2 * input_grid_n_points - 1;
  let evaluation_grid = Array1::linspace(evaluation_grid_min, evaluation_grid_max, evaluation_grid_n_points);

  let f_values = suite.create_f(test_case, &input_grid)?;
  let g_values = suite.create_g(test_case, &input_grid)?;

  let dx = input_grid[1] - input_grid[0];
  let actual_values = algo.convolve(dx, &f_values, &g_values)?;
  let expected_values = suite.analytical_convolution(test_case, &evaluation_grid)?;

  let execution_time = start_time.elapsed().as_secs_f64() * 1000.0;

  let metrics = ConvolutionMetrics::new(&evaluation_grid, &actual_values, &expected_values, execution_time)?;

  Ok(TestResult {
    algorithm: algo.name().to_owned(),
    test_case: test_case.clone(),
    execution_time_ms: execution_time,
    f_x_values: input_grid.clone(),
    f_y_values: f_values,
    g_x_values: input_grid,
    g_y_values: g_values,
    evaluation_grid,
    actual_values,
    expected_values,
    metrics,
    log_scale_actual: None,
    log_scale_expected: None,
    log_scale_error: None,
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
    test_suite_name: suite.test_suite_name().to_owned(),
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

#[derive(Serialize)]
struct ResultsJson<'a, T: TestCase> {
  summary: &'a TestSummary,
  outcomes: &'a [TestRunOutcome<T>],
}

fn save_results_json<T>(output_dir: &str, outcomes: &[TestRunOutcome<T>], summary: &TestSummary) -> Result<(), Report>
where
  T: Serialize + TestCase,
{
  fs::create_dir_all(output_dir)?;

  let results = ResultsJson { summary, outcomes };
  let json_path = format!("{output_dir}/{}_results.json", summary.test_suite_name);
  json_write_file(&json_path, &results, JsonPretty(true))?;
  println!("Saved detailed JSON results to: {json_path}");
  Ok(())
}

pub fn run_multiplication_tests_impl<S>(args: &Args) -> Result<(), Report>
where
  S: MultiplicationTestSuite + Default,
{
  let suite = S::default();

  if args.list_cases {
    list_multiplication_test_cases(&suite);
    return Ok(());
  }

  let output_dir = format!("{}/{}", args.output_dir, suite.test_suite_name());
  let all_test_cases = suite.create_test_cases();
  let test_cases_after_name_filter = filter_multiplication_test_cases(&suite, Some(args.test_cases.as_str()))?;
  let selected_test_cases = filter_by_slowness(&test_cases_after_name_filter, args.slowness);

  print_multiplication_test_configuration(
    suite.test_suite_name(),
    &args.mult_algorithms,
    all_test_cases.len(),
    selected_test_cases.len(),
    args.test_cases != "all",
    args.slowness,
    &output_dir,
  );

  let outcomes = run_all_multiplication_tests(&suite, &all_test_cases, &selected_test_cases, &args.mult_algorithms)?;

  generate_plot_outputs(&output_dir, &outcomes)?;

  generate_tsv_outputs(&output_dir, &outcomes)?;

  let summary = generate_multiplication_summary(&suite, &outcomes, &args.mult_algorithms);

  ConvolutionTestConsole::print_summary(&summary, &outcomes);

  save_results_json(&output_dir, &outcomes, &summary)?;

  println!(
    "{} multiplication test framework completed successfully!",
    suite.test_suite_name()
  );
  println!("Check {output_dir} for detailed results.\n\n\n");

  Ok(())
}

fn list_multiplication_test_cases<S: MultiplicationTestSuite>(suite: &S) {
  println!("Available {} test cases:", suite.test_suite_name());
  for case in suite.create_test_cases() {
    println!(
      "  - {} (slowness: {}) : {}",
      case.name(),
      case.slowness(),
      case.description()
    );
  }
}

fn filter_multiplication_test_cases<S: MultiplicationTestSuite>(
  suite: &S,
  filter: Option<&str>,
) -> Result<Vec<S::TestCase>, Report> {
  filter_test_cases_by_name(suite.create_test_cases(), filter)
}

fn print_multiplication_test_configuration<A: Display>(
  test_suite_name: &str,
  algorithms: &[A],
  total_cases: usize,
  selected_cases: usize,
  has_filter: bool,
  slowness: f64,
  output_dir: &str,
) {
  println!("\n=== {test_suite_name} Multiplication Test Configuration ===");
  println!(
    "Algorithms: {}",
    algorithms.iter().map(|a| a.to_string()).collect::<Vec<_>>().join(", ")
  );
  println!("Test cases: {selected_cases}/{total_cases}");
  if has_filter {
    println!("  (filtered by name)");
  }
  println!("Slowness threshold: {slowness}");
  println!("Output directory: {output_dir}\n");
}

fn run_all_multiplication_tests<S>(
  suite: &S,
  all_test_cases: &[S::TestCase],
  selected_test_cases: &[S::TestCase],
  algorithms: &[MultiplicationAlgorithm],
) -> Result<Vec<TestRunOutcome<S::TestCase>>, Report>
where
  S: MultiplicationTestSuite,
{
  let total_tests = selected_test_cases.len() * algorithms.len();
  let completed = AtomicUsize::new(0);

  println!(
    "\nRunning {} multiplication tests ({} cases x {} algorithms)...\n",
    suite.test_suite_name(),
    selected_test_cases.len(),
    algorithms.len()
  );

  println!(
    "{:<30} {:<25} {:>10} {:>15}",
    "Test Case", "Algorithm", "Time (ms)", "Status"
  );
  println!("{:-<80}", "");

  let selected_names: BTreeSet<_> = selected_test_cases.iter().map(|tc| tc.name()).collect();

  let test_combinations: Vec<_> = all_test_cases
    .iter()
    .flat_map(|test_case| {
      let is_selected = selected_names.contains(test_case.name());
      algorithms
        .iter()
        .map(move |algorithm| (test_case, *algorithm, is_selected))
    })
    .collect();

  let outcomes = test_combinations
    .into_par_iter()
    .filter_map(|(test_case, algorithm, is_selected)| {
      if is_selected {
        Some(execute_single_multiplication_test(
          suite,
          test_case,
          algorithm,
          &completed,
          total_tests,
        ))
      } else {
        println!(
          "{:<30} {:<25} {:>10} {:>15}",
          test_case.name(),
          algorithm,
          "-",
          "SKIPPED"
        );
        None
      }
    })
    .collect::<Vec<_>>();

  let failures = collect_failures(&outcomes);

  if !failures.is_empty() {
    println!("\n=== Errors ===");
    for failure in &failures {
      println!(
        "  {} ({}): {}",
        failure.test_case.name(),
        failure.algorithm,
        failure.error
      );
    }
  }

  Ok(outcomes)
}

fn execute_single_multiplication_test<S>(
  suite: &S,
  test_case: &S::TestCase,
  algorithm: MultiplicationAlgorithm,
  completed: &AtomicUsize,
  total_tests: usize,
) -> TestRunOutcome<S::TestCase>
where
  S: MultiplicationTestSuite,
{
  let start_time = Instant::now();

  let algo = algorithm.instantiate();

  match run_multiplication_test(suite, test_case, &*algo) {
    Ok(result) => {
      let completed_count = completed.fetch_add(1, Ordering::Relaxed) + 1;
      println!(
        "{:<30} {:<25} {:>10.2} {:>15} [{}/{}]",
        test_case.name(),
        algorithm,
        result.execution_time_ms,
        "OK",
        completed_count,
        total_tests
      );
      TestRunOutcome::Success(result)
    },
    Err(error) => {
      let completed_count = completed.fetch_add(1, Ordering::Relaxed) + 1;
      let elapsed_ms = start_time.elapsed().as_secs_f64() * 1000.0;
      println!(
        "{:<30} {:<25} {:>10.2} {:>15} [{}/{}]",
        test_case.name(),
        algorithm,
        elapsed_ms,
        "FAIL",
        completed_count,
        total_tests
      );
      TestRunOutcome::Failure(TestFailure {
        algorithm: algorithm.to_string(),
        test_case: test_case.clone(),
        error: format!("{error:?}"),
        execution_time_ms: elapsed_ms,
      })
    },
  }
}

fn run_multiplication_test<S: MultiplicationTestSuite>(
  suite: &S,
  test_case: &S::TestCase,
  algo: &dyn MultiplyAlgo,
) -> Result<TestResult<S::TestCase>, Report> {
  let start_time = Instant::now();

  let (input_grid_min, input_grid_max) = test_case.input_grid_domain();
  let input_grid_n_points = test_case.input_grid_n_points();
  let input_grid = Array1::linspace(input_grid_min, input_grid_max, input_grid_n_points);

  let f_values = suite.create_f(test_case, &input_grid)?;
  let g_values = suite.create_g(test_case, &input_grid)?;

  let actual_values = algo.multiply(&f_values, &g_values);
  let (expected_shape, expected_log_scale) = suite.analytical_multiplication(test_case, &input_grid)?;

  let expected_values = expected_shape.mapv(|v| v * expected_log_scale.exp());

  let execution_time = start_time.elapsed().as_secs_f64() * 1000.0;

  let metrics = ConvolutionMetrics::new(&input_grid, &actual_values, &expected_values, execution_time)?;

  Ok(TestResult {
    algorithm: algo.name().to_owned(),
    test_case: test_case.clone(),
    execution_time_ms: execution_time,
    f_x_values: input_grid.clone(),
    f_y_values: f_values,
    g_x_values: input_grid.clone(),
    g_y_values: g_values,
    evaluation_grid: input_grid,
    actual_values,
    expected_values,
    metrics,
    log_scale_actual: None,
    log_scale_expected: None,
    log_scale_error: None,
  })
}

fn generate_multiplication_summary_from_outcomes<T: TestCase>(
  test_suite_name: &str,
  outcomes: &[TestRunOutcome<T>],
  algorithms: &[MultiplicationAlgorithm],
) -> TestSummary {
  let successes = collect_successes(outcomes);
  let failures = collect_failures(outcomes);
  let total_execution_time = calculate_total_execution_time(outcomes);
  let algorithm_summaries = build_multiplication_algorithm_summaries(algorithms, &successes, &failures);
  let overall_assessment = assess_overall_performance(&algorithm_summaries);

  TestSummary {
    test_suite_name: test_suite_name.to_owned(),
    total_tests: outcomes.len(),
    total_successes: successes.len(),
    total_failures: failures.len(),
    total_algorithms: algorithms.len(),
    execution_time_total_ms: total_execution_time,
    algorithm_summaries,
    overall_assessment,
  }
}

fn generate_multiplication_summary<S>(
  suite: &S,
  outcomes: &[TestRunOutcome<S::TestCase>],
  algorithms: &[MultiplicationAlgorithm],
) -> TestSummary
where
  S: MultiplicationTestSuite,
{
  generate_multiplication_summary_from_outcomes(suite.test_suite_name(), outcomes, algorithms)
}

fn build_multiplication_algorithm_summaries<T: TestCase>(
  algorithms: &[MultiplicationAlgorithm],
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

      Some(AlgorithmSummary::new_from_name(
        &algorithm_name,
        &algo_successes,
        &algo_failures,
      ))
    })
    .collect()
}

pub fn run_chain_multiplication_tests_impl<S>(args: &Args) -> Result<(), Report>
where
  S: ChainMultiplicationTestSuite + Default,
{
  let suite = S::default();

  if args.list_cases {
    list_chain_multiplication_test_cases(&suite);
    return Ok(());
  }

  let output_dir = format!("{}/{}", args.output_dir, suite.test_suite_name());
  let all_test_cases = suite.create_test_cases();
  let test_cases_after_name_filter = filter_chain_multiplication_test_cases(&suite, Some(args.test_cases.as_str()))?;
  let selected_test_cases = filter_by_slowness(&test_cases_after_name_filter, args.slowness);

  print_multiplication_test_configuration(
    suite.test_suite_name(),
    &args.mult_algorithms,
    all_test_cases.len(),
    selected_test_cases.len(),
    args.test_cases != "all",
    args.slowness,
    &output_dir,
  );

  let outcomes =
    run_all_chain_multiplication_tests(&suite, &all_test_cases, &selected_test_cases, &args.mult_algorithms)?;

  generate_plot_outputs(&output_dir, &outcomes)?;

  generate_tsv_outputs(&output_dir, &outcomes)?;

  let summary = generate_chain_multiplication_summary(&suite, &outcomes, &args.mult_algorithms);

  ConvolutionTestConsole::print_summary(&summary, &outcomes);

  save_results_json(&output_dir, &outcomes, &summary)?;

  println!(
    "{} chain multiplication test framework completed successfully!",
    suite.test_suite_name()
  );
  println!("Check {output_dir} for detailed results.\n\n\n");

  Ok(())
}

fn list_chain_multiplication_test_cases<S: ChainMultiplicationTestSuite>(suite: &S) {
  println!("Available {} test cases:", suite.test_suite_name());
  for case in suite.create_test_cases() {
    println!(
      "  - {} (slowness: {}) : {}",
      case.name(),
      case.slowness(),
      case.description()
    );
  }
}

fn filter_chain_multiplication_test_cases<S: ChainMultiplicationTestSuite>(
  suite: &S,
  filter: Option<&str>,
) -> Result<Vec<S::TestCase>, Report> {
  filter_test_cases_by_name(suite.create_test_cases(), filter)
}

fn run_all_chain_multiplication_tests<S>(
  suite: &S,
  all_test_cases: &[S::TestCase],
  selected_test_cases: &[S::TestCase],
  algorithms: &[MultiplicationAlgorithm],
) -> Result<Vec<TestRunOutcome<S::TestCase>>, Report>
where
  S: ChainMultiplicationTestSuite,
{
  let total_tests = selected_test_cases.len() * algorithms.len();
  let completed = AtomicUsize::new(0);

  println!(
    "\nRunning {} chain multiplication tests ({} cases x {} algorithms)...\n",
    suite.test_suite_name(),
    selected_test_cases.len(),
    algorithms.len()
  );

  println!(
    "{:<30} {:<25} {:>10} {:>15}",
    "Test Case", "Algorithm", "Time (ms)", "Status"
  );
  println!("{:-<80}", "");

  let selected_names: BTreeSet<_> = selected_test_cases.iter().map(|tc| tc.name()).collect();

  let test_combinations: Vec<_> = all_test_cases
    .iter()
    .flat_map(|test_case| {
      let is_selected = selected_names.contains(test_case.name());
      algorithms
        .iter()
        .map(move |algorithm| (test_case, *algorithm, is_selected))
    })
    .collect();

  let outcomes = test_combinations
    .into_par_iter()
    .filter_map(|(test_case, algorithm, is_selected)| {
      if is_selected {
        Some(execute_single_chain_multiplication_test(
          suite,
          test_case,
          algorithm,
          &completed,
          total_tests,
        ))
      } else {
        println!(
          "{:<30} {:<25} {:>10} {:>15}",
          test_case.name(),
          algorithm,
          "-",
          "SKIPPED"
        );
        None
      }
    })
    .collect::<Vec<_>>();

  let failures = collect_failures(&outcomes);

  if !failures.is_empty() {
    println!("\n=== Errors ===");
    for failure in &failures {
      println!(
        "  {} ({}): {}",
        failure.test_case.name(),
        failure.algorithm,
        failure.error
      );
    }
  }

  Ok(outcomes)
}

fn execute_single_chain_multiplication_test<S>(
  suite: &S,
  test_case: &S::TestCase,
  algorithm: MultiplicationAlgorithm,
  completed: &AtomicUsize,
  total_tests: usize,
) -> TestRunOutcome<S::TestCase>
where
  S: ChainMultiplicationTestSuite,
{
  let start_time = Instant::now();

  let algo = algorithm.instantiate();

  match run_chain_multiplication_test(suite, test_case, &*algo) {
    Ok(result) => {
      let completed_count = completed.fetch_add(1, Ordering::Relaxed) + 1;
      println!(
        "{:<30} {:<25} {:>10.2} {:>15} [{}/{}]",
        test_case.name(),
        algorithm,
        result.execution_time_ms,
        "OK",
        completed_count,
        total_tests
      );
      TestRunOutcome::Success(result)
    },
    Err(error) => {
      let completed_count = completed.fetch_add(1, Ordering::Relaxed) + 1;
      let elapsed_ms = start_time.elapsed().as_secs_f64() * 1000.0;
      println!(
        "{:<30} {:<25} {:>10.2} {:>15} [{}/{}]",
        test_case.name(),
        algorithm,
        elapsed_ms,
        "FAIL",
        completed_count,
        total_tests
      );
      TestRunOutcome::Failure(TestFailure {
        algorithm: algorithm.to_string(),
        test_case: test_case.clone(),
        error: format!("{error:?}"),
        execution_time_ms: elapsed_ms,
      })
    },
  }
}

fn run_chain_multiplication_test<S: ChainMultiplicationTestSuite>(
  suite: &S,
  test_case: &S::TestCase,
  algo: &dyn MultiplyAlgo,
) -> Result<TestResult<S::TestCase>, Report> {
  let start_time = Instant::now();

  let (input_grid_min, input_grid_max) = test_case.input_grid_domain();
  let input_grid_n_points = test_case.input_grid_n_points();
  let input_grid = Array1::linspace(input_grid_min, input_grid_max, input_grid_n_points);

  let factors = suite.create_factors(test_case, &input_grid)?;
  let factor_refs: Vec<&Array1<f64>> = factors.iter().collect();

  let (actual_shape, actual_log_scale) = if factors.is_empty() {
    (Array1::ones(input_grid.len()), 0.0)
  } else {
    algo.multiply_many(&factor_refs)
  };
  let (expected_shape, expected_log_scale) = suite.analytical_chain_multiplication(test_case, &input_grid)?;

  let log_scale_error = (actual_log_scale - expected_log_scale).abs();

  let actual_values = actual_shape.mapv(|v| v * actual_log_scale.exp());
  let expected_values = expected_shape.mapv(|v| v * expected_log_scale.exp());

  let execution_time = start_time.elapsed().as_secs_f64() * 1000.0;

  let metrics = ConvolutionMetrics::new(&input_grid, &actual_values, &expected_values, execution_time)?;

  let f_values = factors.first().cloned().unwrap_or_else(|| Array1::ones(input_grid.len()));
  let g_values = factors.get(1).cloned().unwrap_or_else(|| Array1::ones(input_grid.len()));

  Ok(TestResult {
    algorithm: algo.name().to_owned(),
    test_case: test_case.clone(),
    execution_time_ms: execution_time,
    f_x_values: input_grid.clone(),
    f_y_values: f_values,
    g_x_values: input_grid.clone(),
    g_y_values: g_values,
    evaluation_grid: input_grid,
    actual_values,
    expected_values,
    metrics,
    log_scale_actual: Some(actual_log_scale),
    log_scale_expected: Some(expected_log_scale),
    log_scale_error: Some(log_scale_error),
  })
}

fn generate_chain_multiplication_summary<S>(
  suite: &S,
  outcomes: &[TestRunOutcome<S::TestCase>],
  algorithms: &[MultiplicationAlgorithm],
) -> TestSummary
where
  S: ChainMultiplicationTestSuite,
{
  generate_multiplication_summary_from_outcomes(suite.test_suite_name(), outcomes, algorithms)
}
