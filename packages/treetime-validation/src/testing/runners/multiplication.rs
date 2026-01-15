use crate::algorithms::MultiplicationAlgorithm;
use crate::algorithms::MultiplyAlgo;
use crate::testing::framework::results::{TestFailure, TestResult};
use crate::testing::framework::test_case::TestCase;
use crate::testing::metrics::metrics::ValidationMetrics;
use crate::testing::run::Args;
use crate::testing::runners::runner::TestRunner;
use crate::testing::test_suites::test_suites::ChainMultiplicationTestSuite;
use crate::testing::test_suites::test_suites::MultiplicationTestSuite;
use eyre::Report;
use ndarray::Array1;
use std::fmt::Display;
use std::time::Instant;
use treetime_ops::ScaledArray;

pub struct MultiplicationRunner<S: MultiplicationTestSuite>(std::marker::PhantomData<S>);

impl<S: MultiplicationTestSuite + Default> TestRunner for MultiplicationRunner<S> {
  type TestCase = S::TestCase;
  type Algorithm = MultiplicationAlgorithm;
  type Suite = S;

  fn test_suite_name(suite: &Self::Suite) -> &'static str {
    suite.test_suite_name()
  }

  fn create_test_cases(suite: &Self::Suite) -> Vec<Self::TestCase> {
    suite.create_test_cases()
  }

  fn get_algorithms(args: &Args) -> &[Self::Algorithm] {
    &args.mult_algorithms
  }

  fn run_single_test(
    suite: &Self::Suite,
    test_case: &Self::TestCase,
    algorithm: Self::Algorithm,
  ) -> Result<TestResult<Self::TestCase>, Report> {
    let start_time = Instant::now();
    let algo = algorithm.instantiate()?;
    run_multiplication_test(suite, test_case, &*algo, start_time)
  }

  fn print_test_configuration(
    test_suite_name: &str,
    algorithms: &[Self::Algorithm],
    total_available: usize,
    selected_count: usize,
    name_filter_applied: bool,
    slowness_threshold: f64,
    output_dir: &str,
  ) {
    print_multiplication_config(
      test_suite_name,
      algorithms,
      total_available,
      selected_count,
      name_filter_applied,
      slowness_threshold,
      output_dir,
    );
  }

  fn print_header(test_suite_name: &str, test_cases_count: usize, algorithms_count: usize) {
    println!(
      "\nRunning {test_suite_name} multiplication tests ({test_cases_count} cases x {algorithms_count} algorithms)...\n"
    );
  }

  fn print_table_header() {
    println!(
      "{:<30} {:<25} {:>10} {:>15}",
      "Test Case", "Algorithm", "Time (ms)", "Status"
    );
    println!("{:-<80}", "");
  }

  fn print_success_row(result: &TestResult<Self::TestCase>, completed_count: usize, total_tests: usize) {
    println!(
      "{:<30} {:<25} {:>10.2} {:>15} [{}/{}]",
      result.test_case.name(),
      result.algorithm,
      result.execution_time_ms,
      "OK",
      completed_count,
      total_tests
    );
  }

  fn print_failure_row(
    test_case: &Self::TestCase,
    algorithm: Self::Algorithm,
    elapsed_ms: f64,
    completed_count: usize,
    total_tests: usize,
  ) {
    println!(
      "{:<30} {:<25} {:>10.2} {:>15} [{}/{}]",
      test_case.name(),
      algorithm,
      elapsed_ms,
      "FAIL",
      completed_count,
      total_tests
    );
  }

  fn print_skipped_row(test_case: &Self::TestCase, algorithm: Self::Algorithm) {
    println!(
      "{:<30} {:<25} {:>10} {:>15}",
      test_case.name(),
      algorithm,
      "-",
      "SKIPPED"
    );
  }

  fn print_error_summary(failures: &[&TestFailure<Self::TestCase>]) -> Result<(), Report> {
    print_multiplication_errors(failures)
  }
}

pub struct ChainMultiplicationRunner<S: ChainMultiplicationTestSuite>(std::marker::PhantomData<S>);

impl<S: ChainMultiplicationTestSuite + Default> TestRunner for ChainMultiplicationRunner<S> {
  type TestCase = S::TestCase;
  type Algorithm = MultiplicationAlgorithm;
  type Suite = S;

  fn test_suite_name(suite: &Self::Suite) -> &'static str {
    suite.test_suite_name()
  }

  fn create_test_cases(suite: &Self::Suite) -> Vec<Self::TestCase> {
    suite.create_test_cases()
  }

  fn get_algorithms(args: &Args) -> &[Self::Algorithm] {
    &args.mult_algorithms
  }

  fn run_single_test(
    suite: &Self::Suite,
    test_case: &Self::TestCase,
    algorithm: Self::Algorithm,
  ) -> Result<TestResult<Self::TestCase>, Report> {
    let start_time = Instant::now();
    let algo = algorithm.instantiate()?;
    run_chain_multiplication_test(suite, test_case, &*algo, start_time)
  }

  fn print_test_configuration(
    test_suite_name: &str,
    algorithms: &[Self::Algorithm],
    total_available: usize,
    selected_count: usize,
    name_filter_applied: bool,
    slowness_threshold: f64,
    output_dir: &str,
  ) {
    print_multiplication_config(
      test_suite_name,
      algorithms,
      total_available,
      selected_count,
      name_filter_applied,
      slowness_threshold,
      output_dir,
    );
  }

  fn print_header(test_suite_name: &str, test_cases_count: usize, algorithms_count: usize) {
    println!(
      "\nRunning {test_suite_name} chain multiplication tests ({test_cases_count} cases x {algorithms_count} algorithms)...\n"
    );
  }

  fn print_table_header() {
    println!(
      "{:<30} {:<25} {:>10} {:>15}",
      "Test Case", "Algorithm", "Time (ms)", "Status"
    );
    println!("{:-<80}", "");
  }

  fn print_success_row(result: &TestResult<Self::TestCase>, completed_count: usize, total_tests: usize) {
    println!(
      "{:<30} {:<25} {:>10.2} {:>15} [{}/{}]",
      result.test_case.name(),
      result.algorithm,
      result.execution_time_ms,
      "OK",
      completed_count,
      total_tests
    );
  }

  fn print_failure_row(
    test_case: &Self::TestCase,
    algorithm: Self::Algorithm,
    elapsed_ms: f64,
    completed_count: usize,
    total_tests: usize,
  ) {
    println!(
      "{:<30} {:<25} {:>10.2} {:>15} [{}/{}]",
      test_case.name(),
      algorithm,
      elapsed_ms,
      "FAIL",
      completed_count,
      total_tests
    );
  }

  fn print_skipped_row(test_case: &Self::TestCase, algorithm: Self::Algorithm) {
    println!(
      "{:<30} {:<25} {:>10} {:>15}",
      test_case.name(),
      algorithm,
      "-",
      "SKIPPED"
    );
  }

  fn print_error_summary(failures: &[&TestFailure<Self::TestCase>]) -> Result<(), Report> {
    print_multiplication_errors(failures)
  }
}

fn print_multiplication_config<A: Display>(
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

fn print_multiplication_errors<T: TestCase>(failures: &[&TestFailure<T>]) -> Result<(), Report> {
  if failures.is_empty() {
    return Ok(());
  }

  println!("\n=== Errors ===");
  for failure in failures {
    println!(
      "  {} ({}): {}",
      failure.test_case.name(),
      failure.algorithm,
      failure.error
    );
  }
  Ok(())
}

fn run_multiplication_test<S: MultiplicationTestSuite>(
  suite: &S,
  test_case: &S::TestCase,
  algo: &dyn MultiplyAlgo,
  start_time: Instant,
) -> Result<TestResult<S::TestCase>, Report> {
  let (input_grid_min, input_grid_max) = test_case.input_grid_domain();
  let input_grid_n_points = test_case.input_grid_n_points();
  let input_grid = Array1::linspace(input_grid_min, input_grid_max, input_grid_n_points);

  let f_values = suite.create_f(test_case, &input_grid)?;
  let g_values = suite.create_g(test_case, &input_grid)?;

  let actual_values = algo.multiply(&f_values, &g_values);
  let expected = suite.analytical_multiplication(test_case, &input_grid)?;

  let expected_values = expected.normalized.mapv(|v| v * expected.log_scale.exp());

  let execution_time = start_time.elapsed().as_secs_f64() * 1000.0;

  let metrics = ValidationMetrics::new(&input_grid, &actual_values, &expected_values, execution_time)?;

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

fn run_chain_multiplication_test<S: ChainMultiplicationTestSuite>(
  suite: &S,
  test_case: &S::TestCase,
  algo: &dyn MultiplyAlgo,
  start_time: Instant,
) -> Result<TestResult<S::TestCase>, Report> {
  let (input_grid_min, input_grid_max) = test_case.input_grid_domain();
  let input_grid_n_points = test_case.input_grid_n_points();
  let input_grid = Array1::linspace(input_grid_min, input_grid_max, input_grid_n_points);

  let factors = suite.create_factors(test_case, &input_grid)?;
  let factor_refs: Vec<&Array1<f64>> = factors.iter().collect();

  let actual = if factors.is_empty() {
    ScaledArray::new(Array1::ones(input_grid.len()), 0.0)
  } else {
    algo.multiply_many(&factor_refs)
  };
  let expected = suite.analytical_chain_multiplication(test_case, &input_grid)?;

  let log_scale_error = (actual.log_scale - expected.log_scale).abs();

  let actual_values = actual.normalized.mapv(|v| v * actual.log_scale.exp());
  let expected_values = expected.normalized.mapv(|v| v * expected.log_scale.exp());

  let execution_time = start_time.elapsed().as_secs_f64() * 1000.0;

  let metrics = ValidationMetrics::new(&input_grid, &actual_values, &expected_values, execution_time)?;

  let f_values = factors
    .first()
    .cloned()
    .unwrap_or_else(|| Array1::ones(input_grid.len()));
  let g_values = factors
    .get(1)
    .cloned()
    .unwrap_or_else(|| Array1::ones(input_grid.len()));

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
    log_scale_actual: Some(actual.log_scale),
    log_scale_expected: Some(expected.log_scale),
    log_scale_error: Some(log_scale_error),
  })
}
