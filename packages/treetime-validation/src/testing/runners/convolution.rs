use crate::algorithms::Algo;
use crate::algorithms::ConvolutionAlgorithm;
use crate::testing::console::console::ValidationConsole;
use crate::testing::framework::results::{TestFailure, TestResult};
use crate::testing::framework::test_case::TestCase;
use crate::testing::metrics::metrics::ValidationMetrics;
use crate::testing::run::Args;
use crate::testing::runners::runner::TestRunner;
use crate::testing::test_suites::test_suites::TestSuite;
use eyre::Report;
use ndarray::Array1;
use std::time::Instant;

pub struct ConvolutionRunner<S: TestSuite>(std::marker::PhantomData<S>);

impl<S: TestSuite + Default> TestRunner for ConvolutionRunner<S> {
  type TestCase = S::TestCase;
  type Algorithm = ConvolutionAlgorithm;
  type Suite = S;

  fn test_suite_name(suite: &Self::Suite) -> &'static str {
    suite.test_suite_name()
  }

  fn create_test_cases(suite: &Self::Suite) -> Vec<Self::TestCase> {
    suite.create_test_cases()
  }

  fn get_algorithms(args: &Args) -> &[Self::Algorithm] {
    &args.algorithms
  }

  fn run_single_test(
    suite: &Self::Suite,
    test_case: &Self::TestCase,
    algorithm: Self::Algorithm,
  ) -> Result<TestResult<Self::TestCase>, Report> {
    let start_time = Instant::now();
    let algo = algorithm.instantiate()?;
    run_convolution_test(suite, test_case, &*algo, start_time)
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
    ValidationConsole::print_test_configuration(
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
    ValidationConsole::print_header(test_suite_name, test_cases_count, algorithms_count);
  }

  fn print_table_header() {
    ValidationConsole::print_progress_table_header();
  }

  fn print_success_row(result: &TestResult<Self::TestCase>, completed_count: usize, total_tests: usize) {
    ValidationConsole::print_success_row(result, completed_count, total_tests);
  }

  fn print_failure_row(
    test_case: &Self::TestCase,
    algorithm: Self::Algorithm,
    elapsed_ms: f64,
    completed_count: usize,
    total_tests: usize,
  ) {
    ValidationConsole::print_failure_row(test_case, algorithm, elapsed_ms, completed_count, total_tests);
  }

  fn print_skipped_row(test_case: &Self::TestCase, algorithm: Self::Algorithm) {
    ValidationConsole::print_skipped_row(test_case, algorithm);
  }

  fn print_error_summary(failures: &[&TestFailure<Self::TestCase>]) -> Result<(), Report> {
    ValidationConsole::print_error_summary(failures)
  }
}

fn run_convolution_test<S: TestSuite>(
  suite: &S,
  test_case: &S::TestCase,
  algo: &dyn Algo,
  start_time: Instant,
) -> Result<TestResult<S::TestCase>, Report> {
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

  let metrics = ValidationMetrics::new(&evaluation_grid, &actual_values, &expected_values, execution_time)?;

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
