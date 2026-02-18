use crate::testing::console::console::ValidationConsole;
use crate::testing::framework::results::{TestFailure, TestResult, TestRunOutcome};
use crate::testing::framework::summary::{AlgorithmSummary, TestSummary};
use crate::testing::framework::test_case::TestCase;
use crate::testing::framework::tsv_output::generate_tsv_outputs;
use crate::testing::plots::plots::generate_plot_outputs;
use crate::testing::run::Args;
use eyre::Report;
use itertools::Itertools;
use rayon::prelude::*;
use serde::Serialize;
use std::collections::BTreeSet;
use std::fmt::Display;
use std::fs;
use std::time::Instant;
use treetime_io::json::{JsonPretty, json_write_file};
use treetime_utils::make_error;

pub trait TestRunner: Send + Sync + Sized {
  type TestCase: TestCase;
  type Algorithm: Copy + Display + Send + Sync;
  type Suite: Send + Sync;

  fn test_suite_name(suite: &Self::Suite) -> &'static str;

  fn create_test_cases(suite: &Self::Suite) -> Vec<Self::TestCase>;

  fn get_algorithms(args: &Args) -> &[Self::Algorithm];

  fn run_single_test(
    suite: &Self::Suite,
    test_case: &Self::TestCase,
    algorithm: Self::Algorithm,
  ) -> Result<TestResult<Self::TestCase>, Report>;

  fn print_test_configuration(
    test_suite_name: &str,
    algorithms: &[Self::Algorithm],
    total_available: usize,
    selected_count: usize,
    name_filter_applied: bool,
    slowness_threshold: f64,
    output_dir: &str,
  );

  fn print_header(test_cases_count: usize, algorithms_count: usize);

  fn print_table_header();

  fn print_success_row(result: &TestResult<Self::TestCase>);

  fn print_failure_row(test_case: &Self::TestCase, algorithm: Self::Algorithm, elapsed_ms: f64);

  fn print_error_summary(failures: &[&TestFailure<Self::TestCase>]) -> Result<(), Report>;
}

#[allow(clippy::needless_pass_by_value)]
pub fn run_tests_generic<R: TestRunner>(
  args: &Args,
  suite: R::Suite,
  list_cases_fn: impl Fn(&R::Suite),
) -> Result<(), Report>
where
  R::Suite: Default,
{
  if args.list_cases {
    list_cases_fn(&suite);
    return Ok(());
  }

  let test_suite_name = R::test_suite_name(&suite);
  let output_dir = format!("{}/{}", args.output_dir, test_suite_name);
  let all_test_cases = R::create_test_cases(&suite);
  let algorithms = R::get_algorithms(args);

  let test_cases_after_name_filter = filter_test_cases_by_name(all_test_cases.clone(), Some(args.test_cases.as_str()))?;
  let selected_test_cases = filter_by_slowness(&test_cases_after_name_filter, args.slowness);

  R::print_test_configuration(
    test_suite_name,
    algorithms,
    all_test_cases.len(),
    selected_test_cases.len(),
    args.test_cases != "all",
    args.slowness,
    &output_dir,
  );

  let outcomes = run_all_tests::<R>(&suite, &all_test_cases, &selected_test_cases, algorithms, args.verbose)?;

  generate_plot_outputs(&output_dir, &outcomes)?;
  generate_tsv_outputs(&output_dir, &outcomes)?;

  let summary = generate_summary_generic(test_suite_name, &outcomes, algorithms);

  ValidationConsole::print_summary(&summary, &outcomes)?;

  save_results_json(&output_dir, &outcomes, &summary)?;

  Ok(())
}

fn filter_test_cases_by_name<T: TestCase + Clone>(all_cases: Vec<T>, filter: Option<&str>) -> Result<Vec<T>, Report> {
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

fn filter_by_slowness<T: TestCase>(test_cases: &[T], threshold: f64) -> Vec<T> {
  test_cases
    .iter()
    .filter(|case| case.slowness() <= threshold)
    .cloned()
    .collect()
}

fn run_all_tests<R: TestRunner>(
  suite: &R::Suite,
  all_test_cases: &[R::TestCase],
  selected_test_cases: &[R::TestCase],
  algorithms: &[R::Algorithm],
  verbose: bool,
) -> Result<Vec<TestRunOutcome<R::TestCase>>, Report> {
  R::print_header(selected_test_cases.len(), algorithms.len());
  R::print_table_header();

  let selected: Vec<_> = selected_test_cases
    .iter()
    .flat_map(|test_case| algorithms.iter().map(move |algorithm| (test_case, *algorithm)))
    .collect();

  let outcomes = selected
    .into_par_iter()
    .map(|(test_case, algorithm)| execute_single_test::<R>(suite, test_case, algorithm, verbose))
    .collect::<Vec<_>>();

  let failures = collect_failures(&outcomes);
  R::print_error_summary(&failures)?;

  Ok(outcomes)
}

fn execute_single_test<R: TestRunner>(
  suite: &R::Suite,
  test_case: &R::TestCase,
  algorithm: R::Algorithm,
  verbose: bool,
) -> TestRunOutcome<R::TestCase> {
  let start_time = Instant::now();

  match R::run_single_test(suite, test_case, algorithm) {
    Ok(result) => {
      R::print_success_row(&result);
      if verbose {
        ValidationConsole::print_verbose_details(&result);
      }
      TestRunOutcome::Success(Box::new(result))
    },
    Err(error) => {
      let elapsed_ms = start_time.elapsed().as_secs_f64() * 1000.0;
      R::print_failure_row(test_case, algorithm, elapsed_ms);
      TestRunOutcome::Failure(TestFailure {
        algorithm: algorithm.to_string(),
        test_case: test_case.clone(),
        error: format!("{error:?}"),
        execution_time_ms: elapsed_ms,
      })
    },
  }
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

fn collect_successes<T: TestCase>(outcomes: &[TestRunOutcome<T>]) -> Vec<&TestResult<T>> {
  outcomes
    .iter()
    .filter_map(|outcome| match outcome {
      TestRunOutcome::Success(result) => Some(result.as_ref()),
      TestRunOutcome::Failure(_) => None,
    })
    .collect()
}

fn calculate_total_execution_time<T: TestCase>(outcomes: &[TestRunOutcome<T>]) -> f64 {
  outcomes
    .iter()
    .map(|outcome| match outcome {
      TestRunOutcome::Success(result) => result.as_ref().execution_time_ms,
      TestRunOutcome::Failure(failure) => failure.execution_time_ms,
    })
    .sum()
}

pub fn generate_summary_generic<T: TestCase, A: Display>(
  test_suite_name: &str,
  outcomes: &[TestRunOutcome<T>],
  algorithms: &[A],
) -> TestSummary {
  let successes = collect_successes(outcomes);
  let failures = collect_failures(outcomes);
  let total_execution_time = calculate_total_execution_time(outcomes);
  let algorithm_summaries = build_algorithm_summaries_generic(algorithms, &successes, &failures);

  TestSummary {
    test_suite_name: test_suite_name.to_owned(),
    total_tests: outcomes.len(),
    total_successes: successes.len(),
    total_failures: failures.len(),
    total_algorithms: algorithms.len(),
    execution_time_total_ms: total_execution_time,
    algorithm_summaries,
  }
}

fn build_algorithm_summaries_generic<T: TestCase, A: Display>(
  algorithms: &[A],
  successes: &[&TestResult<T>],
  failures: &[&TestFailure<T>],
) -> Vec<AlgorithmSummary> {
  algorithms
    .iter()
    .filter_map(|algorithm| {
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
  json_write_file(&json_path, &results, JsonPretty(true))
}
