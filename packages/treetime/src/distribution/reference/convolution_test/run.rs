use crate::distribution::reference::convolution_test::algorithm_summary::AlgorithmSummary;
use crate::distribution::reference::convolution_test::algorithms::ConvolutionAlgorithm;
use crate::distribution::reference::convolution_test::console::ConvolutionTestConsole;
use crate::distribution::reference::convolution_test::framework::results::{TestFailure, TestRunOutcome};
use crate::distribution::reference::convolution_test::framework::summary::TestSummary;
use crate::distribution::reference::convolution_test::framework::test_case::TestCase;
use crate::distribution::reference::convolution_test::functions::exponential::ExponentialConvInput;
use crate::distribution::reference::convolution_test::functions::functions::FunctionType;
use crate::distribution::reference::convolution_test::functions::gaussian::GaussianConvInput;
use crate::distribution::reference::convolution_test::plots::plots::generate_plot_outputs;
use crate::distribution::reference::convolution_test::traits::ConvInput;
use crate::io::json::{JsonPretty, json_write_file};
use clap::Parser;
use eyre::Report;
use itertools::Itertools;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::fs;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

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
  #[arg(long, default_value = "tmp/convolution_test")]
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

fn run_convolution_tests_impl<I>(args: &Args) -> Result<(), Report>
where
  I: ConvInput,
{
  let input = I::new(&args.test_cases)?;

  if args.list_cases {
    input.list_test_cases();
    return Ok(());
  }

  let output_dir = format!("{}/{}", args.output_dir, input.function_type());
  let test_cases = input.filter_test_cases(Some(args.test_cases.as_str()))?;

  if args.verbose {
    println!("Test Configuration:");
    println!("  Function type: {}", input.function_type());
    println!("  Algorithms: {:?}", args.algorithms);
    println!("  Test cases: {} selected", test_cases.len());
    println!("  Output directory: {output_dir}\n");
  }

  let outcomes = run_all_tests(&input, &test_cases, &args.algorithms)?;

  generate_plot_outputs(&output_dir, &outcomes)?;

  let summary = generate_summary(&input, &outcomes, &args.algorithms);

  ConvolutionTestConsole::print_summary(&summary, &outcomes);

  save_results_json(&output_dir, &outcomes, &summary)?;

  println!(
    "{} convolution test framework completed successfully!",
    input.function_type()
  );
  println!("Check {output_dir} for detailed results.");

  Ok(())
}

fn run_all_tests<I>(
  input: &I,
  test_cases: &[I::TestCase],
  algorithms: &[ConvolutionAlgorithm],
) -> Result<Vec<TestRunOutcome<I::TestCase>>, Report>
where
  I: ConvInput,
{
  let total_tests = test_cases.len() * algorithms.len();
  let completed = AtomicUsize::new(0);

  ConvolutionTestConsole::print_header(
    input.function_type(),
    test_cases.len(),
    algorithms.len(),
  );

  ConvolutionTestConsole::print_progress_table_header();

  let test_combinations: Vec<_> = test_cases
    .iter()
    .cartesian_product(algorithms)
    .collect();

  let outcomes: Vec<TestRunOutcome<I::TestCase>> = test_combinations
    .into_par_iter()
    .map(|(test_case, &algorithm)| {
      let start_time = Instant::now();

      match input.run_test(test_case, algorithm) {
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

fn generate_summary<I>(
  input: &I,
  outcomes: &[TestRunOutcome<I::TestCase>],
  algorithms: &[ConvolutionAlgorithm],
) -> TestSummary
where
  I: ConvInput,
{
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

  for &algorithm in algorithms {
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
    function_type: input.function_type().to_owned(),
    total_tests: outcomes.len(),
    total_successes: successes.len(),
    total_failures: failures.len(),
    total_algorithms: algorithms.len(),
    execution_time_total_ms: total_execution_time,
    algorithm_summaries,
    overall_assessment,
  }
}

fn save_results_json<T>(
  output_dir: &str,
  outcomes: &[TestRunOutcome<T>],
  summary: &TestSummary,
) -> Result<(), Report>
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
