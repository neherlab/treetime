use crate::distribution::reference::convolution_test::algorithm_summary::AlgorithmSummary;
use crate::distribution::reference::convolution_test::algorithms::ConvolutionAlgorithm;
use crate::distribution::reference::convolution_test::console::ConvolutionTestConsole;
use crate::distribution::reference::convolution_test::framework::results::{TestFailure, TestRunOutcome};
use crate::distribution::reference::convolution_test::framework::runner::ConvolutionTestRunner;
use crate::distribution::reference::convolution_test::framework::summary::TestSummary;
use crate::distribution::reference::convolution_test::framework::test_case::TestCase;
use crate::io::json::{JsonPretty, json_write_file};
use eyre::Report;
use itertools::Itertools;
use rayon::prelude::*;
use serde::Serialize;
use std::fs;
use std::marker::PhantomData;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

pub struct ConvolutionTestFramework<T: TestCase, R: ConvolutionTestRunner<T>> {
  pub runner: R,
  pub algorithms: Vec<ConvolutionAlgorithm>,
  pub output_dir: String,
  pub _phantom: PhantomData<T>,
}

impl<T: TestCase, R: ConvolutionTestRunner<T>> ConvolutionTestFramework<T, R> {
  pub fn new(runner: R, output_dir: String) -> Self {
    Self {
      runner,
      algorithms: ConvolutionAlgorithm::all(),
      output_dir,
      _phantom: PhantomData,
    }
  }

  pub fn set_algorithms(&mut self, algorithms: Vec<ConvolutionAlgorithm>) {
    self.algorithms = algorithms;
  }

  pub fn run_all_tests(&self) -> Result<Vec<TestRunOutcome<T>>, Report> {
    let total_tests = self.runner.test_cases().len() * self.algorithms.len();
    let completed = AtomicUsize::new(0);

    ConvolutionTestConsole::print_header(
      self.runner.function_type(),
      self.runner.test_cases().len(),
      self.algorithms.len(),
    );

    ConvolutionTestConsole::print_progress_table_header();

    let test_combinations: Vec<_> = self
      .runner
      .test_cases()
      .iter()
      .cartesian_product(&self.algorithms)
      .collect();

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

  pub fn print_summary(&self, summary: &TestSummary, outcomes: &[TestRunOutcome<T>]) {
    ConvolutionTestConsole::print_summary(summary, outcomes);
  }
}

#[derive(Serialize)]
struct FullResults<'a, T: TestCase> {
  summary: &'a TestSummary,
  outcomes: &'a [TestRunOutcome<T>],
}

impl<T: TestCase, R: ConvolutionTestRunner<T>> ConvolutionTestFramework<T, R> {
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
