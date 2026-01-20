use crate::testing::framework::results::{TestResult, TestRunOutcome};
use crate::testing::framework::test_case::TestCase;
use csv::WriterBuilder;
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use std::fs;
use treetime_utils::make_error;

pub fn generate_tsv_outputs<T>(output_dir: &str, outcomes: &[TestRunOutcome<T>]) -> Result<(), Report>
where
  T: TestCase,
{
  for outcome in outcomes {
    if let TestRunOutcome::Success(result) = outcome {
      let result = result.as_ref();
      let algorithm_dir = format!(
        "{}/{}/{}",
        output_dir,
        result.test_case.name(),
        result.algorithm.to_lowercase()
      );
      fs::create_dir_all(&algorithm_dir)?;

      write_reference_functions_tsv(result, &algorithm_dir)?;
      write_convolution_results_tsv(result, &algorithm_dir)?;
    }
  }

  write_comparative_summary_tsv(output_dir, outcomes)?;

  Ok(())
}

fn write_reference_functions_tsv<T>(result: &TestResult<T>, output_dir: &str) -> Result<(), Report>
where
  T: TestCase,
{
  let tsv_path = format!("{output_dir}/reference_functions.tsv");

  if result.f_x_values.len() != result.f_y_values.len() {
    return make_error!(
      "f_x_values length ({}) != f_y_values length ({})",
      result.f_x_values.len(),
      result.f_y_values.len()
    );
  }

  if result.g_x_values.len() != result.g_y_values.len() {
    return make_error!(
      "g_x_values length ({}) != g_y_values length ({})",
      result.g_x_values.len(),
      result.g_y_values.len()
    );
  }

  let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(tsv_path)?;

  writer.write_record(["x", "f(x)", "g(x)"])?;

  for i in 0..result.f_x_values.len() {
    let x = result.f_x_values[i];
    let f_val = result.f_y_values[i];
    let g_val = result.g_y_values[i];
    writer.write_record(&[x.to_string(), f_val.to_string(), g_val.to_string()])?;
  }

  writer.flush()?;
  Ok(())
}

fn write_convolution_results_tsv<T>(result: &TestResult<T>, output_dir: &str) -> Result<(), Report>
where
  T: TestCase,
{
  let tsv_path = format!("{output_dir}/convolution_results.tsv");

  if result.evaluation_grid.len() != result.expected_values.len() {
    return make_error!(
      "evaluation_grid length ({}) != expected_values length ({})",
      result.evaluation_grid.len(),
      result.expected_values.len()
    );
  }

  if result.evaluation_grid.len() != result.actual_values.len() {
    return make_error!(
      "evaluation_grid length ({}) != actual_values length ({})",
      result.evaluation_grid.len(),
      result.actual_values.len()
    );
  }

  let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(tsv_path)?;

  writer.write_record(["t", "expected", "actual", "absolute_error", "relative_error"])?;

  for i in 0..result.evaluation_grid.len() {
    let t = result.evaluation_grid[i];
    let expected = result.expected_values[i];
    let actual = result.actual_values[i];
    let abs_error = (actual - expected).abs();
    let rel_error = if expected.abs() > 1e-15 {
      abs_error / expected.abs()
    } else {
      0.0
    };

    writer.write_record(&[
      t.to_string(),
      expected.to_string(),
      actual.to_string(),
      abs_error.to_string(),
      rel_error.to_string(),
    ])?;
  }

  writer.flush()?;
  Ok(())
}

fn write_comparative_summary_tsv<T>(output_dir: &str, outcomes: &[TestRunOutcome<T>]) -> Result<(), Report>
where
  T: TestCase,
{
  let successes: Vec<_> = outcomes
    .iter()
    .filter_map(|outcome| match outcome {
      TestRunOutcome::Success(result) => Some(result.as_ref()),
      TestRunOutcome::Failure(_) => None,
    })
    .collect();

  if successes.is_empty() {
    return Ok(());
  }

  let algorithms: Vec<_> = successes
    .iter()
    .map(|r| r.algorithm.as_str())
    .unique()
    .sorted()
    .collect();

  let grouped_by_test_case: BTreeMap<&str, Vec<&TestResult<T>>> = successes
    .iter()
    .copied()
    .into_group_map_by(|r| r.test_case.name())
    .into_iter()
    .collect();

  fs::create_dir_all(output_dir)?;
  let tsv_path = format!("{output_dir}/comparative_summary.tsv");

  let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(&tsv_path)?;

  let mut headers = vec!["test_case".to_owned()];
  for algo in &algorithms {
    headers.push(format!("{algo}_time_ms"));
    headers.push(format!("{algo}_r2"));
    headers.push(format!("{algo}_rmse"));
    headers.push(format!("{algo}_correlation"));
    headers.push(format!("{algo}_mass_error"));
  }
  writer.write_record(&headers)?;

  for (test_case_name, results) in &grouped_by_test_case {
    let results_by_algo: BTreeMap<&str, &TestResult<T>> =
      results.iter().map(|r| (r.algorithm.as_str(), *r)).collect();

    let mut row = vec![(*test_case_name).to_owned()];

    for algo in &algorithms {
      if let Some(result) = results_by_algo.get(algo) {
        let metrics = &result.metrics.aggregate.domain_agreement.quality_metrics;
        row.push(result.execution_time_ms.to_string());
        row.push(metrics.r_squared.to_string());
        row.push(metrics.rmse.to_string());
        row.push(metrics.correlation.to_string());
        row.push(metrics.mass_error.to_string());
      } else {
        row.push(String::new());
        row.push(String::new());
        row.push(String::new());
        row.push(String::new());
        row.push(String::new());
      }
    }

    writer.write_record(&row)?;
  }

  writer.flush()?;
  Ok(())
}
