use crate::testing::framework::results::{TestResult, TestRunOutcome};
use crate::testing::framework::test_case::TestCase;
use csv::WriterBuilder;
use eyre::Report;
use std::fs;
use treetime_utils::make_error;

pub fn generate_tsv_outputs<T>(output_dir: &str, outcomes: &[TestRunOutcome<T>]) -> Result<(), Report>
where
  T: TestCase,
{
  for outcome in outcomes {
    if let TestRunOutcome::Success(result) = outcome {
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
