use crate::distribution::reference::convolution_test::framework::{TestCase, TestResult};
use crate::io::csv::CsvStructFileWriter;
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::fs;

/// Trait for flattening test results to TSV-compatible format
pub trait ToFlatResult {
  type FlatResult: Serialize;

  /// Convert test result to flat structure for TSV output
  fn to_flat_result(&self) -> Self::FlatResult;
}

/// Generic output utilities for test frameworks
#[derive(Clone, Debug)]
pub struct TestOutputWriter {
  output_dir: String,
}

impl TestOutputWriter {
  pub fn new(output_dir: String) -> Self {
    Self { output_dir }
  }

  /// Save results to TSV file for analysis in Python/R
  pub fn save_results_tsv<F: Serialize>(&self, flat_results: &[F]) -> Result<(), Report> {
    fs::create_dir_all(&self.output_dir)?;

    let output_dir = &self.output_dir;
    let tsv_path = format!("{output_dir}/convolution_test_results.tsv");
    let mut writer = CsvStructFileWriter::new(&tsv_path, b'\t')?;

    for flat_result in flat_results {
      writer.write(flat_result)?;
    }

    println!("Saved TSV results to: {tsv_path}");
    Ok(())
  }
}

/// Base flat result structure with common fields
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BaseFlatResult {
  pub test_case_name: String,
  pub algorithm: String,
  pub dx: f64,
  pub grid_points: usize,
  pub execution_time_ms: f64,
  pub r_squared: f64,
  pub max_abs_error: f64,
  pub mean_abs_error: f64,
  pub max_rel_error_percent: f64,
  pub mean_rel_error_percent: f64,
  pub rmse: f64,
  pub correlation: f64,
  pub mass_conservation_error: f64,
  pub peak_error_x: f64,
  pub peak_error_abs: f64,
  pub tolerance_strict_abs: f64,
  pub tolerance_moderate_abs: f64,
  pub tolerance_loose_abs: f64,
  pub tolerance_strict_rel: f64,
  pub tolerance_moderate_rel: f64,
  pub tolerance_loose_rel: f64,
  pub stress_type: String,
  pub overall_assessment: String,
}

impl BaseFlatResult {
  /// Create base flat result from any test result
  pub fn from_test_result<T: TestCase>(result: &TestResult<T>) -> Self {
    Self {
      test_case_name: result.test_case.name().to_owned(),
      algorithm: result.algorithm.to_string(),
      dx: result.test_case.dx(),
      grid_points: result.evaluation_grid.len(),
      execution_time_ms: result.execution_time_ms,
      r_squared: result.metrics.quality_metrics.r_squared,
      max_abs_error: result.metrics.abs_error_stats.max,
      mean_abs_error: result.metrics.abs_error_stats.mean,
      max_rel_error_percent: result.metrics.rel_error_stats.max * 100.0,
      mean_rel_error_percent: result.metrics.rel_error_stats.mean * 100.0,
      rmse: result.metrics.quality_metrics.rmse,
      correlation: result.metrics.quality_metrics.correlation,
      mass_conservation_error: result.metrics.quality_metrics.mass_error,
      peak_error_x: result.metrics.max_error_location.x_value,
      peak_error_abs: result.metrics.abs_error_stats.max,
      tolerance_strict_abs: result.metrics.abs_tolerance_fraction(0),
      tolerance_moderate_abs: result.metrics.abs_tolerance_fraction(1),
      tolerance_loose_abs: result.metrics.abs_tolerance_fraction(2),
      tolerance_strict_rel: result.metrics.rel_tolerance_fraction(0),
      tolerance_moderate_rel: result.metrics.rel_tolerance_fraction(1),
      tolerance_loose_rel: result.metrics.rel_tolerance_fraction(2),
      stress_type: result.test_case.stress_type().to_owned(),
      overall_assessment: result.metrics.overall_assessment().to_string(),
    }
  }
}
