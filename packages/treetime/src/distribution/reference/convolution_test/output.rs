use crate::distribution::reference::convolution_test::framework::{TestCase, TestResult};
use crate::io::csv::CsvStructFileWriter;
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::fs;

/// Row structure for pointwise error data export
#[derive(Serialize)]
struct PointwiseErrorRow {
  x: f64,
  abs_error: f64,
  rel_error: f64,
  signed_error: f64,
  log_error: f64,
  d1_error: f64,
  d2_error: f64,
  symmetry_residual: f64,
  monotonicity_violation: f64,
  peak_region_error: f64,
  tail_region_error: f64,
  cumulative_error: f64,
  support_coverage_mismatch: f64,
  tolerance_pass_strict: f64,
  tolerance_pass_moderate: f64,
  tolerance_pass_loose: f64,
}

/// Row structure for spatial metrics data export
#[derive(Serialize)]
struct SpatialMetricsRow {
  x: f64,
  sliding_rms: f64,
  sliding_max: f64,
}

/// Row structure for error histogram data export
#[derive(Serialize)]
struct HistogramRow {
  bin_center: f64,
  bin_count: usize,
}

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

  /// Save pointwise array data grouped by evaluation grid to TSV files
  pub fn save_pointwise_arrays<T: TestCase>(&self, result: &TestResult<T>, output_dir: &str) -> Result<(), Report> {
    fs::create_dir_all(output_dir)?;

    // Group 1: Evaluation grid-based arrays (main errors)
    self.save_evaluation_grid_arrays(result, output_dir)?;

    // Group 2: Sliding window arrays (may have different grid due to boundary handling)
    self.save_spatial_arrays(result, output_dir)?;

    Ok(())
  }

  /// Save arrays aligned with evaluation grid
  fn save_evaluation_grid_arrays<T: TestCase>(&self, result: &TestResult<T>, output_dir: &str) -> Result<(), Report> {
    let pw = &result.metrics.pointwise;
    let sp = &result.metrics.spatial;
    let tsv_path = format!("{output_dir}/pointwise_errors.tsv");
    let mut writer = CsvStructFileWriter::new(&tsv_path, b'\t')?;

    for i in 0..pw.total_points {
      writer.write(&PointwiseErrorRow {
        x: result.evaluation_grid[i],
        abs_error: pw.errors.absolute[i],
        rel_error: pw.errors.relative[i],
        signed_error: pw.errors.signed[i],
        log_error: pw.errors.logarithmic[i],
        d1_error: pw.structural.first_derivative[i],
        d2_error: pw.structural.second_derivative[i],
        symmetry_residual: pw.structural.symmetry_residual[i],
        monotonicity_violation: pw.structural.monotonicity_violations[i],
        peak_region_error: sp.regional.peak_region_errors[i],
        tail_region_error: sp.regional.tail_region_errors[i],
        cumulative_error: sp.cumulative.cumulative_error[i],
        support_coverage_mismatch: pw.tolerance.support_coverage_mask[i],
        tolerance_pass_strict: pw.tolerance.pass_masks[0][i],
        tolerance_pass_moderate: pw.tolerance.pass_masks[1][i],
        tolerance_pass_loose: pw.tolerance.pass_masks[2][i],
      })?;
    }

    Ok(())
  }

  /// Save spatial (sliding window) arrays
  fn save_spatial_arrays<T: TestCase>(&self, result: &TestResult<T>, output_dir: &str) -> Result<(), Report> {
    let sp = &result.metrics.spatial;
    let tsv_path = format!("{output_dir}/pointwise_spatial.tsv");
    let mut writer = CsvStructFileWriter::new(&tsv_path, b'\t')?;

    for i in 0..sp.total_points {
      writer.write(&SpatialMetricsRow {
        x: result.evaluation_grid[i],
        sliding_rms: sp.windowed.sliding_rms[i],
        sliding_max: sp.windowed.sliding_max[i],
      })?;
    }

    Ok(())
  }

  /// Save error histogram data
  pub fn save_error_histogram<T: TestCase>(&self, result: &TestResult<T>, output_dir: &str) -> Result<(), Report> {
    let dist = &result.metrics.distribution;
    let tsv_path = format!("{output_dir}/pointwise_histogram.tsv");
    let mut writer = CsvStructFileWriter::new(&tsv_path, b'\t')?;

    for (center, count) in dist
      .histograms
      .abs_error_histogram
      .bin_centers
      .iter()
      .zip(dist.histograms.abs_error_histogram.bin_counts.iter())
    {
      writer.write(&HistogramRow {
        bin_center: *center,
        bin_count: *count,
      })?;
    }

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
  // Pointwise error summaries
  pub pw_abs_mean: f64,
  pub pw_abs_max: f64,
  pub pw_abs_std: f64,
  pub pw_rel_mean: f64,
  pub pw_rel_max: f64,
  pub pw_rel_median: f64,
  pub pw_signed_bias: f64,
  pub pw_log_valid_count: usize,
  pub pw_log_max: f64,
  // Structural error summaries
  pub pw_d1_max: f64,
  pub pw_d1_mean: f64,
  pub pw_d2_max: f64,
  pub pw_d2_mean: f64,
  pub pw_symmetry_max: f64,
  pub pw_monotonicity_violations: usize,
  // Spatial summaries
  pub pw_peak_region_count: usize,
  pub pw_peak_region_mean: f64,
  pub pw_peak_region_max: f64,
  pub pw_tail_region_count: usize,
  pub pw_tail_region_mean: f64,
  pub pw_tail_region_max: f64,
  pub pw_sliding_rms_max: f64,
  pub pw_sliding_max_max: f64,
  pub pw_cumulative_final: f64,
  pub pw_cumulative_max_abs: f64,
  // Distribution summaries
  pub pw_support_mismatch_count: usize,
  pub pw_tolerance_pass_strict: f64,
  pub pw_tolerance_pass_moderate: f64,
  pub pw_tolerance_pass_loose: f64,
  pub pw_error_median: f64,
  pub pw_error_q95: f64,
  pub pw_error_q99: f64,
}

impl BaseFlatResult {
  /// Create base flat result from any test result
  pub fn from_test_result<T: TestCase>(result: &TestResult<T>) -> Self {
    let aggregate = &result.metrics.aggregate;
    let pointwise = &result.metrics.pointwise;
    let spatial = &result.metrics.spatial;
    let distribution = &result.metrics.distribution;
    Self {
      test_case_name: result.test_case.name().to_owned(),
      algorithm: result.algorithm.to_string(),
      dx: result.test_case.dx(),
      grid_points: result.evaluation_grid.len(),
      execution_time_ms: result.execution_time_ms,
      r_squared: aggregate.domain_agreement.quality_metrics.r_squared,
      max_abs_error: aggregate.domain_agreement.abs_error_stats.max,
      mean_abs_error: aggregate.domain_agreement.abs_error_stats.mean,
      max_rel_error_percent: aggregate.domain_agreement.rel_error_stats.max * 100.0,
      mean_rel_error_percent: aggregate.domain_agreement.rel_error_stats.mean * 100.0,
      rmse: aggregate.domain_agreement.quality_metrics.rmse,
      correlation: aggregate.domain_agreement.quality_metrics.correlation,
      mass_conservation_error: aggregate.domain_agreement.quality_metrics.mass_error,
      peak_error_x: aggregate.domain_agreement.max_error_location.x_value,
      peak_error_abs: aggregate.domain_agreement.abs_error_stats.max,
      tolerance_strict_abs: aggregate.domain_agreement.abs_tolerance_fraction(0),
      tolerance_moderate_abs: aggregate.domain_agreement.abs_tolerance_fraction(1),
      tolerance_loose_abs: aggregate.domain_agreement.abs_tolerance_fraction(2),
      tolerance_strict_rel: aggregate.domain_agreement.rel_tolerance_fraction(0),
      tolerance_moderate_rel: aggregate.domain_agreement.rel_tolerance_fraction(1),
      tolerance_loose_rel: aggregate.domain_agreement.rel_tolerance_fraction(2),
      stress_type: result.test_case.stress_type().to_owned(),
      overall_assessment: aggregate.domain_agreement.overall_assessment().to_string(),
      // Pointwise error summaries
      pw_abs_mean: pointwise.errors.summary.abs_mean,
      pw_abs_max: pointwise.errors.summary.abs_max,
      pw_abs_std: pointwise.errors.summary.abs_std,
      pw_rel_mean: pointwise.errors.summary.rel_mean,
      pw_rel_max: pointwise.errors.summary.rel_max,
      pw_rel_median: pointwise.errors.summary.rel_median,
      pw_signed_bias: pointwise.errors.summary.signed_bias,
      pw_log_valid_count: pointwise.errors.summary.log_valid_count,
      pw_log_max: pointwise.errors.summary.log_max,
      // Structural error summaries
      pw_d1_max: pointwise.structural.summary.d1_max,
      pw_d1_mean: pointwise.structural.summary.d1_mean,
      pw_d2_max: pointwise.structural.summary.d2_max,
      pw_d2_mean: pointwise.structural.summary.d2_mean,
      pw_symmetry_max: pointwise.structural.summary.symmetry_max,
      pw_monotonicity_violations: pointwise.structural.summary.monotonicity_violation_count,
      // Spatial summaries
      pw_peak_region_count: spatial.regional.summary.peak_region.point_count,
      pw_peak_region_mean: spatial.regional.summary.peak_region.mean_error,
      pw_peak_region_max: spatial.regional.summary.peak_region.max_error,
      pw_tail_region_count: spatial.regional.summary.tail_region.point_count,
      pw_tail_region_mean: spatial.regional.summary.tail_region.mean_error,
      pw_tail_region_max: spatial.regional.summary.tail_region.max_error,
      pw_sliding_rms_max: spatial.windowed.summary.sliding_rms_max,
      pw_sliding_max_max: spatial.windowed.summary.sliding_max_max,
      pw_cumulative_final: spatial.cumulative.summary.final_value,
      pw_cumulative_max_abs: spatial.cumulative.summary.max_abs,
      // Distribution summaries
      pw_support_mismatch_count: pointwise.tolerance.summary.support_mismatch_count,
      pw_tolerance_pass_strict: pointwise.tolerance.summary.pass_fractions[0],
      pw_tolerance_pass_moderate: pointwise.tolerance.summary.pass_fractions[1],
      pw_tolerance_pass_loose: pointwise.tolerance.summary.pass_fractions[2],
      pw_error_median: distribution.statistics.abs_error_stats.median,
      pw_error_q95: distribution.statistics.abs_error_stats.q95,
      pw_error_q99: distribution.statistics.abs_error_stats.q99,
    }
  }
}
