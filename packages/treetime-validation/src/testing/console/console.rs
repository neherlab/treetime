use crate::testing::framework::results::{TestFailure, TestResult, TestRunOutcome};
use crate::testing::framework::summary::TestSummary;
use crate::testing::framework::test_case::TestCase;
use eyre::Report;
use itertools::izip;
use std::fmt::Display;
use treetime_io::json::{JsonPretty, json_write_str};
use treetime_utils::float_fmt::float_to_significant_digits;

const VERBOSE_LABEL_WIDTH: usize = 32;

/// Console display functionality for convolution test framework
pub struct ValidationConsole;

impl ValidationConsole {
  /// Print test configuration
  pub fn print_test_configuration<A: Display>(
    test_suite_name: &str,
    algorithms: &[A],
    total_available: usize,
    selected_count: usize,
    name_filter_applied: bool,
    slowness_threshold: f64,
    output_dir: &str,
  ) {
    println!("# {}\n", test_suite_name.to_uppercase());
    println!("## Configuration\n");
    println!("- Test suite: {test_suite_name}");
    println!(
      "- Algorithms: {}",
      algorithms.iter().map(|a| a.to_string()).collect::<Vec<_>>().join(", ")
    );
    println!("- Test cases: running {selected_count}/{total_available} test cases");
    if name_filter_applied {
      println!("  - --test-cases name filter applied");
    }
    if slowness_threshold < 1.0 {
      println!("  - --slowness threshold: {slowness_threshold}");
    }
    let skipped = total_available - selected_count;
    if skipped > 0 {
      println!("  - ⚠️ {skipped} test cases will be skipped");
    }
    println!("- Output directory: {output_dir}");
    println!();
  }

  /// Print test framework header
  pub fn print_header(test_cases_count: usize, algorithms_count: usize) {
    let total_tests = test_cases_count * algorithms_count;
    println!("## Execution\n");
    println!("Running {test_cases_count} test cases with {algorithms_count} algorithms ({total_tests} total)\n",);
  }

  /// Print progress table header
  pub fn print_progress_table_header() {
    println!(
      "| S  | {:^25} | {:^30} | {:^8} | {:^10} | {:^8} | {:^9} |",
      "Algorithm", "Test Case", "Time, ms", "R² err ppm", "RMSE", "Corr err"
    );
    println!(
      "|----|{:-<27}|{:-<32}|{:->10}|{:->12}|{:->10}|{:->11}|",
      "", "", "", "", "", ""
    );
  }

  /// Print success row in progress table
  pub fn print_success_row<T: TestCase>(result: &TestResult<T>) {
    let elapsed_ms = result.execution_time_ms;
    let r_squared = result.metrics.aggregate.domain_agreement.quality_metrics.r_squared;
    let r2_error_ppm = (1.0 - r_squared) * 1_000_000.0;
    let rmse = result.metrics.aggregate.domain_agreement.quality_metrics.rmse;
    let correlation = result.metrics.aggregate.domain_agreement.quality_metrics.correlation;
    let correlation_error_ppm = (1.0 - correlation) * 1_000_000.0;
    println!(
      "| ok | {:<25} | {:<30} | {:>8.1} | {:>10} | {:>8} | {:>9} |",
      result.algorithm,
      result.test_case.name(),
      elapsed_ms,
      float_to_significant_digits(r2_error_ppm, 3),
      float_to_significant_digits(rmse, 3),
      float_to_significant_digits(correlation_error_ppm, 3)
    );
  }

  /// Print failure row in progress table
  pub fn print_failure_row<T: TestCase, A: Display>(test_case: &T, algorithm: A, elapsed_ms: f64) {
    println!(
      "| !! | {:<25} | {:<30} | {:>8.1} | {:>10} | {:>8} | {:>9} |",
      format!("{algorithm}"),
      test_case.name(),
      elapsed_ms,
      "FAILED",
      "FAILED",
      "FAILED"
    );
  }

  /// Print verbose details for a single test result
  pub fn print_verbose_details<T: TestCase>(result: &TestResult<T>) {
    let w = VERBOSE_LABEL_WIDTH;
    let m = &result.metrics;

    println!();
    println!(
      "    --- Verbose Details: {} + {} ---",
      result.test_case.name(),
      result.algorithm
    );

    println!("    Test Case Info:");
    println!("      {:<w$} {}", "Description:", result.test_case.description());
    println!("      {:<w$} {}", "Stress Type:", result.test_case.stress_type());
    println!(
      "      {:<w$} {}",
      "Analytical Caution:",
      result.test_case.analytical_caution()
    );
    let (domain_min, domain_max) = result.test_case.input_grid_domain();
    println!(
      "      {:<w$} [{}, {}] ({} points)",
      "Grid Domain:",
      float_to_significant_digits(domain_min, 4),
      float_to_significant_digits(domain_max, 4),
      result.test_case.input_grid_n_points()
    );

    println!("    Quality Metrics:");
    let qm = &m.aggregate.domain_agreement.quality_metrics;
    println!(
      "      {:<w$} {} (error: {} ppm)",
      "R-squared:",
      float_to_significant_digits(qm.r_squared, 6),
      float_to_significant_digits((1.0 - qm.r_squared) * 1_000_000.0, 3)
    );
    println!(
      "      {:<w$} {} (error: {} ppm)",
      "Correlation:",
      float_to_significant_digits(qm.correlation, 6),
      float_to_significant_digits((1.0 - qm.correlation) * 1_000_000.0, 3)
    );
    println!("      {:<w$} {}", "RMSE:", float_to_significant_digits(qm.rmse, 6));
    println!(
      "      {:<w$} {}",
      "Mass Error:",
      float_to_significant_digits(qm.mass_error, 6)
    );
    println!(
      "      {:<w$} {}",
      "Relative L2 Error:",
      float_to_significant_digits(qm.rel_l2_error, 6)
    );

    println!("    Peak Metrics:");
    let pm = &m.aggregate.domain_agreement.peak_metrics;
    println!(
      "      {:<w$} {}",
      "Location Error:",
      float_to_significant_digits(pm.location_error, 6)
    );
    println!(
      "      {:<w$} {}",
      "Value Error:",
      float_to_significant_digits(pm.value_error, 6)
    );

    println!("    Error Statistics:");
    let abs_err = &m.aggregate.domain_agreement.abs_error_stats;
    let rel_err = &m.aggregate.domain_agreement.rel_error_stats;
    println!(
      "      {:<w$} max={}, mean={}, std={}",
      "Absolute Error:",
      float_to_significant_digits(abs_err.max, 6),
      float_to_significant_digits(abs_err.mean, 6),
      float_to_significant_digits(abs_err.std, 6)
    );
    println!(
      "      {:<w$} max={}, mean={}, mape={}%",
      "Relative Error:",
      float_to_significant_digits(rel_err.max, 6),
      float_to_significant_digits(rel_err.mean, 6),
      float_to_significant_digits(rel_err.mape, 3)
    );

    println!("    Pointwise Metrics:");
    let pw = &m.pointwise;
    println!(
      "      {:<w$} {} (dx={})",
      "Total Points:",
      pw.total_points,
      float_to_significant_digits(pw.dx, 6)
    );
    let pws = &pw.errors.summary;
    println!(
      "      {:<w$} max={}, mean={}, std={}",
      "Pointwise Abs Error:",
      float_to_significant_digits(pws.abs_max, 6),
      float_to_significant_digits(pws.abs_mean, 6),
      float_to_significant_digits(pws.abs_std, 6)
    );
    println!(
      "      {:<w$} max={}, mean={}, median={}",
      "Pointwise Rel Error:",
      float_to_significant_digits(pws.rel_max, 6),
      float_to_significant_digits(pws.rel_mean, 6),
      float_to_significant_digits(pws.rel_median, 6)
    );
    println!(
      "      {:<w$} {}",
      "Signed Bias:",
      float_to_significant_digits(pws.signed_bias, 6)
    );
    println!(
      "      {:<w$} {}",
      "Log Max Error:",
      float_to_significant_digits(pws.log_max, 6)
    );

    println!("    Structural Metrics:");
    let st = &pw.structural.summary;
    println!(
      "      {:<w$} max={}, mean={}",
      "1st Derivative Error:",
      float_to_significant_digits(st.d1_max, 6),
      float_to_significant_digits(st.d1_mean, 6)
    );
    println!(
      "      {:<w$} max={}, mean={}",
      "2nd Derivative Error:",
      float_to_significant_digits(st.d2_max, 6),
      float_to_significant_digits(st.d2_mean, 6)
    );
    println!(
      "      {:<w$} {}",
      "Symmetry Max:",
      float_to_significant_digits(st.symmetry_max, 6)
    );
    println!(
      "      {:<w$} {}",
      "Monotonicity Violations:", st.monotonicity_violation_count
    );

    println!("    Tolerance Compliance:");
    let tol = &pw.tolerance.summary;
    let tolerance_labels = ["strict", "moderate", "loose"];
    for (label, fraction) in izip!(tolerance_labels.iter(), tol.pass_fractions.iter()) {
      println!("      {:<w$} {:.2}%", format!("Pass rate ({label}):"), fraction * 100.0);
    }
    println!("      {:<w$} {}", "Support Mismatches:", tol.support_mismatch_count);

    println!("    Spatial Metrics:");
    let peak_region = &m.spatial.regional.summary.peak_region;
    let tail_region = &m.spatial.regional.summary.tail_region;
    println!(
      "      {:<w$} mean={}, max={}",
      "Peak Region Error:",
      float_to_significant_digits(peak_region.mean_error, 6),
      float_to_significant_digits(peak_region.max_error, 6)
    );
    println!(
      "      {:<w$} mean={}, max={}",
      "Tail Region Error:",
      float_to_significant_digits(tail_region.mean_error, 6),
      float_to_significant_digits(tail_region.max_error, 6)
    );

    let cum = &m.spatial.cumulative.summary;
    println!(
      "      {:<w$} final={}, max_abs={}",
      "Cumulative Error:",
      float_to_significant_digits(cum.final_value, 6),
      float_to_significant_digits(cum.max_abs, 6)
    );

    let win = &m.spatial.windowed.summary;
    println!(
      "      {:<w$} rms_max={}, max_max={}",
      "Windowed Error:",
      float_to_significant_digits(win.sliding_rms_max, 6),
      float_to_significant_digits(win.sliding_max_max, 6)
    );

    println!("    Performance:");
    println!(
      "      {:<w$} {}ms",
      "Execution Time:",
      float_to_significant_digits(result.execution_time_ms, 4)
    );
    println!();
  }

  /// Print error summary section
  pub fn print_error_summary<T: TestCase>(failures: &[&TestFailure<T>]) -> Result<(), Report> {
    if failures.is_empty() {
      return Ok(());
    }

    println!("\n\n## Errors\n");
    for failure in failures {
      let test_case_json = json_write_str(&failure.test_case, JsonPretty(true))?;
      println!("❌ **ERROR**: {} + {}\n", failure.test_case.name(), failure.algorithm);
      println!("```");
      println!("{}", failure.error);
      println!("```\n");
      println!("**Test case:**");
      println!("```json\n{test_case_json}\n```\n");
    }
    Ok(())
  }

  /// Print comprehensive summary to console
  pub fn print_summary<T: TestCase>(summary: &TestSummary, outcomes: &[TestRunOutcome<T>]) {
    println!("\n## Results\n");

    Self::print_overall_statistics(summary);
    Self::print_per_test_case_comparison(outcomes);
    Self::print_unified_metrics_table(summary, outcomes);
  }

  /// Print overall statistics section
  fn print_overall_statistics(summary: &TestSummary) {
    let total_tests = summary.total_tests;
    let total_successes = summary.total_successes;
    let total_failures = summary.total_failures;
    let total_algorithms = summary.total_algorithms;
    let execution_time_total_ms = summary.execution_time_total_ms;
    let test_suite = &summary.test_suite_name;

    println!("- Test suite: {test_suite}");
    println!("- Total tests run: {total_tests}");
    println!("- Successful runs: {total_successes}");
    println!("- Failed runs: {total_failures}");
    println!("- Total algorithms: {total_algorithms}");
    println!("- Total execution time: {execution_time_total_ms:.1}ms\n");
  }
}
