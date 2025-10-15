use crate::distribution::reference::convolution_test::algorithms::ConvolutionAlgorithm;
use crate::distribution::reference::convolution_test::framework::{
  TestCase, TestFailure, TestResult, TestRunOutcome, TestSummary,
};
use crate::io::json::{JsonPretty, json_write_str};
use crate::utils::float_fmt::float_to_significant_digits;
use crate::utils::iterator::mean_by_key::MeanByKey;
use eyre::Report;
use itertools::Itertools;
use ordered_float::OrderedFloat;

/// Console display functionality for convolution test framework
pub struct ConvolutionTestConsole;

impl ConvolutionTestConsole {
  /// Print test framework header
  pub fn print_header(function_type: &str, test_cases_count: usize, algorithms_count: usize) {
    let total_tests = test_cases_count * algorithms_count;
    println!(
      "=== {} Convolution Test Framework ===",
      function_type.to_uppercase()
    );
    println!(
      "Running {test_cases_count} test cases with {algorithms_count} algorithms ({total_tests} total)\n",
    );
  }

  /// Print progress table header
  pub fn print_progress_table_header() {
    println!(
      "| {:^12} | {:^2} | {:^12} | {:^30} | {:^8} | {:^10} | {:^8} | {:^8} |",
      "Progress", "S", "Algorithm", "Test Case", "Time, ms", "R² err ppm", "RMSE", "Corr err"
    );
    println!(
      "|{:-<14}|{:-^4}|{:-<14}|{:-<32}|{:->10}|{:->12}|{:->10}|{:->10}|",
      "", "", "", "", "", "", "", ""
    );
  }

  /// Print success row in progress table
  pub fn print_success_row<T: TestCase>(
    result: &TestResult<T>,
    completed_count: usize,
    total_tests: usize,
  ) {
    let elapsed_ms = result.execution_time_ms;
    let r_squared = result.metrics.aggregate.domain_agreement.quality_metrics.r_squared;
    let r2_error_ppm = (1.0 - r_squared) * 1_000_000.0;
    let rmse = result.metrics.aggregate.domain_agreement.quality_metrics.rmse;
    let correlation = result.metrics.aggregate.domain_agreement.quality_metrics.correlation;
    let correlation_error_ppm = (1.0 - correlation) * 1_000_000.0;
    let progress = format!("[{completed_count}/{total_tests}]");
    println!(
      "| {:<12} | {:^1} | {:<12} | {:<30} | {:>8.1} | {:>10} | {:>8} | {:>8} |",
      progress,
      "✅",
      format!("{}", result.algorithm),
      result.test_case.name(),
      elapsed_ms,
      float_to_significant_digits(r2_error_ppm, 3),
      float_to_significant_digits(rmse, 3),
      float_to_significant_digits(correlation_error_ppm, 3)
    );
  }

  /// Print failure row in progress table
  pub fn print_failure_row<T: TestCase>(
    test_case: &T,
    algorithm: ConvolutionAlgorithm,
    elapsed_ms: f64,
    completed_count: usize,
    total_tests: usize,
  ) {
    let progress = format!("[{completed_count}/{total_tests}]");
    println!(
      "| {:<12} | {:^1} | {:<12} | {:<30} | {:>8.1} | {:>10} | {:>8} | {:>8} |",
      progress,
      "❌",
      format!("{algorithm}"),
      test_case.name(),
      elapsed_ms,
      "FAILED",
      "FAILED",
      "FAILED"
    );
  }

  /// Print error summary section
  pub fn print_error_summary<T: TestCase>(failures: &[&TestFailure<T>]) -> Result<(), Report> {
    if failures.is_empty() {
      return Ok(());
    }

    println!("\n\n=== ERRORS ENCOUNTERED ===");
    for failure in failures {
      let test_case_json = json_write_str(&failure.test_case, JsonPretty(true))?;
      println!(
        "\n❌ ERROR: {} + {}: {}\n{test_case_json}",
        failure.test_case.name(),
        failure.algorithm,
        failure.error
      );
    }
    println!();
    Ok(())
  }

  /// Print comprehensive summary to console
  pub fn print_summary<T: TestCase>(summary: &TestSummary, outcomes: &[TestRunOutcome<T>]) {
    println!(
      "=== {} CONVOLUTION TEST SUMMARY ===\n",
      summary.function_type.to_uppercase()
    );

    Self::print_overall_statistics(summary);
    Self::print_unified_metrics_table(summary, outcomes);
  }

  /// Print overall statistics section
  fn print_overall_statistics(summary: &TestSummary) {
    println!("Overall Statistics:");
    let function_type = &summary.function_type;
    let total_tests = summary.total_tests;
    let total_successes = summary.total_successes;
    let total_failures = summary.total_failures;
    let total_algorithms = summary.total_algorithms;
    let execution_time_total_ms = summary.execution_time_total_ms;
    let overall_assessment = &summary.overall_assessment;
    println!("  Function type: {function_type}");
    println!("  Total tests run: {total_tests}");
    println!("  Successful runs: {total_successes}");
    println!("  Failed runs: {total_failures}");
    println!("  Total algorithms: {total_algorithms}");
    println!("  Total execution time: {execution_time_total_ms:.1}ms");
    println!("  Assessment: {overall_assessment}\n");
  }

  /// Print comprehensive metrics table
  fn print_unified_metrics_table<T: TestCase>(summary: &TestSummary, outcomes: &[TestRunOutcome<T>]) {
    let successes: Vec<_> = outcomes
      .iter()
      .filter_map(|outcome| match outcome {
        TestRunOutcome::Success(result) => Some(result),
        TestRunOutcome::Failure(_) => None,
      })
      .collect();

    if successes.is_empty() {
      return;
    }

    println!("Comprehensive Metrics by Algorithm:");

    let grouped_by_algorithm = successes.into_iter().into_group_map_by(|result| result.algorithm);
    let algorithms: Vec<_> = grouped_by_algorithm.keys().sorted().copied().collect();
    let all_metrics = Self::compute_all_metrics(summary, &grouped_by_algorithm);

    Self::print_metrics_table_headers(&algorithms, &all_metrics);
    Self::print_metrics_table_rows(&algorithms, &all_metrics);
  }

  /// Compute all metrics for the table
  fn compute_all_metrics<T: TestCase>(
    summary: &TestSummary,
    grouped_by_algorithm: &std::collections::HashMap<ConvolutionAlgorithm, Vec<&TestResult<T>>>,
  ) -> std::collections::HashMap<ConvolutionAlgorithm, std::collections::HashMap<&'static str, String>> {
    let mut all_metrics = std::collections::HashMap::new();

    for (algorithm, algorithm_results) in grouped_by_algorithm {
      let algo_summary = summary
        .algorithm_summaries
        .iter()
        .find(|s| s.algorithm == *algorithm)
        .unwrap();

      let correlation_mean: f64 = algorithm_results
        .iter()
        .mean_by_key(|r| r.metrics.aggregate.domain_agreement.quality_metrics.correlation);

      let rmse_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.aggregate.domain_agreement.quality_metrics.rmse))
        .map_or(0.0, |r| r.metrics.aggregate.domain_agreement.quality_metrics.rmse);

      let mass_error_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.aggregate.domain_agreement.quality_metrics.mass_error.abs()))
        .map_or(0.0, |r| {
          r.metrics.aggregate.domain_agreement.quality_metrics.mass_error.abs()
        });

      let rel_l2_error_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.aggregate.domain_agreement.quality_metrics.rel_l2_error))
        .map_or(0.0, |r| {
          r.metrics.aggregate.domain_agreement.quality_metrics.rel_l2_error
        });

      let agg_abs_mean: f64 = algorithm_results
        .iter()
        .mean_by_key(|r| r.metrics.aggregate.domain_agreement.abs_error_stats.mean);

      let agg_rel_mean: f64 = algorithm_results
        .iter()
        .mean_by_key(|r| r.metrics.aggregate.domain_agreement.rel_error_stats.mean);

      let pw_abs_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.pointwise.errors.summary.abs_max))
        .map_or(0.0, |r| r.metrics.pointwise.errors.summary.abs_max);

      let pw_abs_mean = algorithm_results
        .iter()
        .map(|r| r.metrics.pointwise.errors.summary.abs_mean)
        .sum::<f64>()
        / algorithm_results.len() as f64;

      let pw_abs_std = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.pointwise.errors.summary.abs_std))
        .map_or(0.0, |r| r.metrics.pointwise.errors.summary.abs_std);

      let pw_rel_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.pointwise.errors.summary.rel_max))
        .map_or(0.0, |r| r.metrics.pointwise.errors.summary.rel_max);

      let pw_rel_mean = algorithm_results
        .iter()
        .map(|r| r.metrics.pointwise.errors.summary.rel_mean)
        .sum::<f64>()
        / algorithm_results.len() as f64;

      let pw_rel_median = algorithm_results
        .iter()
        .map(|r| r.metrics.pointwise.errors.summary.rel_median)
        .sum::<f64>()
        / algorithm_results.len() as f64;

      let pw_signed_bias_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.pointwise.errors.summary.signed_bias.abs()))
        .map_or(0.0, |r| r.metrics.pointwise.errors.summary.signed_bias.abs());

      let pw_log_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.pointwise.errors.summary.log_max))
        .map_or(0.0, |r| r.metrics.pointwise.errors.summary.log_max);

      let d1_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.pointwise.structural.summary.d1_max))
        .map_or(0.0, |r| r.metrics.pointwise.structural.summary.d1_max);

      let d1_mean = algorithm_results
        .iter()
        .map(|r| r.metrics.pointwise.structural.summary.d1_mean)
        .sum::<f64>()
        / algorithm_results.len() as f64;

      let d2_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.pointwise.structural.summary.d2_max))
        .map_or(0.0, |r| r.metrics.pointwise.structural.summary.d2_max);

      let d2_mean = algorithm_results
        .iter()
        .map(|r| r.metrics.pointwise.structural.summary.d2_mean)
        .sum::<f64>()
        / algorithm_results.len() as f64;

      let symmetry_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.pointwise.structural.summary.symmetry_max))
        .map_or(0.0, |r| r.metrics.pointwise.structural.summary.symmetry_max);

      let monotonicity_violations_max = algorithm_results
        .iter()
        .map(|r| r.metrics.pointwise.structural.summary.monotonicity_violation_count)
        .max()
        .unwrap_or(0);

      let peak_mean = algorithm_results
        .iter()
        .map(|r| r.metrics.spatial.regional.summary.peak_region.mean_error)
        .sum::<f64>()
        / algorithm_results.len() as f64;

      let peak_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.spatial.regional.summary.peak_region.max_error))
        .map_or(0.0, |r| r.metrics.spatial.regional.summary.peak_region.max_error);

      let tail_mean = algorithm_results
        .iter()
        .map(|r| r.metrics.spatial.regional.summary.tail_region.mean_error)
        .sum::<f64>()
        / algorithm_results.len() as f64;

      let tail_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.spatial.regional.summary.tail_region.max_error))
        .map_or(0.0, |r| r.metrics.spatial.regional.summary.tail_region.max_error);

      let cumulative_final = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.spatial.cumulative.summary.final_value.abs()))
        .map_or(0.0, |r| r.metrics.spatial.cumulative.summary.final_value.abs());

      let cumulative_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.spatial.cumulative.summary.max_abs))
        .map_or(0.0, |r| r.metrics.spatial.cumulative.summary.max_abs);

      let sliding_rms_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.spatial.windowed.summary.sliding_rms_max))
        .map_or(0.0, |r| r.metrics.spatial.windowed.summary.sliding_rms_max);

      let sliding_max_max = algorithm_results
        .iter()
        .max_by_key(|r| OrderedFloat(r.metrics.spatial.windowed.summary.sliding_max_max))
        .map_or(0.0, |r| r.metrics.spatial.windowed.summary.sliding_max_max);

      let tolerance_strict = algorithm_results
        .iter()
        .min_by_key(|r| OrderedFloat(r.metrics.pointwise.tolerance.summary.pass_fractions[0] * 100.0))
        .map_or(0.0, |r| r.metrics.pointwise.tolerance.summary.pass_fractions[0] * 100.0);

      let tolerance_moderate = algorithm_results
        .iter()
        .min_by_key(|r| OrderedFloat(r.metrics.pointwise.tolerance.summary.pass_fractions[1] * 100.0))
        .map_or(0.0, |r| r.metrics.pointwise.tolerance.summary.pass_fractions[1] * 100.0);

      let tolerance_loose = algorithm_results
        .iter()
        .min_by_key(|r| OrderedFloat(r.metrics.pointwise.tolerance.summary.pass_fractions[2] * 100.0))
        .map_or(0.0, |r| r.metrics.pointwise.tolerance.summary.pass_fractions[2] * 100.0);

      let support_mismatch_max = algorithm_results
        .iter()
        .map(|r| r.metrics.pointwise.tolerance.summary.support_mismatch_count)
        .max()
        .unwrap_or(0);

      let metrics = vec![
        ("test_count", format!("{}", algo_summary.test_cases_count)),
        ("exec_time", float_to_significant_digits(algo_summary.execution_time_total_ms, 3)),
        ("success_rate", format!("{:.1}%", algo_summary.success_rate * 100.0)),
        ("error_count", format!("{}", algo_summary.error_failures)),
        ("r2_error_ppm_max", float_to_significant_digits((1.0 - algo_summary.r2_min) * 1_000_000.0, 3)),
        ("r2_error_ppm_mean", float_to_significant_digits((1.0 - algo_summary.r2_mean) * 1_000_000.0, 3)),
        ("r2_error_ppm_min", float_to_significant_digits((1.0 - algo_summary.r2_max) * 1_000_000.0, 3)),
        ("correlation_error_ppm", float_to_significant_digits((1.0 - correlation_mean) * 1_000_000.0, 3)),
        ("rmse_max", float_to_significant_digits(rmse_max, 3)),
        ("mass_error_max", float_to_significant_digits(mass_error_max, 3)),
        ("rel_l2_error_max", float_to_significant_digits(rel_l2_error_max, 3)),
        (
          "agg_abs_max",
          float_to_significant_digits(algo_summary.max_abs_error_overall, 3),
        ),
        ("agg_abs_mean", float_to_significant_digits(agg_abs_mean, 3)),
        ("agg_rel_max", float_to_significant_digits(algo_summary.max_rel_error_overall, 3)),
        ("agg_rel_mean", float_to_significant_digits(agg_rel_mean, 3)),
        ("pw_abs_max", float_to_significant_digits(pw_abs_max, 3)),
        ("pw_abs_mean", float_to_significant_digits(pw_abs_mean, 3)),
        ("pw_abs_std", float_to_significant_digits(pw_abs_std, 3)),
        ("pw_rel_max", float_to_significant_digits(pw_rel_max, 3)),
        ("pw_rel_mean", float_to_significant_digits(pw_rel_mean, 3)),
        ("pw_rel_median", float_to_significant_digits(pw_rel_median, 3)),
        ("pw_signed_bias_max", float_to_significant_digits(pw_signed_bias_max, 3)),
        ("pw_log_max", float_to_significant_digits(pw_log_max, 3)),
        ("d1_max", float_to_significant_digits(d1_max, 3)),
        ("d1_mean", float_to_significant_digits(d1_mean, 3)),
        ("d2_max", float_to_significant_digits(d2_max, 3)),
        ("d2_mean", float_to_significant_digits(d2_mean, 3)),
        ("symmetry_max", float_to_significant_digits(symmetry_max, 3)),
        ("monotonicity_violations", format!("{monotonicity_violations_max}")),
        ("peak_mean", float_to_significant_digits(peak_mean, 3)),
        ("peak_max", float_to_significant_digits(peak_max, 3)),
        ("tail_mean", float_to_significant_digits(tail_mean, 3)),
        ("tail_max", float_to_significant_digits(tail_max, 3)),
        ("cumulative_final", float_to_significant_digits(cumulative_final, 3)),
        ("cumulative_max", float_to_significant_digits(cumulative_max, 3)),
        ("sliding_rms_max", float_to_significant_digits(sliding_rms_max, 3)),
        ("sliding_max_max", float_to_significant_digits(sliding_max_max, 3)),
        ("tolerance_strict", format!("{tolerance_strict:.1}%")),
        ("tolerance_moderate", format!("{tolerance_moderate:.1}%")),
        ("tolerance_loose", format!("{tolerance_loose:.1}%")),
        ("support_mismatch", format!("{support_mismatch_max}")),
      ];

      all_metrics.insert(
        *algorithm,
        metrics.into_iter().collect::<std::collections::HashMap<_, _>>(),
      );
    }

    all_metrics
  }

  /// Print metrics table headers
  fn print_metrics_table_headers(
    algorithms: &[ConvolutionAlgorithm],
    all_metrics: &std::collections::HashMap<ConvolutionAlgorithm, std::collections::HashMap<&'static str, String>>,
  ) {
    let metric_col_width = "Moderate tolerance pass (min%)".len();
    let mut algo_col_widths = std::collections::HashMap::new();

    for algo in algorithms {
      let name_w = algo.to_string().len();
      let max_value_w = all_metrics[algo].values().map(|s| s.len()).max().unwrap_or(0);
      algo_col_widths.insert(*algo, name_w.max(max_value_w));
    }

    print!("| {:^width$} |", "Metric", width = metric_col_width);
    for algo in algorithms {
      print!(" {:^width$} |", algo.to_string(), width = algo_col_widths[algo]);
    }
    println!();

    print!("|{:-^width$}|", "", width = metric_col_width + 2);
    for algo in algorithms {
      print!("{:-^width$}|", "", width = algo_col_widths[algo] + 2);
    }
    println!();
  }

  /// Print metrics table rows
  fn print_metrics_table_rows(
    algorithms: &[ConvolutionAlgorithm],
    all_metrics: &std::collections::HashMap<ConvolutionAlgorithm, std::collections::HashMap<&'static str, String>>,
  ) {
    let metric_col_width = "Moderate tolerance pass (min%)".len();
    let mut algo_col_widths = std::collections::HashMap::new();

    for algo in algorithms {
      let name_w = algo.to_string().len();
      let max_value_w = all_metrics[algo].values().map(|s| s.len()).max().unwrap_or(0);
      algo_col_widths.insert(*algo, name_w.max(max_value_w));
    }

    let categories = vec![
      (
        "PERFORMANCE",
        vec![
          ("Test count", "test_count"),
          ("Execution time (ms)", "exec_time"),
          ("Success rate", "success_rate"),
          ("Error count", "error_count"),
        ],
      ),
      (
        "QUALITY METRICS",
        vec![
          ("R² error max (ppm)", "r2_error_ppm_max"),
          ("R² error min (ppm)", "r2_error_ppm_min"),
          ("R² error mean (ppm)", "r2_error_ppm_mean"),
          ("Correlation error (ppm)", "correlation_error_ppm"),
          ("RMSE (max)", "rmse_max"),
          ("Mass conservation error (max)", "mass_error_max"),
          ("Relative L2 error (max)", "rel_l2_error_max"),
        ],
      ),
      (
        "AGGREGATE ERRORS",
        vec![
          ("Absolute max", "agg_abs_max"),
          ("Absolute mean", "agg_abs_mean"),
          ("Relative max", "agg_rel_max"),
          ("Relative mean", "agg_rel_mean"),
        ],
      ),
      (
        "POINTWISE ERRORS",
        vec![
          ("Absolute max", "pw_abs_max"),
          ("Absolute mean", "pw_abs_mean"),
          ("Absolute std dev", "pw_abs_std"),
          ("Relative max", "pw_rel_max"),
          ("Relative mean", "pw_rel_mean"),
          ("Relative median", "pw_rel_median"),
          ("Signed bias (max)", "pw_signed_bias_max"),
          ("Logarithmic max", "pw_log_max"),
        ],
      ),
      (
        "STRUCTURAL",
        vec![
          ("1st derivative error (max)", "d1_max"),
          ("1st derivative error (mean)", "d1_mean"),
          ("2nd derivative error (max)", "d2_max"),
          ("2nd derivative error (mean)", "d2_mean"),
          ("Symmetry residual (max)", "symmetry_max"),
          ("Monotonicity violations (max)", "monotonicity_violations"),
        ],
      ),
      (
        "REGIONAL",
        vec![
          ("Peak region error (mean)", "peak_mean"),
          ("Peak region error (max)", "peak_max"),
          ("Tail region error (mean)", "tail_mean"),
          ("Tail region error (max)", "tail_max"),
        ],
      ),
      (
        "CUMULATIVE & WINDOWED",
        vec![
          ("Cumulative error (final)", "cumulative_final"),
          ("Cumulative error (max abs)", "cumulative_max"),
          ("Sliding RMS error (max)", "sliding_rms_max"),
          ("Sliding max error (max)", "sliding_max_max"),
        ],
      ),
      (
        "TOLERANCE COMPLIANCE",
        vec![
          ("Strict tolerance pass (min%)", "tolerance_strict"),
          ("Moderate tolerance pass (min%)", "tolerance_moderate"),
          ("Loose tolerance pass (min%)", "tolerance_loose"),
          ("Support mismatches (max)", "support_mismatch"),
        ],
      ),
    ];

    for (category, metrics) in categories {
      print!("| {:width$} |", "", width = metric_col_width);
      for algo in algorithms {
        print!(" {:width$} |", "", width = algo_col_widths[algo]);
      }
      println!();

      print!("| {category:<metric_col_width$} |");
      for algo in algorithms {
        print!(" {:width$} |", "", width = algo_col_widths[algo]);
      }
      println!();

      for (metric_name, metric_key) in metrics {
        print!("| {metric_name:<metric_col_width$} |");
        for algo in algorithms {
          let value = &all_metrics[algo][metric_key];
          print!(" {:>width$} |", value, width = algo_col_widths[algo]);
        }
        println!();
      }
    }

    println!();
  }
}
