use crate::testing::framework::results::TestRunOutcome;
use crate::testing::framework::summary::TestSummary;
use crate::testing::framework::test_case::TestCase;
use itertools::Itertools;
use std::collections::BTreeMap;
use treetime_utils::float_fmt::float_to_significant_digits;

use crate::testing::console::console::ValidationConsole;

#[allow(clippy::multiple_inherent_impl)]
impl ValidationConsole {
  /// Print per-test-case algorithm comparison table
  pub(crate) fn print_per_test_case_comparison<T: TestCase>(outcomes: &[TestRunOutcome<T>]) {
    let successes: Vec<_> = outcomes
      .iter()
      .filter_map(|outcome| match outcome {
        TestRunOutcome::Success(result) => Some(result.as_ref()),
        TestRunOutcome::Failure(_) => None,
      })
      .collect();

    if successes.is_empty() {
      return;
    }

    let algorithms: Vec<_> = successes
      .iter()
      .map(|r| r.algorithm.clone())
      .sorted()
      .dedup()
      .collect_vec();

    if algorithms.len() < 2 {
      return;
    }

    let grouped_by_test_case: BTreeMap<_, Vec<_>> = successes
      .into_iter()
      .into_group_map_by(|result| result.test_case.name().to_owned())
      .into_iter()
      .collect();

    println!("### By Test Case\n");

    let test_case_col_width = grouped_by_test_case
      .keys()
      .map(|name| name.len())
      .max()
      .map_or(20, |w| w.max(20));
    let algo_col_width = algorithms.iter().map(|a| a.len()).max().map_or(12, |w| w.max(12));

    print!("| {:^test_case_col_width$} |", "Test Case");
    for algo in &algorithms {
      print!(" {algo:^algo_col_width$} |");
    }
    println!();

    print!("|{:-^width$}|", "", width = test_case_col_width + 2);
    for _ in &algorithms {
      print!("{:-^width$}|", "", width = algo_col_width + 2);
    }
    println!();

    for (test_case_name, results) in &grouped_by_test_case {
      let results_by_algo: BTreeMap<_, _> = results.iter().map(|r| (r.algorithm.clone(), *r)).collect();

      print!("| {test_case_name:<test_case_col_width$} |");

      for algo in &algorithms {
        if let Some(result) = results_by_algo.get(algo) {
          let r2 = result.metrics.aggregate.domain_agreement.quality_metrics.r_squared;
          let r2_error_ppm = (1.0 - r2) * 1_000_000.0;
          print!(
            " {:>width$} |",
            float_to_significant_digits(r2_error_ppm, 3),
            width = algo_col_width
          );
        } else {
          print!(" {:^width$} |", "N/A", width = algo_col_width);
        }
      }

      println!();
    }

    println!("\n(Values shown: R^2 error in ppm - lower is better)\n");
  }

  /// Print comprehensive metrics table
  pub(crate) fn print_unified_metrics_table<T: TestCase>(summary: &TestSummary, outcomes: &[TestRunOutcome<T>]) {
    let successes: Vec<_> = outcomes
      .iter()
      .filter_map(|outcome| match outcome {
        TestRunOutcome::Success(result) => Some(result.as_ref()),
        TestRunOutcome::Failure(_) => None,
      })
      .collect();

    if successes.is_empty() {
      return;
    }

    println!("### By Algorithm\n");

    let grouped_by_algorithm: BTreeMap<_, _> = successes
      .into_iter()
      .into_group_map_by(|result| result.algorithm.clone())
      .into_iter()
      .collect();
    let algorithms: Vec<_> = grouped_by_algorithm.keys().cloned().collect();
    let all_metrics = Self::compute_all_metrics(summary, &grouped_by_algorithm);

    Self::print_metrics_table_headers(&algorithms, &all_metrics);
    Self::print_metrics_table_rows(&algorithms, &all_metrics);
  }

  /// Print metrics table headers
  fn print_metrics_table_headers(
    algorithms: &[String],
    all_metrics: &BTreeMap<String, BTreeMap<&'static str, String>>,
  ) {
    let metric_col_width = "Moderate tolerance pass (min%)".len();
    let mut algo_col_widths = BTreeMap::new();

    for algo in algorithms {
      let name_w = algo.len();
      let max_value_w = all_metrics[algo].values().map(|s| s.len()).max().unwrap_or(0);
      algo_col_widths.insert(algo.clone(), name_w.max(max_value_w));
    }

    print!("| {:^width$} |", "Metric", width = metric_col_width);
    for algo in algorithms {
      print!(" {:^width$} |", algo, width = algo_col_widths[algo]);
    }
    println!();

    print!("|{:-^width$}|", "", width = metric_col_width + 2);
    for algo in algorithms {
      print!("{:-^width$}|", "", width = algo_col_widths[algo] + 2);
    }
    println!();
  }

  /// Print metrics table rows
  fn print_metrics_table_rows(algorithms: &[String], all_metrics: &BTreeMap<String, BTreeMap<&'static str, String>>) {
    let metric_col_width = "Moderate tolerance pass (min%)".len();
    let mut algo_col_widths = BTreeMap::new();

    for algo in algorithms {
      let name_w = algo.len();
      let max_value_w = all_metrics[algo].values().map(|s| s.len()).max().unwrap_or(0);
      algo_col_widths.insert(algo.clone(), name_w.max(max_value_w));
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
      print!(
        "| **{category}** {:<width$} |",
        "",
        width = metric_col_width - category.len() - 5
      );
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
