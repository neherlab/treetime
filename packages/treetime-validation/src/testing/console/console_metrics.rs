use crate::testing::framework::results::TestResult;
use crate::testing::framework::summary::TestSummary;
use crate::testing::framework::test_case::TestCase;
use ordered_float::OrderedFloat;
use std::collections::BTreeMap;
use treetime_utils::float_fmt::float_to_significant_digits;
use treetime_utils::iterator::mean_by_key::MeanByKey;

use crate::testing::console::console::ValidationConsole;

#[allow(clippy::multiple_inherent_impl)]
impl ValidationConsole {
  /// Compute all metrics for the table
  pub(crate) fn compute_all_metrics<T: TestCase>(
    summary: &TestSummary,
    grouped_by_algorithm: &BTreeMap<String, Vec<&TestResult<T>>>,
  ) -> BTreeMap<String, BTreeMap<&'static str, String>> {
    let mut all_metrics = BTreeMap::new();

    for (algorithm, algorithm_results) in grouped_by_algorithm {
      let algo_summary = summary
        .algorithm_summaries
        .iter()
        .find(|s| s.algorithm_name == *algorithm)
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
        (
          "exec_time",
          float_to_significant_digits(algo_summary.execution_time_total_ms, 3),
        ),
        ("success_rate", format!("{:.1}%", algo_summary.success_rate * 100.0)),
        ("error_count", format!("{}", algo_summary.error_failures)),
        (
          "r2_error_ppm_max",
          float_to_significant_digits((1.0 - algo_summary.r2_min) * 1_000_000.0, 3),
        ),
        (
          "r2_error_ppm_mean",
          float_to_significant_digits((1.0 - algo_summary.r2_mean) * 1_000_000.0, 3),
        ),
        (
          "r2_error_ppm_min",
          float_to_significant_digits((1.0 - algo_summary.r2_max) * 1_000_000.0, 3),
        ),
        (
          "correlation_error_ppm",
          float_to_significant_digits((1.0 - correlation_mean) * 1_000_000.0, 3),
        ),
        ("rmse_max", float_to_significant_digits(rmse_max, 3)),
        ("mass_error_max", float_to_significant_digits(mass_error_max, 3)),
        ("rel_l2_error_max", float_to_significant_digits(rel_l2_error_max, 3)),
        (
          "agg_abs_max",
          float_to_significant_digits(algo_summary.max_abs_error_overall, 3),
        ),
        ("agg_abs_mean", float_to_significant_digits(agg_abs_mean, 3)),
        (
          "agg_rel_max",
          float_to_significant_digits(algo_summary.max_rel_error_overall, 3),
        ),
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

      all_metrics.insert(algorithm.clone(), metrics.into_iter().collect::<BTreeMap<_, _>>());
    }

    all_metrics
  }
}
