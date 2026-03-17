#[cfg(test)]
mod tests {
  use crate::commands::clock::clock_model::{ClockModel, ClockModelStats, RegressionStats};
  use crate::commands::timetree::output::confidence::{
    combine_confidence, date_uncertainty_due_to_rate, determine_rate_std, extract_confidence_intervals,
    quantile_to_zscore,
  };
  use crate::representation::partition::timetree::GraphTimetree;
  use crate::representation::payload::timetree::NodeTimetree;
  use approx::assert_relative_eq;
  use ndarray::array;
  use rstest::rstest;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_graph::node::Named;

  fn make_node(name: Option<&str>, time: Option<f64>, dist: Option<Arc<Distribution>>) -> NodeTimetree {
    let mut node = NodeTimetree::default();
    node.base.set_name(name);
    node.time = time;
    node.time_distribution = dist;
    node
  }

  fn make_node_with_rate_dates(
    name: &str,
    time: f64,
    dist: Option<Arc<Distribution>>,
    rate_dates: [f64; 3],
  ) -> NodeTimetree {
    let mut node = make_node(Some(name), Some(time), dist);
    node.rate_susceptibility_dates = Some(rate_dates);
    node
  }

  // --- extract_confidence_intervals ---

  #[test]
  fn test_extract_confidence_intervals_skips_unnamed_nodes() {
    let mut graph = GraphTimetree::new();
    graph.add_node(make_node(None, Some(2020.0), None));
    graph.add_node(make_node(Some("named"), Some(2021.0), None));
    graph.build().unwrap();

    let intervals = extract_confidence_intervals(&graph);
    assert_eq!(intervals.len(), 1);
    assert_eq!(intervals[0].name, "named");
  }

  #[test]
  fn test_extract_confidence_intervals_skips_nodes_without_time() {
    let mut graph = GraphTimetree::new();
    graph.add_node(make_node(Some("no_time"), None, None));
    graph.add_node(make_node(Some("has_time"), Some(2021.0), None));
    graph.build().unwrap();

    let intervals = extract_confidence_intervals(&graph);
    assert_eq!(intervals.len(), 1);
    assert_eq!(intervals[0].name, "has_time");
  }

  #[test]
  fn test_extract_confidence_intervals_uses_date_as_fallback() {
    let mut graph = GraphTimetree::new();
    graph.add_node(make_node(Some("node_a"), Some(2020.5), None));
    graph.build().unwrap();

    let intervals = extract_confidence_intervals(&graph);
    assert_eq!(intervals.len(), 1);
    assert_relative_eq!(intervals[0].date, 2020.5);
    assert_relative_eq!(intervals[0].lower, 2020.5);
    assert_relative_eq!(intervals[0].upper, 2020.5);
  }

  #[test]
  fn test_extract_confidence_intervals_with_distribution() {
    let mut graph = GraphTimetree::new();
    let dist = Arc::new(Distribution::range((2019.0, 2021.0), 1.0));
    graph.add_node(make_node(Some("node_a"), Some(2020.0), Some(dist)));
    graph.build().unwrap();

    let intervals = extract_confidence_intervals(&graph);
    assert_eq!(intervals.len(), 1);
    assert_relative_eq!(intervals[0].date, 2020.0);
    // 90% CI from uniform [2019, 2021]: 0.05 * 2 + 2019 = 2019.1, 0.95 * 2 + 2019 = 2020.9
    assert_relative_eq!(intervals[0].lower, 2019.1, epsilon = 1e-10);
    assert_relative_eq!(intervals[0].upper, 2020.9, epsilon = 1e-10);
  }

  #[test]
  fn test_extract_confidence_intervals_sorted_by_name() {
    let mut graph = GraphTimetree::new();
    graph.add_node(make_node(Some("zebra"), Some(2020.0), None));
    graph.add_node(make_node(Some("alpha"), Some(2021.0), None));
    graph.add_node(make_node(Some("middle"), Some(2022.0), None));
    graph.build().unwrap();

    let intervals = extract_confidence_intervals(&graph);
    assert_eq!(intervals.len(), 3);
    assert_eq!(intervals[0].name, "alpha");
    assert_eq!(intervals[1].name, "middle");
    assert_eq!(intervals[2].name, "zebra");
  }

  #[test]
  fn test_extract_confidence_intervals_rate_only() {
    // Rate susceptibility data but no marginal distribution
    let mut graph = GraphTimetree::new();
    graph.add_node(make_node_with_rate_dates(
      "node_a",
      2010.0,
      None,
      [2009.0, 2010.0, 2011.0],
    ));
    graph.build().unwrap();

    let intervals = extract_confidence_intervals(&graph);
    assert_eq!(intervals.len(), 1);
    assert_relative_eq!(intervals[0].date, 2010.0);
    // z-score at 0.05 = -1.644854, at 0.95 = +1.644854
    // lower = 2010 + (-1.644854) * |2009 - 2010| = 2008.355146
    // upper = 2010 + 1.644854 * |2011 - 2010| = 2011.644854
    assert_relative_eq!(intervals[0].lower, 2008.355146, epsilon = 1e-4);
    assert_relative_eq!(intervals[0].upper, 2011.644854, epsilon = 1e-4);
  }

  #[test]
  fn test_extract_confidence_intervals_combined_wider_than_either() {
    // Both marginal distribution and rate susceptibility data present.
    // The quadrature combination must be wider than either source alone.
    let mut graph = GraphTimetree::new();
    let dist = Arc::new(Distribution::range((2008.0, 2012.0), 1.0));
    graph.add_node(make_node_with_rate_dates(
      "node_a",
      2010.0,
      Some(dist),
      [2009.0, 2010.0, 2011.0],
    ));
    graph.build().unwrap();

    let intervals = extract_confidence_intervals(&graph);
    assert_eq!(intervals.len(), 1);
    // Mutation CI from uniform [2008, 2012]: 90% = [2008.2, 2011.8]
    // Rate CI at 90%: [2008.355, 2011.645]
    // Combined via quadrature: strictly wider than either
    assert!(
      intervals[0].lower < 2008.2,
      "combined lower should be below mutation-only lower"
    );
    assert!(
      intervals[0].upper > 2011.8,
      "combined upper should be above mutation-only upper"
    );
  }

  // --- combine_confidence ---

  #[test]
  fn test_combine_confidence_no_contributions() {
    let result = combine_confidence(10.0, (5.0, 15.0), None, None);
    assert_relative_eq!(result.0, 5.0);
    assert_relative_eq!(result.1, 15.0);
  }

  #[test]
  fn test_combine_confidence_single_contribution() {
    let result = combine_confidence(10.0, (0.0, 20.0), Some((8.0, 12.0)), None);
    assert_relative_eq!(result.0, 8.0);
    assert_relative_eq!(result.1, 12.0);

    let result = combine_confidence(10.0, (0.0, 20.0), None, Some((7.0, 13.0)));
    assert_relative_eq!(result.0, 7.0);
    assert_relative_eq!(result.1, 13.0);
  }

  #[test]
  fn test_combine_confidence_quadrature() {
    // c1: (8, 12) -> deviations of 2 from center 10
    // c2: (7, 13) -> deviations of 3 from center 10
    // Combined: sqrt(2^2 + 3^2) = sqrt(13) = 3.606
    let result = combine_confidence(10.0, (0.0, 20.0), Some((8.0, 12.0)), Some((7.0, 13.0)));
    let expected_dev = (4.0_f64 + 9.0).sqrt();
    assert_relative_eq!(result.0, 10.0 - expected_dev, epsilon = 1e-10);
    assert_relative_eq!(result.1, 10.0 + expected_dev, epsilon = 1e-10);
  }

  #[test]
  fn test_combine_confidence_clipped_to_limits() {
    // Large contributions that exceed limits
    let result = combine_confidence(10.0, (5.0, 15.0), Some((0.0, 20.0)), Some((0.0, 20.0)));
    // Quadrature would give sqrt(100 + 100) = 14.14 deviation
    // But limits clip to (5, 15)
    assert_relative_eq!(result.0, 5.0);
    assert_relative_eq!(result.1, 15.0);
  }

  // --- quantile_to_zscore ---
  // Probit function: z = sqrt(2) * erf_inv(2p - 1)
  // Standard values: p=0.025 -> z=-1.959964, p=0.5 -> z=0, p=0.975 -> z=1.959964

  #[rustfmt::skip]
  #[rstest]
  #[case::lower_2_5pct(0.025,  -1.959964)]
  #[case::lower_5pct(  0.05,   -1.644854)]
  #[case::median(       0.5,    0.0      )]
  #[case::upper_95pct(  0.95,   1.644854 )]
  #[case::upper_97_5pct(0.975,  1.959964 )]
  #[case::boundary_zero(0.0,    0.0      )]
  #[case::boundary_one( 1.0,    0.0      )]
  #[trace]
  fn test_quantile_to_zscore(#[case] p: f64, #[case] expected: f64) {
    let z = quantile_to_zscore(p);
    assert_relative_eq!(z, expected, epsilon = 1e-4);
  }

  // --- date_uncertainty_due_to_rate ---
  // Converts [lower_date, center_date, upper_date] + quantile interval to CI.
  // ci_lower = center + z(p_lo) * |lower - center|
  // ci_upper = center + z(p_hi) * |upper - center|

  #[rustfmt::skip]
  #[rstest]
  #[case::symmetric_1sigma(   [9.0, 10.0, 11.0],  (0.025, 0.975), (8.040036, 11.959964))]
  #[case::symmetric_narrow(   [9.5, 10.0, 10.5],  (0.025, 0.975), (9.020018, 10.979982))]
  #[case::asymmetric(         [8.0, 10.0, 11.0],  (0.025, 0.975), (6.080072, 11.959964))]
  #[case::boundary_quantiles( [9.0, 10.0, 11.0],  (0.0,   1.0),   (10.0,     10.0)     )]
  #[case::equal_dates(        [10.0, 10.0, 10.0], (0.025, 0.975), (10.0,     10.0)     )]
  #[trace]
  fn test_date_uncertainty_due_to_rate(
    #[case] dates: [f64; 3],
    #[case] interval: (f64, f64),
    #[case] (expected_lower, expected_upper): (f64, f64),
  ) {
    let (lower, upper) = date_uncertainty_due_to_rate(dates, interval);
    assert_relative_eq!(lower, expected_lower, epsilon = 1e-4);
    assert_relative_eq!(upper, expected_upper, epsilon = 1e-4);
  }

  // --- determine_rate_std ---

  #[test]
  fn test_determine_rate_std_explicit_clock_std_dev() {
    let clock_model = ClockModel::for_testing(0.003, 0.0);
    let result = determine_rate_std(Some(0.001), false, &clock_model).unwrap();
    assert_relative_eq!(result.unwrap(), 0.001);
  }

  #[test]
  fn test_determine_rate_std_rejects_negative() {
    let clock_model = ClockModel::for_testing(0.003, 0.0);
    let result = determine_rate_std(Some(-0.001), false, &clock_model);
    let _unused = result.unwrap_err();
  }

  #[test]
  fn test_determine_rate_std_rejects_zero() {
    let clock_model = ClockModel::for_testing(0.003, 0.0);
    let result = determine_rate_std(Some(0.0), false, &clock_model);
    let _unused = result.unwrap_err();
  }

  #[test]
  fn test_determine_rate_std_none_without_covariation() {
    let clock_model = ClockModel::for_testing(0.003, 0.0);
    let result = determine_rate_std(None, false, &clock_model).unwrap();
    assert!(result.is_none());
  }

  #[test]
  fn test_determine_rate_std_from_covariance_matrix() {
    // cov[0,0] = 1e-6, so rate_std = 1e-3
    let clock_model = ClockModel::for_testing_with_stats(
      0.003,
      0.0,
      ClockModelStats::Estimated(RegressionStats {
        chisq: 0.0,
        r_val: 0.9,
        hessian: array![[1.0, 0.0], [0.0, 1.0]],
        cov: array![[1e-6, 0.0], [0.0, 1.0]],
      }),
    );
    let result = determine_rate_std(None, true, &clock_model).unwrap();
    assert_relative_eq!(result.unwrap(), 1e-3, epsilon = 1e-10);
  }

  #[test]
  fn test_determine_rate_std_none_for_fixed_clock_with_covariation() {
    let clock_model = ClockModel::for_testing(0.003, 0.0);
    let result = determine_rate_std(None, true, &clock_model).unwrap();
    assert!(result.is_none());
  }
}
