#[cfg(test)]
mod tests {
  use crate::commands::timetree::output::confidence::extract_confidence_intervals;
  use crate::representation::partition::timetree::GraphTimetree;
  use crate::representation::payload::timetree::NodeTimetree;
  use approx::assert_relative_eq;
  use helpers::{make_node, make_node_with_rate_dates};
  use ndarray::Array1;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_graph::node::Named;

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

  #[test]
  fn test_extract_confidence_intervals_clamps_when_date_outside_rate_ci() {
    // When the final marginal pass date differs from the rate susceptibility
    // central date, the raw rate CI may not bracket the point estimate.
    // The postcondition clamp ensures lower <= date <= upper.
    let mut graph = GraphTimetree::new();
    // date = 2020.5 (final pass), rate susceptibility centered on 2020.0
    // with small variation [2019.9, 2020.0, 2020.1].
    // Rate CI at 90%: 2020.0 +/- 1.645 * 0.1 = [2019.836, 2020.164]
    // date = 2020.5 > 2020.164, so upper must be clamped to date.
    graph.add_node(make_node_with_rate_dates(
      "node_a",
      2020.5,
      None,
      [2019.9, 2020.0, 2020.1],
    ));
    graph.build().unwrap();

    let intervals = extract_confidence_intervals(&graph);
    assert_eq!(intervals.len(), 1);
    // Postcondition holds: lower <= date <= upper
    assert!(intervals[0].lower <= intervals[0].date);
    assert!(intervals[0].date <= intervals[0].upper);
    // Upper was clamped to date since raw rate CI upper (2020.164) < date (2020.5)
    assert_relative_eq!(intervals[0].upper, 2020.5);
    // Lower stays at raw rate CI lower (unclamped, already below date)
    assert_relative_eq!(intervals[0].lower, 2019.8355, epsilon = 1e-3);
  }

  #[test]
  fn test_extract_confidence_intervals_clamps_when_date_below_rate_ci() {
    // Mirror case: date below the raw rate CI lower bound.
    let mut graph = GraphTimetree::new();
    // date = 2019.5 (final pass), rate susceptibility centered on 2020.0
    // Rate CI at 90%: [2019.836, 2020.164]
    // date = 2019.5 < 2019.836, so lower must be clamped to date.
    graph.add_node(make_node_with_rate_dates(
      "node_a",
      2019.5,
      None,
      [2019.9, 2020.0, 2020.1],
    ));
    graph.build().unwrap();

    let intervals = extract_confidence_intervals(&graph);
    assert_eq!(intervals.len(), 1);
    assert!(intervals[0].lower <= intervals[0].date);
    assert!(intervals[0].date <= intervals[0].upper);
    // Lower was clamped to date since raw rate CI lower (2019.836) > date (2019.5)
    assert_relative_eq!(intervals[0].lower, 2019.5);
    // Upper stays at raw rate CI upper (unclamped, already above date)
    assert_relative_eq!(intervals[0].upper, 2020.1645, epsilon = 1e-3);
  }

  // v0 uses get_max_posterior_region(fraction=0.9): highest posterior density region,
  // the NARROWEST interval containing 90% probability mass.
  // For symmetric distributions, HPD equals equal-tailed CI.
  // For skewed distributions (nodes near tree boundaries), HPD is narrower and
  // centered on the peak.

  #[test]
  fn test_extract_confidence_intervals_skewed_distribution_hpd() {
    // Discretized exponential distribution: P(t) = exp(-t) on [0, 10].
    // Peak at t=0, long right tail.
    //
    // Analytical CDF: F(t) = 1 - exp(-t)
    // Equal-tailed 90% CI: [quantile(0.05), quantile(0.95)]
    //   = [-ln(0.95), -ln(0.05)] = [0.0513, 2.9957]
    //   width = 2.9444
    //
    // HPD 90% region: the shortest interval [0, h] such that F(h) - F(0) = 0.9
    //   F(h) = 0.9 => h = -ln(0.1) = 2.3026
    //   HPD = [0, 2.3026], width = 2.3026 (22% narrower)
    let n_points = 500;
    let x_min = 0.0;
    let dx = 10.0 / (n_points as f64 - 1.0);
    let y = Array1::from_shape_fn(n_points, |i| (-(x_min + i as f64 * dx)).exp());

    let dist_fn = treetime_distribution::DistributionFunction::from_start_dx_values(x_min, dx, y).unwrap();
    let dist = Distribution::Function(dist_fn);
    let peak_time = dist.likely_time().unwrap();

    let mut graph = GraphTimetree::new();
    let mut node = NodeTimetree::default();
    node.base.set_name(Some("skewed"));
    node.time = Some(peak_time);
    node.time_distribution = Some(Arc::new(dist));
    graph.add_node(node);
    graph.build().unwrap();

    let intervals = extract_confidence_intervals(&graph);
    assert_eq!(intervals.len(), 1);

    // v0 HPD bounds: [0, 2.3026] (narrowest 90% interval around peak)
    // Allow tolerance for grid discretization (dt ~ 0.02)
    let hpd_lower = 0.0;
    let hpd_upper = (0.1_f64).ln().abs(); // -ln(0.1) = 2.3026
    assert_relative_eq!(intervals[0].lower, hpd_lower, epsilon = dx);
    assert_relative_eq!(intervals[0].upper, hpd_upper, epsilon = dx);
  }

  mod helpers {
    use crate::representation::payload::timetree::NodeTimetree;
    use std::sync::Arc;
    use treetime_distribution::Distribution;
    use treetime_graph::node::Named;

    pub fn make_node(name: Option<&str>, time: Option<f64>, dist: Option<Arc<Distribution>>) -> NodeTimetree {
      let mut node = NodeTimetree::default();
      node.base.set_name(name);
      node.time = time;
      node.time_distribution = dist;
      node
    }

    pub fn make_node_with_rate_dates(
      name: &str,
      time: f64,
      dist: Option<Arc<Distribution>>,
      rate_dates: [f64; 3],
    ) -> NodeTimetree {
      let mut node = make_node(Some(name), Some(time), dist);
      node.rate_susceptibility_dates = Some(rate_dates);
      node
    }
  }
}
