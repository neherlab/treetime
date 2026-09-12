#[cfg(test)]
mod tests {
  use crate::partition::timetree::partition::GraphTimetree;
  use crate::timetree::confidence::{extract_confidence_intervals, write_confidence_intervals};
  use approx::assert_relative_eq;
  use helpers::add_named;
  use maplit::btreemap;
  use ndarray::Array1;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_distribution::Distribution;

  #[test]
  fn test_extract_confidence_intervals_includes_unnamed_nodes() {
    let mut graph = GraphTimetree::new();
    let mut names = BTreeMap::new();
    let unnamed_key = add_named(&mut graph, &mut names, None);
    let named_key = add_named(&mut graph, &mut names, Some("named"));
    graph.build().unwrap();

    let state = helpers::state(
      &graph,
      &[(unnamed_key, Some(2020.0), None), (named_key, Some(2021.0), None)],
    );
    let intervals = extract_confidence_intervals(&graph, &state, &BTreeMap::new(), &names);
    assert_eq!(intervals.len(), 2);
    // Unnamed node has empty name but valid key
    assert_eq!(intervals[0].name, "");
    assert_eq!(intervals[1].name, "named");
  }

  #[test]
  fn test_extract_confidence_intervals_skips_nodes_without_time() {
    let mut graph = GraphTimetree::new();
    let mut names = BTreeMap::new();
    let no_time_key = add_named(&mut graph, &mut names, Some("no_time"));
    let has_time_key = add_named(&mut graph, &mut names, Some("has_time"));
    graph.build().unwrap();

    let state = helpers::state(&graph, &[(no_time_key, None, None), (has_time_key, Some(2021.0), None)]);
    let intervals = extract_confidence_intervals(&graph, &state, &BTreeMap::new(), &names);
    assert_eq!(intervals.len(), 1);
    assert_eq!(intervals[0].name, "has_time");
  }

  #[test]
  fn test_extract_confidence_intervals_uses_date_as_fallback() {
    let mut graph = GraphTimetree::new();
    let mut names = BTreeMap::new();
    let key = add_named(&mut graph, &mut names, Some("node_a"));
    graph.build().unwrap();

    let state = helpers::state(&graph, &[(key, Some(2020.5), None)]);
    let intervals = extract_confidence_intervals(&graph, &state, &BTreeMap::new(), &names);
    assert_eq!(intervals.len(), 1);
    assert_relative_eq!(intervals[0].date, 2020.5);
    assert_relative_eq!(intervals[0].lower, 2020.5);
    assert_relative_eq!(intervals[0].upper, 2020.5);
  }

  // Ignored: the marginal-posterior HPD contribution is disabled in `extract_confidence_intervals`
  // until a NegLog-aware HPD region lands. Without it a node whose only source is a time
  // distribution falls back to the point estimate, so this interval collapses to `[date, date]`.
  #[ignore = "marginal-posterior HPD disabled pending NegLog-aware HPD"]
  #[test]
  fn test_extract_confidence_intervals_with_distribution() {
    let mut graph = GraphTimetree::new();
    let mut names = BTreeMap::new();
    let dist = Arc::new(Distribution::range((2019.0, 2021.0), 0.0));
    let key = add_named(&mut graph, &mut names, Some("node_a"));
    graph.build().unwrap();

    let state = helpers::state(&graph, &[(key, Some(2020.0), Some(dist))]);
    let intervals = extract_confidence_intervals(&graph, &state, &BTreeMap::new(), &names);
    assert_eq!(intervals.len(), 1);
    assert_relative_eq!(intervals[0].date, 2020.0);
    // 90% CI from uniform [2019, 2021]: 0.05 * 2 + 2019 = 2019.1, 0.95 * 2 + 2019 = 2020.9
    assert_relative_eq!(intervals[0].lower, 2019.1, epsilon = 1e-10);
    assert_relative_eq!(intervals[0].upper, 2020.9, epsilon = 1e-10);
  }

  #[test]
  fn test_extract_confidence_intervals_sorted_by_key() {
    let mut graph = GraphTimetree::new();
    let mut names = BTreeMap::new();
    // Insertion order: zebra (key 0), alpha (key 1), middle (key 2)
    let zebra_key = add_named(&mut graph, &mut names, Some("zebra"));
    let alpha_key = add_named(&mut graph, &mut names, Some("alpha"));
    let middle_key = add_named(&mut graph, &mut names, Some("middle"));
    graph.build().unwrap();

    let state = helpers::state(
      &graph,
      &[
        (zebra_key, Some(2020.0), None),
        (alpha_key, Some(2021.0), None),
        (middle_key, Some(2022.0), None),
      ],
    );
    let intervals = extract_confidence_intervals(&graph, &state, &BTreeMap::new(), &names);
    assert_eq!(intervals.len(), 3);
    // Sorted by GraphNodeKey (insertion order), not alphabetical
    assert_eq!(intervals[0].name, "zebra");
    assert_eq!(intervals[1].name, "alpha");
    assert_eq!(intervals[2].name, "middle");
  }

  #[test]
  fn test_extract_confidence_intervals_rate_only() {
    // Rate susceptibility data but no marginal distribution
    let mut graph = GraphTimetree::new();
    let mut names = BTreeMap::new();
    let key = add_named(&mut graph, &mut names, Some("node_a"));
    graph.build().unwrap();
    let rate_map = btreemap! { key => [2009.0, 2010.0, 2011.0] };

    let state = helpers::state(&graph, &[(key, Some(2010.0), None)]);
    let intervals = extract_confidence_intervals(&graph, &state, &rate_map, &names);
    assert_eq!(intervals.len(), 1);
    assert_relative_eq!(intervals[0].date, 2010.0);
    // z-score at 0.05 = -1.644854, at 0.95 = +1.644854
    // lower = 2010 + (-1.644854) * |2009 - 2010| = 2008.355146
    // upper = 2010 + 1.644854 * |2011 - 2010| = 2011.644854
    assert_relative_eq!(intervals[0].lower, 2008.355146, epsilon = 1e-4);
    assert_relative_eq!(intervals[0].upper, 2011.644854, epsilon = 1e-4);
  }

  // Ignored: needs the marginal-posterior HPD contribution, which is disabled in
  // `extract_confidence_intervals` until a NegLog-aware HPD region lands. Without it only the rate
  // contribution remains, so the interval cannot be wider than the mutation source alone.
  #[ignore = "marginal-posterior HPD disabled pending NegLog-aware HPD"]
  #[test]
  fn test_extract_confidence_intervals_combined_wider_than_either() {
    // Both marginal distribution and rate susceptibility data present.
    // The quadrature combination must be wider than either source alone.
    let mut graph = GraphTimetree::new();
    let mut names = BTreeMap::new();
    let dist = Arc::new(Distribution::range((2008.0, 2012.0), 0.0));
    let key = add_named(&mut graph, &mut names, Some("node_a"));
    graph.build().unwrap();
    let rate_map = btreemap! { key => [2009.0, 2010.0, 2011.0] };

    let state = helpers::state(&graph, &[(key, Some(2010.0), Some(dist))]);
    let intervals = extract_confidence_intervals(&graph, &state, &rate_map, &names);
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
    let mut names = BTreeMap::new();
    // date = 2020.5 (final pass), rate susceptibility centered on 2020.0
    // with small variation [2019.9, 2020.0, 2020.1].
    // Rate CI at 90%: 2020.0 +/- 1.645 * 0.1 = [2019.836, 2020.164]
    // date = 2020.5 > 2020.164, so upper must be clamped to date.
    let key = add_named(&mut graph, &mut names, Some("node_a"));
    graph.build().unwrap();
    let rate_map = btreemap! { key => [2019.9, 2020.0, 2020.1] };

    let state = helpers::state(&graph, &[(key, Some(2020.5), None)]);
    let intervals = extract_confidence_intervals(&graph, &state, &rate_map, &names);
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
    let mut names = BTreeMap::new();
    // date = 2019.5 (final pass), rate susceptibility centered on 2020.0
    // Rate CI at 90%: [2019.836, 2020.164]
    // date = 2019.5 < 2019.836, so lower must be clamped to date.
    let key = add_named(&mut graph, &mut names, Some("node_a"));
    graph.build().unwrap();
    let rate_map = btreemap! { key => [2019.9, 2020.0, 2020.1] };

    let state = helpers::state(&graph, &[(key, Some(2019.5), None)]);
    let intervals = extract_confidence_intervals(&graph, &state, &rate_map, &names);
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

  // Ignored: exercises the marginal-posterior HPD region directly, which is disabled in
  // `extract_confidence_intervals` until a NegLog-aware HPD lands. The distribution now stores
  // neg-log ordinates so the test is ready to re-enable once that HPD path returns.
  #[ignore = "marginal-posterior HPD disabled pending NegLog-aware HPD"]
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
    // P(t) = exp(-t) stored on the neg-log axis: the ordinate is `-ln P(t) = t`, so the peak
    // (minimum ordinate) sits at t = 0 and the long right tail rises linearly.
    let y = Array1::from_shape_fn(n_points, |i| x_min + i as f64 * dx);

    let dist_fn = treetime_distribution::DistributionFunction::from_start_dx_values(x_min, dx, y).unwrap();
    let dist = Distribution::Function(dist_fn);
    let peak_time = dist.likely_time().unwrap();

    let mut graph = GraphTimetree::new();
    let mut names = BTreeMap::new();
    let key = add_named(&mut graph, &mut names, Some("skewed"));
    graph.build().unwrap();

    let state = helpers::state(&graph, &[(key, Some(peak_time), Some(Arc::new(dist)))]);
    let intervals = extract_confidence_intervals(&graph, &state, &BTreeMap::new(), &names);
    assert_eq!(intervals.len(), 1);

    // v0 HPD bounds: [0, 2.3026] (narrowest 90% interval around peak)
    let hpd_lower = 0.0;
    let hpd_upper = (0.1_f64).ln().abs(); // -ln(0.1) = 2.3026
    // Measured errors: lower=0.0, upper=3.9e-4.
    assert_relative_eq!(intervals[0].lower, hpd_lower, epsilon = 1e-3);
    assert_relative_eq!(intervals[0].upper, hpd_upper, epsilon = 1e-3);
  }

  #[test]
  fn test_write_confidence_intervals_omits_internal_key_column() {
    // The confidence TSV mirrors the augur node-data contract: columns are
    // name, date, lower, upper. The graph node key is internal (serde-skipped)
    // and must never leak as a serialized column.
    let mut graph = GraphTimetree::new();
    let mut names = BTreeMap::new();
    let key = add_named(&mut graph, &mut names, Some("named"));
    graph.build().unwrap();
    let state = helpers::state(&graph, &[(key, Some(2020.0), None)]);
    let intervals = extract_confidence_intervals(&graph, &state, &BTreeMap::new(), &names);

    let mut buf = Vec::new();
    write_confidence_intervals(&intervals, &mut buf).unwrap();
    let output = String::from_utf8(buf).unwrap();

    let header = output.lines().next().unwrap();
    assert_eq!(header, "name\tdate\tlower\tupper");
  }

  mod helpers {
    use crate::partition::timetree::partition::GraphTimetree;
    use crate::timetree::timetree_state::TimetreeState;
    use std::collections::BTreeMap;
    use std::sync::Arc;
    use treetime_distribution::{Distribution, NegLog};
    use treetime_graph::node::GraphNodeKey;

    /// Per-node committed time and time distribution the confidence extraction reads, keyed by node.
    pub type NodeTimeEntry = (GraphNodeKey, Option<f64>, Option<Arc<Distribution<NegLog>>>);

    /// Date state built from the test nodes' committed times and distributions as values, so
    /// `extract_confidence_intervals` reads them from the state rather than the graph payload.
    pub fn state(graph: &GraphTimetree, entries: &[NodeTimeEntry]) -> TimetreeState {
      let mut state = TimetreeState::new(graph);
      for (key, time, dist) in entries {
        let node = state.node_mut(*key);
        node.time = *time;
        node.time_distribution = dist.clone();
      }
      state
    }

    pub fn add_named(
      graph: &mut GraphTimetree,
      names: &mut BTreeMap<GraphNodeKey, Option<String>>,
      name: Option<&str>,
    ) -> GraphNodeKey {
      let key = graph.add_node();
      names.insert(key, name.map(|n| n.to_owned()));
      key
    }
  }
}
