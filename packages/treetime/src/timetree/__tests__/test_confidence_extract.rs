#![allow(
  clippy::as_conversions,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

#[cfg(test)]
mod tests {
  use crate::timetree::confidence::extract_confidence_intervals;
  use approx::assert_relative_eq;
  use helpers::add_named;
  use maplit::btreemap;
  use ndarray::Array1;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_graph::graph::Graph;

  #[test]
  fn test_extract_confidence_intervals_includes_unnamed_nodes() {
    let mut graph = Graph::new();
    let mut names = BTreeMap::new();
    let unnamed_key = add_named(&mut graph, &mut names, None);
    let named_key = add_named(&mut graph, &mut names, Some("named"));
    graph.build().unwrap();

    let state = helpers::state(
      &graph,
      &[(unnamed_key, Some(2020.0), None), (named_key, Some(2021.0), None)],
    );
    let intervals = extract_confidence_intervals(&graph, &state, &BTreeMap::new(), &names);
    assert_eq!(2, intervals.len());
    assert_eq!("", intervals[0].name);
    assert_eq!("named", intervals[1].name);
  }

  #[test]
  fn test_extract_confidence_intervals_skips_nodes_without_time() {
    let mut graph = Graph::new();
    let mut names = BTreeMap::new();
    let no_time_key = add_named(&mut graph, &mut names, Some("no_time"));
    let has_time_key = add_named(&mut graph, &mut names, Some("has_time"));
    graph.build().unwrap();

    let state = helpers::state(&graph, &[(no_time_key, None, None), (has_time_key, Some(2021.0), None)]);
    let intervals = extract_confidence_intervals(&graph, &state, &BTreeMap::new(), &names);
    assert_eq!(1, intervals.len());
    assert_eq!("has_time", intervals[0].name);
  }

  #[test]
  fn test_extract_confidence_intervals_uses_date_as_fallback() {
    let mut graph = Graph::new();
    let mut names = BTreeMap::new();
    let key = add_named(&mut graph, &mut names, Some("node_a"));
    graph.build().unwrap();

    let state = helpers::state(&graph, &[(key, Some(2020.5), None)]);
    let intervals = extract_confidence_intervals(&graph, &state, &BTreeMap::new(), &names);
    assert_eq!(1, intervals.len());
    assert_relative_eq!(intervals[0].date, 2020.5);
    assert_relative_eq!(intervals[0].lower, 2020.5);
    assert_relative_eq!(intervals[0].upper, 2020.5);
  }

  #[ignore = "marginal-posterior HPD disabled pending NegLog-aware HPD"]
  #[test]
  fn test_extract_confidence_intervals_with_distribution() {
    let mut graph = Graph::new();
    let mut names = BTreeMap::new();
    let dist = Arc::new(Distribution::range((2019.0, 2021.0), 0.0));
    let key = add_named(&mut graph, &mut names, Some("node_a"));
    graph.build().unwrap();

    let state = helpers::state(&graph, &[(key, Some(2020.0), Some(dist))]);
    let intervals = extract_confidence_intervals(&graph, &state, &BTreeMap::new(), &names);
    assert_eq!(1, intervals.len());
    assert_relative_eq!(intervals[0].date, 2020.0);
    assert_relative_eq!(intervals[0].lower, 2019.1, epsilon = 1e-10);
    assert_relative_eq!(intervals[0].upper, 2020.9, epsilon = 1e-10);
  }

  #[test]
  fn test_extract_confidence_intervals_sorted_by_key() {
    let mut graph = Graph::new();
    let mut names = BTreeMap::new();
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
    assert_eq!(3, intervals.len());
    assert_eq!("zebra", intervals[0].name);
    assert_eq!("alpha", intervals[1].name);
    assert_eq!("middle", intervals[2].name);
  }

  #[test]
  fn test_extract_confidence_intervals_rate_only() {
    let mut graph = Graph::new();
    let mut names = BTreeMap::new();
    let key = add_named(&mut graph, &mut names, Some("node_a"));
    graph.build().unwrap();
    let rate_map = btreemap! { key => [2009.0, 2010.0, 2011.0] };

    let state = helpers::state(&graph, &[(key, Some(2010.0), None)]);
    let intervals = extract_confidence_intervals(&graph, &state, &rate_map, &names);
    assert_eq!(1, intervals.len());
    assert_relative_eq!(intervals[0].date, 2010.0);
    assert_relative_eq!(intervals[0].lower, 2008.355146, epsilon = 1e-4);
    assert_relative_eq!(intervals[0].upper, 2011.644854, epsilon = 1e-4);
  }

  #[ignore = "marginal-posterior HPD disabled pending NegLog-aware HPD"]
  #[test]
  fn test_extract_confidence_intervals_combined_wider_than_either() {
    let mut graph = Graph::new();
    let mut names = BTreeMap::new();
    let dist = Arc::new(Distribution::range((2008.0, 2012.0), 0.0));
    let key = add_named(&mut graph, &mut names, Some("node_a"));
    graph.build().unwrap();
    let rate_map = btreemap! { key => [2009.0, 2010.0, 2011.0] };

    let state = helpers::state(&graph, &[(key, Some(2010.0), Some(dist))]);
    let intervals = extract_confidence_intervals(&graph, &state, &rate_map, &names);
    assert_eq!(1, intervals.len());
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
    let mut graph = Graph::new();
    let mut names = BTreeMap::new();
    let key = add_named(&mut graph, &mut names, Some("node_a"));
    graph.build().unwrap();
    let rate_map = btreemap! { key => [2019.9, 2020.0, 2020.1] };

    let state = helpers::state(&graph, &[(key, Some(2020.5), None)]);
    let intervals = extract_confidence_intervals(&graph, &state, &rate_map, &names);
    assert_eq!(1, intervals.len());
    assert!(intervals[0].lower <= intervals[0].date);
    assert!(intervals[0].date <= intervals[0].upper);
    assert_relative_eq!(intervals[0].upper, 2020.5);
    assert_relative_eq!(intervals[0].lower, 2019.8355, epsilon = 1e-3);
  }

  #[test]
  fn test_extract_confidence_intervals_clamps_when_date_below_rate_ci() {
    let mut graph = Graph::new();
    let mut names = BTreeMap::new();
    let key = add_named(&mut graph, &mut names, Some("node_a"));
    graph.build().unwrap();
    let rate_map = btreemap! { key => [2019.9, 2020.0, 2020.1] };

    let state = helpers::state(&graph, &[(key, Some(2019.5), None)]);
    let intervals = extract_confidence_intervals(&graph, &state, &rate_map, &names);
    assert_eq!(1, intervals.len());
    assert!(intervals[0].lower <= intervals[0].date);
    assert!(intervals[0].date <= intervals[0].upper);
    assert_relative_eq!(intervals[0].lower, 2019.5);
    assert_relative_eq!(intervals[0].upper, 2020.1645, epsilon = 1e-3);
  }

  #[ignore = "marginal-posterior HPD disabled pending NegLog-aware HPD"]
  #[test]
  fn test_extract_confidence_intervals_skewed_distribution_hpd() {
    let n_points = 500;
    let x_min = 0.0;
    let dx = 10.0 / (n_points as f64 - 1.0);
    let y = Array1::from_shape_fn(n_points, |i| x_min + i as f64 * dx);

    let dist_fn = treetime_distribution::DistributionFunction::from_start_dx_values(x_min, dx, y).unwrap();
    let dist = Distribution::Function(dist_fn);
    let peak_time = dist.likely_time().unwrap().unwrap();

    let mut graph = Graph::new();
    let mut names = BTreeMap::new();
    let key = add_named(&mut graph, &mut names, Some("skewed"));
    graph.build().unwrap();

    let state = helpers::state(&graph, &[(key, Some(peak_time), Some(Arc::new(dist)))]);
    let intervals = extract_confidence_intervals(&graph, &state, &BTreeMap::new(), &names);
    assert_eq!(1, intervals.len());

    let hpd_lower = 0.0;
    let hpd_upper = (0.1_f64).ln().abs();
    assert_relative_eq!(intervals[0].lower, hpd_lower, epsilon = 1e-3);
    assert_relative_eq!(intervals[0].upper, hpd_upper, epsilon = 1e-3);
  }

  mod helpers {
    use crate::timetree::timetree_state::TimetreeState;
    use std::collections::BTreeMap;
    use std::sync::Arc;
    use treetime_distribution::{Distribution, NegLog};
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;

    pub(super) type NodeTimeEntry = (GraphNodeKey, Option<f64>, Option<Arc<Distribution<NegLog>>>);

    pub(super) fn state(graph: &Graph, entries: &[NodeTimeEntry]) -> TimetreeState {
      let mut state = TimetreeState::new(graph);
      for (key, time, dist) in entries {
        let node = state.node_mut(*key);
        node.time = *time;
        node.time_distribution = dist.clone();
      }
      state
    }

    pub(super) fn add_named(
      graph: &mut Graph,
      names: &mut BTreeMap<GraphNodeKey, Option<String>>,
      name: Option<&str>,
    ) -> GraphNodeKey {
      let key = graph.add_node();
      names.insert(key, name.map(|n| n.to_owned()));
      key
    }
  }
}
