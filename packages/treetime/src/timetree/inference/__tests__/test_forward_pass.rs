#![allow(
  clippy::wildcard_enum_match_arm,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::DateConstraints;
  use crate::pretty_assert_ulps_eq;
  use crate::test_utils::find_node_key_by_name;
  use crate::timetree::inference::forward_pass::{propagate_distributions_forward, set_likely_time};
  use crate::timetree::timetree_state::{DateNodeState, TimetreeState};
  use eyre::Report;
  use ndarray::Array1;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_distribution::{Distribution, NegLog};
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read_str;

  type TestGraph = Graph;

  #[test]
  fn test_forward_pass_set_likely_time_empty_distribution_returns_none() {
    let mut node = node_with_distribution(Some(Distribution::empty()));
    assert_eq!(None, set_likely_time(&mut node, None));
    assert_eq!(None, node.time);
  }

  #[test]
  fn test_forward_pass_set_likely_time_missing_distribution_returns_none() {
    let mut node = node_with_distribution(None);
    assert_eq!(None, set_likely_time(&mut node, None));
    assert_eq!(None, node.time);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::no_parent(      None,      5.0)]
  #[case::parent_earlier( Some(3.0), 5.0)]
  #[case::parent_clamps(  Some(8.0), 8.0)]
  #[trace]
  fn test_forward_pass_set_likely_time_uses_distribution_peak(
    #[case] parent_time: Option<f64>,
    #[case] expected: f64,
  ) {
    let mut node = node_with_distribution(Some(Distribution::point(5.0, 1.0)));
    let assigned = set_likely_time(&mut node, parent_time).expect("a time should be assigned");
    pretty_assert_ulps_eq!(assigned, expected, max_ulps = 4);
    let committed = node.time.expect("node time should be committed");
    pretty_assert_ulps_eq!(committed, expected, max_ulps = 4);
  }

  #[test]
  fn test_forward_pass_leaves_internal_node_with_empty_distribution_undated() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:2.5)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;

    let root_key = find_node_key_by_name(&graph, &names, "root").expect("root not found");
    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");

    let mut constraints = DateConstraints::default();
    let mut state = TimetreeState::new(&graph);
    state.node_mut(root_key).time_distribution = Some(Arc::new(Distribution::empty()));
    set_date(&mut constraints, &mut state, leaf_key, Distribution::point(2013.0, 0.0));

    let state = run_forward_pass(&graph, &constraints, &names, state)?;

    let root_time = node_time(&state, root_key);

    assert_eq!(None, root_time);
    let leaf_time = node_time(&state, leaf_key).expect("leaf A should keep its observed date");
    pretty_assert_ulps_eq!(leaf_time, 2013.0, max_ulps = 4);

    Ok(())
  }

  #[test]
  fn test_forward_pass_refines_uncertain_leaf_date_from_parent() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let root_key = find_node_key_by_name(&graph, &names, "root").expect("root not found");
    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");

    let mut constraints = DateConstraints::default();
    let mut state = TimetreeState::new(&graph);
    set_time_distribution(&mut state, root_key, Distribution::point(2009.0, 0.0));
    set_date(
      &mut constraints,
      &mut state,
      leaf_key,
      Distribution::range((2009.5, 2013.5), 0.0),
    );
    set_branch_length_distribution(&graph, &mut state, leaf_key, 1.0);

    let state = run_forward_pass(&graph, &constraints, &names, state)?;

    let leaf_time = node_time(&state, leaf_key).expect("leaf A should be dated");
    pretty_assert_ulps_eq!(leaf_time, 2010.0, max_ulps = 4);

    let leaf_dist = leaf_time_distribution(&state, leaf_key).expect("leaf A should have a distribution");
    assert_eq!(Distribution::point(2010.0, 0.0), leaf_dist);

    Ok(())
  }

  #[test]
  fn test_forward_pass_keeps_exact_leaf_date_earlier_than_its_parent() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let root_key = find_node_key_by_name(&graph, &names, "root").expect("root not found");
    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");

    let mut constraints = DateConstraints::default();
    let mut state = TimetreeState::new(&graph);
    set_time_distribution(&mut state, root_key, Distribution::point(2009.0, 0.0));
    set_date(&mut constraints, &mut state, leaf_key, Distribution::point(2008.0, 0.0));
    set_branch_length_distribution(&graph, &mut state, leaf_key, 1.0);

    let state = run_forward_pass(&graph, &constraints, &names, state)?;

    let leaf_time = node_time(&state, leaf_key).expect("leaf A should keep its observed date");
    pretty_assert_ulps_eq!(leaf_time, 2008.0, max_ulps = 4);

    let leaf_dist = leaf_time_distribution(&state, leaf_key).expect("leaf A should have a distribution");
    assert_eq!(Distribution::point(2008.0, 0.0), leaf_dist);

    Ok(())
  }

  #[test]
  fn test_forward_pass_clamps_uncertain_leaf_date_to_parent_time() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let root_key = find_node_key_by_name(&graph, &names, "root").expect("root not found");
    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");

    let mut constraints = DateConstraints::default();
    let mut state = TimetreeState::new(&graph);
    set_time_distribution(&mut state, root_key, Distribution::point(2009.0, 0.0));
    set_date(
      &mut constraints,
      &mut state,
      leaf_key,
      Distribution::range((2005.0, 2007.0), 0.0),
    );

    let state = run_forward_pass(&graph, &constraints, &names, state)?;

    let leaf_time = node_time(&state, leaf_key).expect("leaf A should be dated");
    pretty_assert_ulps_eq!(leaf_time, 2009.0, max_ulps = 4);

    Ok(())
  }

  #[test]
  fn test_forward_pass_keeps_uncertain_leaf_date_the_tree_contradicts() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let root_key = find_node_key_by_name(&graph, &names, "root").expect("root not found");
    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");

    let mut constraints = DateConstraints::default();
    let mut state = TimetreeState::new(&graph);
    set_time_distribution(&mut state, root_key, Distribution::point(2009.0, 0.0));
    let given = Distribution::range((2005.0, 2007.0), 0.0);
    set_date(&mut constraints, &mut state, leaf_key, given.clone());
    set_branch_length_distribution(&graph, &mut state, leaf_key, 1.0);

    let state = run_forward_pass(&graph, &constraints, &names, state)?;

    let leaf_dist = leaf_time_distribution(&state, leaf_key).expect("leaf A should keep its given date");
    assert_eq!(given, leaf_dist);

    let leaf_time = node_time(&state, leaf_key).expect("leaf A should be dated");
    pretty_assert_ulps_eq!(leaf_time, 2009.0, max_ulps = 4);

    Ok(())
  }

  #[test]
  fn test_forward_pass_refinement_is_unchanged_by_the_reachable_window() -> Result<(), Report> {
    let wide_parent = {
      let t = Array1::linspace(1950.0, 2050.0, 4001);
      let y = t.mapv(|t: f64| 0.5 * ((t - 2009.0) / 3.0).powi(2));
      Distribution::function(t, y)?
    };

    let refine = |parent: Distribution<NegLog>| -> Result<f64, Report> {
      let nwk_parsed = nwk_read_str("(A:1.0)root;")?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let root_key = find_node_key_by_name(&graph, &names, "root").expect("root not found");
      let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
      let mut constraints = DateConstraints::default();
      let mut state = TimetreeState::new(&graph);
      set_time_distribution(&mut state, root_key, parent);
      set_date(
        &mut constraints,
        &mut state,
        leaf_key,
        Distribution::range((2010.5, 2010.6), 0.0),
      );
      set_branch_length_distribution(&graph, &mut state, leaf_key, 1.0);
      let state = run_forward_pass(&graph, &constraints, &names, state)?;
      node_time(&state, leaf_key).ok_or_else(|| eyre::eyre!("leaf A should be dated"))
    };

    let windowed = match &wide_parent {
      Distribution::Function(f) => Distribution::Function(f.resample_range_dx((2009.4, 2009.7), f.dx())?),
      other => other.clone(),
    };

    pretty_assert_ulps_eq!(refine(wide_parent)?, refine(windowed)?, max_ulps = 4);

    Ok(())
  }

  #[test]
  fn test_forward_pass_refines_a_date_range_narrower_than_the_parent_grid() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let root_key = find_node_key_by_name(&graph, &names, "root").expect("root not found");
    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");

    let coarse_parent = {
      let t = Array1::linspace(1900.0, 2020.0, 121);
      let y = t.mapv(|t: f64| 0.5 * ((t - 2009.0) / 5.0).powi(2));
      Distribution::function(t, y)?
    };
    let mut constraints = DateConstraints::default();
    let mut state = TimetreeState::new(&graph);
    set_time_distribution(&mut state, root_key, coarse_parent);
    set_date(
      &mut constraints,
      &mut state,
      leaf_key,
      Distribution::range((2010.500, 2010.503), 0.0),
    );
    set_branch_length_distribution(&graph, &mut state, leaf_key, 1.0);

    let state = run_forward_pass(&graph, &constraints, &names, state)?;

    let leaf_time = node_time(&state, leaf_key).expect("leaf A should be dated");
    assert!(
      (2010.500..=2010.503).contains(&leaf_time),
      "leaf A should be dated within the day it was given, got {leaf_time}"
    );

    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn run_forward_pass(
      graph: &TestGraph,
      constraints: &DateConstraints,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      mut state: TimetreeState,
    ) -> Result<TimetreeState, Report> {
      propagate_distributions_forward(graph, constraints, names, &mut state)?;
      Ok(state)
    }

    pub(super) fn node_with_distribution(distribution: Option<Distribution<NegLog>>) -> DateNodeState {
      DateNodeState {
        time_distribution: distribution.map(Arc::new),
        ..DateNodeState::default()
      }
    }

    pub(super) fn set_date(
      constraints: &mut DateConstraints,
      state: &mut TimetreeState,
      key: GraphNodeKey,
      dist: Distribution<NegLog>,
    ) {
      let dist = Arc::new(dist);
      constraints.date_constraints.insert(key, Some(Arc::clone(&dist)));
      state.node_mut(key).time_distribution = Some(dist);
    }

    pub(super) fn set_time_distribution(state: &mut TimetreeState, key: GraphNodeKey, dist: Distribution<NegLog>) {
      state.node_mut(key).time_distribution = Some(Arc::new(dist));
    }

    pub(super) fn set_branch_length_distribution(
      graph: &TestGraph,
      state: &mut TimetreeState,
      target_key: GraphNodeKey,
      branch_length: f64,
    ) {
      for edge in graph.get_edges() {
        if edge.target() == target_key {
          state.edge_mut(edge.key()).branch_length_distribution =
            Some(Arc::new(Distribution::point(branch_length, 0.0)));
        }
      }
    }

    pub(super) fn node_time(state: &TimetreeState, key: GraphNodeKey) -> Option<f64> {
      state.nodes.get(&key).and_then(|node| node.time)
    }

    pub(super) fn leaf_time_distribution(state: &TimetreeState, key: GraphNodeKey) -> Option<Distribution<NegLog>> {
      state
        .nodes
        .get(&key)
        .and_then(|node| node.time_distribution.as_ref())
        .map(|dist| dist.as_ref().clone())
    }
  }

  use helpers::*;
}
