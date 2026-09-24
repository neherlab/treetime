#![allow(
  clippy::as_conversions,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

#[cfg(test)]
mod tests {
  use crate::pretty_assert_ulps_eq;
  use crate::test_utils::find_node_key_by_name;
  use crate::timetree::optimization::relaxed_clock::apply_relaxed_clock;
  use crate::timetree::timetree_state::TimetreeState;
  use eyre::Report;
  use rstest::rstest;
  use serde::Deserialize;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::io::json::json_read_str;
  use treetime_utils::pretty_assert_map_ulps_eq;

  use helpers::{GmInput, GmOutput, build_deep_tree, build_simple_tree, compute_variance, seed_state_scaled};

  #[rustfmt::skip]
  #[rstest]
  #[case::flu_rate_binary("flu_rate_binary")]
  #[case::flu_rate_deep(  "flu_rate_deep")]
  #[trace]
  fn test_gm_relaxed_clock_matches_v0(#[case] case_name: &str) -> Result<(), Report> {
    let inputs: Vec<GmInput> = json_read_str(include_str!("__fixtures__/gm_relaxed_clock_inputs.json"))?;
    let outputs: Vec<GmOutput> = json_read_str(include_str!("__fixtures__/gm_relaxed_clock_outputs.json"))?;
    let input = inputs
      .iter()
      .find(|input| input.name == case_name)
      .ok_or_else(|| eyre::eyre!("Golden-master input case {case_name} not found"))?;
    let expected = outputs
      .iter()
      .find(|output| output.name == case_name)
      .ok_or_else(|| eyre::eyre!("Golden-master output case {case_name} not found"))?;
    let nwk_parsed = nwk_read_str(&input.newick)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let mut state = TimetreeState::new(&graph);
    for (name, branch) in &input.branches {
      let node_key = find_node_key_by_name(&graph, &names, name).ok_or_else(|| eyre::eyre!("Node {name} not found"))?;
      let node = graph.get_node(node_key).ok_or_else(|| eyre::eyre!("Node {name} not found"))?;
      let (_, edge) = graph
        .parents_of(node)
        .next()
        .ok_or_else(|| eyre::eyre!("Parent edge for {name} not found"))?;
      state.edge_mut(edge.key()).time_length = Some(branch.clock_length / input.clock_rate);
    }

    apply_relaxed_clock(
      &graph,
      &branch_lengths,
      &[input.slack, input.coupling],
      input.one_mutation,
      input.clock_rate,
      &mut state,
    )?;

    let actual = input
      .branches
      .keys()
      .map(|name| {
        let node_key = find_node_key_by_name(&graph, &names, name).ok_or_else(|| eyre::eyre!("Node {name} not found"))?;
        let node = graph.get_node(node_key).ok_or_else(|| eyre::eyre!("Node {name} not found"))?;
        let (_, edge) = graph
          .parents_of(node)
          .next()
          .ok_or_else(|| eyre::eyre!("Parent edge for {name} not found"))?;
        let gamma = state.edge(edge.key()).gamma;
        Ok((name.to_owned(), gamma))
      })
      .collect::<Result<BTreeMap<_, _>, Report>>()?;

    pretty_assert_map_ulps_eq!(expected.gammas, actual, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_relaxed_clock_default_params_produce_reasonable_gamma() -> Result<(), Report> {
    let (graph, branch_lengths) = build_simple_tree()?;
    let one_mutation = 0.01;
    let params = [1.0, 1.0];

    let mut state = seed_state_scaled(&graph, &branch_lengths, 100.0);
    apply_relaxed_clock(&graph, &branch_lengths, &params, one_mutation, 1.0, &mut state)?;

    for edge in graph.get_edges() {
      let gamma = state.edge(edge.key()).gamma;

      assert!(gamma >= 0.1, "gamma={gamma} should be >= 0.1 (minimum bound)");
      assert!(gamma < 10.0, "gamma={gamma} should be reasonable (< 10.0)");
    }

    Ok(())
  }

  #[test]
  fn test_relaxed_clock_all_gamma_above_minimum() -> Result<(), Report> {
    let (graph, branch_lengths) = build_deep_tree()?;
    let one_mutation = 0.001;
    let params = [1.0, 1.0];

    let mut state = seed_state_scaled(&graph, &branch_lengths, 100.0);
    apply_relaxed_clock(&graph, &branch_lengths, &params, one_mutation, 1.0, &mut state)?;

    for edge in graph.get_edges() {
      let gamma = state.edge(edge.key()).gamma;
      assert!(gamma >= 0.1, "gamma={gamma} must be >= 0.1 (algorithm minimum bound)");
    }

    Ok(())
  }

  #[test]
  fn test_relaxed_clock_uniform_branches_produce_similar_gamma() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.1,C:0.1)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let one_mutation = 0.01;
    let params = [1.0, 1.0];
    let mut state = TimetreeState::new(&graph);
    for edge in graph.get_edges() {
      state.edge_mut(edge.key()).time_length = Some(10.0);
    }
    apply_relaxed_clock(&graph, &branch_lengths, &params, one_mutation, 1.0, &mut state)?;

    let gammas: Vec<f64> = graph.get_edges().map(|e| state.edge(e.key()).gamma).collect();

    let mean_gamma: f64 = gammas.iter().sum::<f64>() / gammas.len() as f64;
    for gamma in &gammas {
      pretty_assert_ulps_eq!(*gamma, mean_gamma, max_ulps = 1000);
    }

    Ok(())
  }

  #[test]
  fn test_relaxed_clock_empty_params_uses_defaults() -> Result<(), Report> {
    let (graph, branch_lengths) = build_simple_tree()?;
    let one_mutation = 0.01;

    let mut state = seed_state_scaled(&graph, &branch_lengths, 100.0);
    apply_relaxed_clock(&graph, &branch_lengths, &[], one_mutation, 1.0, &mut state)?;

    for edge in graph.get_edges() {
      let gamma = state.edge(edge.key()).gamma;
      assert!(gamma >= 0.1, "gamma should be >= minimum bound");
    }

    Ok(())
  }

  #[test]
  fn test_relaxed_clock_gamma_stored_in_edges() -> Result<(), Report> {
    let (graph, branch_lengths) = build_simple_tree()?;
    let one_mutation = 0.01;
    let params = [1.0, 1.0];

    let mut state = seed_state_scaled(&graph, &branch_lengths, 100.0);

    for edge in graph.get_edges() {
      let gamma = state.edge(edge.key()).gamma;
      pretty_assert_ulps_eq!(gamma, 1.0, max_ulps = 4);
    }

    apply_relaxed_clock(&graph, &branch_lengths, &params, one_mutation, 1.0, &mut state)?;

    let mut any_changed = false;
    for edge in graph.get_edges() {
      let gamma = state.edge(edge.key()).gamma;
      if (gamma - 1.0).abs() > 1e-6 {
        any_changed = true;
      }
    }

    assert!(any_changed, "At least one gamma should differ from default 1.0");

    Ok(())
  }

  #[test]
  fn test_relaxed_clock_high_slack_pulls_toward_one() -> Result<(), Report> {
    let (graph, branch_lengths) = build_deep_tree()?;
    let one_mutation = 0.01;

    let params_low = [1.0, 1.0];
    let mut state = seed_state_scaled(&graph, &branch_lengths, 100.0);
    apply_relaxed_clock(&graph, &branch_lengths, &params_low, one_mutation, 1.0, &mut state)?;

    let gammas_low: Vec<f64> = graph.get_edges().map(|e| state.edge(e.key()).gamma).collect();
    let deviation_low: f64 = gammas_low.iter().map(|g| (g - 1.0).abs()).sum();

    let params_high = [100.0, 1.0];
    apply_relaxed_clock(&graph, &branch_lengths, &params_high, one_mutation, 1.0, &mut state)?;

    let gammas_high: Vec<f64> = graph.get_edges().map(|e| state.edge(e.key()).gamma).collect();
    let deviation_high: f64 = gammas_high.iter().map(|g| (g - 1.0).abs()).sum();

    assert!(
      deviation_high <= deviation_low,
      "High slack should pull gammas toward 1.0: high_dev={deviation_high}, low_dev={deviation_low}"
    );

    Ok(())
  }

  #[test]
  fn test_relaxed_clock_high_coupling_reduces_variation() -> Result<(), Report> {
    let (graph, branch_lengths) = build_deep_tree()?;
    let one_mutation = 0.01;

    let params_low = [1.0, 0.1];
    let mut state = seed_state_scaled(&graph, &branch_lengths, 100.0);
    apply_relaxed_clock(&graph, &branch_lengths, &params_low, one_mutation, 1.0, &mut state)?;

    let gammas_low: Vec<f64> = graph.get_edges().map(|e| state.edge(e.key()).gamma).collect();
    let variance_low = compute_variance(&gammas_low);

    let params_high = [1.0, 10.0];
    apply_relaxed_clock(&graph, &branch_lengths, &params_high, one_mutation, 1.0, &mut state)?;

    let gammas_high: Vec<f64> = graph.get_edges().map(|e| state.edge(e.key()).gamma).collect();
    let variance_high = compute_variance(&gammas_high);

    assert!(
      variance_high <= variance_low,
      "High coupling should reduce gamma variance: high={variance_high}, low={variance_low}"
    );

    Ok(())
  }

  #[test]
  fn test_relaxed_clock_one_mutation_affects_gamma() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.01,B:0.02)AB:0.015,(C:0.005,D:0.01)CD:0.008)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let params = [1.0, 1.0];

    let one_mutation_single = 0.001;
    let mut state = seed_state_scaled(&graph, &branch_lengths, 0.8);
    apply_relaxed_clock(&graph, &branch_lengths, &params, one_mutation_single, 1.0, &mut state)?;

    let gammas_single: Vec<f64> = graph.get_edges().map(|e| state.edge(e.key()).gamma).collect();

    let one_mutation_multi = 0.0001;
    apply_relaxed_clock(&graph, &branch_lengths, &params, one_mutation_multi, 1.0, &mut state)?;

    let gammas_multi: Vec<f64> = graph.get_edges().map(|e| state.edge(e.key()).gamma).collect();

    let any_differ = gammas_single
      .iter()
      .zip(&gammas_multi)
      .any(|(s, m)| (s - m).abs() > 1e-10);

    assert!(
      any_differ,
      "Different one_mutation values should produce different gammas: single={gammas_single:?}, multi={gammas_multi:?}"
    );

    Ok(())
  }

  #[test]
  fn test_relaxed_clock_handles_tiny_one_mutation() -> Result<(), Report> {
    let (graph, branch_lengths) = build_simple_tree()?;
    let params = [1.0, 1.0];

    let tiny_one_mutation = 1e-15;
    let mut state = seed_state_scaled(&graph, &branch_lengths, 100.0);
    apply_relaxed_clock(&graph, &branch_lengths, &params, tiny_one_mutation, 1.0, &mut state)?;

    for edge in graph.get_edges() {
      let gamma = state.edge(edge.key()).gamma;
      assert!(gamma.is_finite(), "gamma must be finite, got {gamma}");
      assert!(gamma >= 0.1, "gamma must respect minimum bound, got {gamma}");
    }

    Ok(())
  }

  #[test]
  fn test_relaxed_clock_root_has_branch_penalty() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let one_mutation = 0.01;
    let params = [1.0, 1.0];
    let mut state = seed_state_scaled(&graph, &branch_lengths, 100.0);
    apply_relaxed_clock(&graph, &branch_lengths, &params, one_mutation, 1.0, &mut state)?;

    let root_edge_gamma = graph
      .get_edges()
      .collect::<Vec<_>>()
      .first()
      .map(|e| state.edge(e.key()).gamma);

    if let Some(gamma) = root_edge_gamma {
      assert!(gamma >= 0.1, "Root gamma should respect minimum bound: {gamma}");
      assert!(gamma < 10.0, "Root gamma should be reasonable: {gamma}");
    }

    Ok(())
  }

  #[rstest]
  #[trace]
  fn test_relaxed_clock_childless_root_gamma_equals_one(
    #[values(0.1, 1.0, 10.0, 100.0)] slack: f64,
    #[values(1e-6, 0.001, 0.01, 0.1, 1.0)] one_mutation: f64,
  ) -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let params = [slack, 1.0];
    let mut state = TimetreeState::new(&graph);
    apply_relaxed_clock(&graph, &branch_lengths, &params, one_mutation, 1.0, &mut state)?;

    for edge in graph.get_edges() {
      let gamma = state.edge(edge.key()).gamma;
      pretty_assert_ulps_eq!(gamma, 1.0, max_ulps = 4);
    }

    Ok(())
  }

  #[test]
  fn test_relaxed_clock_childless_root_gamma_stored() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.01)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let one_mutation = 0.01;
    let params = [1.0, 1.0];
    let mut state = seed_state_scaled(&graph, &branch_lengths, 1.0);
    apply_relaxed_clock(&graph, &branch_lengths, &params, one_mutation, 1.0, &mut state)?;

    let gamma = graph
      .get_edges()
      .collect::<Vec<_>>()
      .first()
      .map_or(1.0, |e| state.edge(e.key()).gamma);

    pretty_assert_ulps_eq!(gamma, 1.0, max_ulps = 100);

    Ok(())
  }

  mod helpers {
    use super::*;

    #[derive(Deserialize)]
    pub(super) struct GmInput {
      pub name: String,
      pub newick: String,
      pub clock_rate: f64,
      pub one_mutation: f64,
      pub slack: f64,
      pub coupling: f64,
      pub branches: BTreeMap<String, GmBranchInput>,
    }

    #[derive(Deserialize)]
    pub(super) struct GmBranchInput {
      pub clock_length: f64,
    }

    #[derive(Deserialize)]
    pub(super) struct GmOutput {
      pub name: String,
      pub gammas: BTreeMap<String, f64>,
    }

    pub(super) fn build_simple_tree() -> Result<(Graph, BTreeMap<GraphEdgeKey, Option<f64>>), Report> {
      let nwk_parsed = nwk_read_str("(A:0.1,B:0.2)root:0.0;")?;
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      Ok((graph, branch_lengths))
    }

    pub(super) fn build_deep_tree() -> Result<(Graph, BTreeMap<GraphEdgeKey, Option<f64>>), Report> {
      let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.15,(C:0.05,D:0.1)CD:0.08)root:0.0;")?;
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      Ok((graph, branch_lengths))
    }

    pub(super) fn seed_state_scaled(
      graph: &Graph,
      branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
      factor: f64,
    ) -> TimetreeState {
      let mut state = TimetreeState::new(graph);
      for (&edge_key, &branch_length) in branch_lengths {
        state.edge_mut(edge_key).time_length = Some(branch_length.unwrap_or(0.0) * factor);
      }
      state
    }

    pub(super) fn compute_variance(values: &[f64]) -> f64 {
      let n = values.len() as f64;
      let mean = values.iter().sum::<f64>() / n;
      values.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n
    }
  }
}
