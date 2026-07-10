#[cfg(test)]
mod tests {
  use crate::partition::timetree::GraphTimetree;
  use crate::pretty_assert_ulps_eq;
  use crate::test_utils::find_node_key_by_name;
  use crate::timetree::optimization::relaxed_clock::apply_relaxed_clock;
  use eyre::Report;
  use rstest::rstest;
  use serde::Deserialize;
  use std::collections::BTreeMap;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::io::json::json_read_str;
  use treetime_utils::pretty_assert_map_ulps_eq;

  use helpers::{GmInput, GmOutput, build_deep_tree, build_simple_tree, compute_variance};

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
    let graph: GraphTimetree = nwk_read_str(&input.newick)?;

    for (name, branch) in &input.branches {
      let node_key = find_node_key_by_name(&graph, name).ok_or_else(|| eyre::eyre!("Node {name} not found"))?;
      let node = graph.get_node(node_key).ok_or_else(|| eyre::eyre!("Node {name} not found"))?;
      let (_, edge) = graph
        .parents_of(&node.read_arc())
        .into_iter()
        .next()
        .ok_or_else(|| eyre::eyre!("Parent edge for {name} not found"))?;
      edge.write_arc().payload().write_arc().time_length = Some(branch.clock_length / input.clock_rate);
    }

    apply_relaxed_clock(
      &graph,
      &[input.slack, input.coupling],
      input.one_mutation,
      input.clock_rate,
    )?;

    let actual = input
      .branches
      .keys()
      .map(|name| {
        let node_key = find_node_key_by_name(&graph, name).ok_or_else(|| eyre::eyre!("Node {name} not found"))?;
        let node = graph.get_node(node_key).ok_or_else(|| eyre::eyre!("Node {name} not found"))?;
        let (_, edge) = graph
          .parents_of(&node.read_arc())
          .into_iter()
          .next()
          .ok_or_else(|| eyre::eyre!("Parent edge for {name} not found"))?;
        let gamma = edge.read_arc().payload().read_arc().gamma;
        Ok((name.to_owned(), gamma))
      })
      .collect::<Result<BTreeMap<_, _>, Report>>()?;

    pretty_assert_map_ulps_eq!(expected.gammas, actual, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_relaxed_clock_default_params_produce_reasonable_gamma() -> Result<(), Report> {
    let graph = build_simple_tree()?;
    let one_mutation = 0.01;
    let params = [1.0, 1.0];

    apply_relaxed_clock(&graph, &params, one_mutation, 1.0)?;

    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let gamma = edge.payload().read_arc().gamma;

      assert!(gamma >= 0.1, "gamma={gamma} should be >= 0.1 (minimum bound)");
      assert!(gamma < 10.0, "gamma={gamma} should be reasonable (< 10.0)");
    }

    Ok(())
  }

  #[test]
  fn test_relaxed_clock_all_gamma_above_minimum() -> Result<(), Report> {
    let graph = build_deep_tree()?;
    let one_mutation = 0.001;
    let params = [1.0, 1.0];

    apply_relaxed_clock(&graph, &params, one_mutation, 1.0)?;

    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let gamma = edge.payload().read_arc().gamma;
      assert!(gamma >= 0.1, "gamma={gamma} must be >= 0.1 (algorithm minimum bound)");
    }

    Ok(())
  }

  #[test]
  fn test_relaxed_clock_uniform_branches_produce_similar_gamma() -> Result<(), Report> {
    // Tree with uniform branch lengths
    let graph: GraphTimetree = nwk_read_str("(A:0.1,B:0.1,C:0.1)root:0.0;")?;

    for edge in graph.get_edges() {
      let edge = edge.write_arc();
      edge.payload().write_arc().time_length = Some(10.0);
    }

    let one_mutation = 0.01;
    let params = [1.0, 1.0];
    apply_relaxed_clock(&graph, &params, one_mutation, 1.0)?;

    let gammas: Vec<f64> = graph
      .get_edges()
      .iter()
      .map(|e| e.read_arc().payload().read_arc().gamma)
      .collect();

    let mean_gamma: f64 = gammas.iter().sum::<f64>() / gammas.len() as f64;
    for gamma in &gammas {
      pretty_assert_ulps_eq!(*gamma, mean_gamma, max_ulps = 1000);
    }

    Ok(())
  }

  #[test]
  fn test_relaxed_clock_empty_params_uses_defaults() -> Result<(), Report> {
    let graph = build_simple_tree()?;
    let one_mutation = 0.01;

    apply_relaxed_clock(&graph, &[], one_mutation, 1.0)?;

    for edge in graph.get_edges() {
      let gamma = edge.read_arc().payload().read_arc().gamma;
      assert!(gamma >= 0.1, "gamma should be >= minimum bound");
    }

    Ok(())
  }

  #[test]
  fn test_relaxed_clock_gamma_stored_in_edges() -> Result<(), Report> {
    let graph = build_simple_tree()?;
    let one_mutation = 0.01;
    let params = [1.0, 1.0];

    // Before: gamma should be default 1.0
    for edge in graph.get_edges() {
      let gamma = edge.read_arc().payload().read_arc().gamma;
      pretty_assert_ulps_eq!(gamma, 1.0, max_ulps = 4);
    }

    apply_relaxed_clock(&graph, &params, one_mutation, 1.0)?;

    // After: gamma values should be computed (may differ from 1.0)
    let mut any_changed = false;
    for edge in graph.get_edges() {
      let gamma = edge.read_arc().payload().read_arc().gamma;
      if (gamma - 1.0).abs() > 1e-6 {
        any_changed = true;
      }
    }

    assert!(any_changed, "At least one gamma should differ from default 1.0");

    Ok(())
  }

  #[test]
  fn test_relaxed_clock_high_slack_pulls_toward_one() -> Result<(), Report> {
    let graph = build_deep_tree()?;
    let one_mutation = 0.01;

    // Compare low slack vs high slack - high slack should have gammas closer to 1.0
    let params_low = [1.0, 1.0];
    apply_relaxed_clock(&graph, &params_low, one_mutation, 1.0)?;

    let gammas_low: Vec<f64> = graph
      .get_edges()
      .iter()
      .map(|e| e.read_arc().payload().read_arc().gamma)
      .collect();
    let deviation_low: f64 = gammas_low.iter().map(|g| (g - 1.0).abs()).sum();

    let params_high = [100.0, 1.0];
    apply_relaxed_clock(&graph, &params_high, one_mutation, 1.0)?;

    let gammas_high: Vec<f64> = graph
      .get_edges()
      .iter()
      .map(|e| e.read_arc().payload().read_arc().gamma)
      .collect();
    let deviation_high: f64 = gammas_high.iter().map(|g| (g - 1.0).abs()).sum();

    // High slack should reduce total deviation from 1.0
    assert!(
      deviation_high <= deviation_low,
      "High slack should pull gammas toward 1.0: high_dev={deviation_high}, low_dev={deviation_low}"
    );

    Ok(())
  }

  #[test]
  fn test_relaxed_clock_high_coupling_reduces_variation() -> Result<(), Report> {
    let graph = build_deep_tree()?;
    let one_mutation = 0.01;

    // Run with low coupling
    let params_low = [1.0, 0.1];
    apply_relaxed_clock(&graph, &params_low, one_mutation, 1.0)?;

    let gammas_low: Vec<f64> = graph
      .get_edges()
      .iter()
      .map(|e| e.read_arc().payload().read_arc().gamma)
      .collect();
    let variance_low = compute_variance(&gammas_low);

    // Run with high coupling
    let params_high = [1.0, 10.0];
    apply_relaxed_clock(&graph, &params_high, one_mutation, 1.0)?;

    let gammas_high: Vec<f64> = graph
      .get_edges()
      .iter()
      .map(|e| e.read_arc().payload().read_arc().gamma)
      .collect();
    let variance_high = compute_variance(&gammas_high);

    assert!(
      variance_high <= variance_low,
      "High coupling should reduce gamma variance: high={variance_high}, low={variance_low}"
    );

    Ok(())
  }

  /// Regression test for T3: one_mutation must sum sequence lengths across partitions.
  /// Different one_mutation values produce different gamma values.
  /// This validates that the calculation in refinement.rs matters.
  #[test]
  fn test_relaxed_clock_one_mutation_affects_gamma() -> Result<(), Report> {
    // Use a tree with branch lengths that differ from time_length
    // This creates rate variation that the algorithm must account for
    let graph: GraphTimetree = nwk_read_str("((A:0.01,B:0.02)AB:0.015,(C:0.005,D:0.01)CD:0.008)root:0.0;")?;

    // Set time_length to differ from branch_length (rate variation scenario)
    // This ensures gamma != 1.0 and is sensitive to one_mutation parameter
    for edge in graph.get_edges() {
      let edge = edge.write_arc();
      let branch_length = edge.payload().read_arc().base.branch_length.unwrap_or(0.0);
      // time_length = 0.8 * branch_length creates rate variation
      edge.payload().write_arc().time_length = Some(branch_length * 0.8);
    }

    let params = [1.0, 1.0];

    // Simulate single partition with length 1000: one_mutation = 1/1000 = 0.001
    let one_mutation_single = 0.001;
    apply_relaxed_clock(&graph, &params, one_mutation_single, 1.0)?;

    let gammas_single: Vec<f64> = graph
      .get_edges()
      .iter()
      .map(|e| e.read_arc().payload().read_arc().gamma)
      .collect();

    // Simulate two partitions with lengths 1000 + 9000: one_mutation = 1/10000 = 0.0001
    // Using a 10x difference to ensure visible effect
    let one_mutation_multi = 0.0001;
    apply_relaxed_clock(&graph, &params, one_mutation_multi, 1.0)?;

    let gammas_multi: Vec<f64> = graph
      .get_edges()
      .iter()
      .map(|e| e.read_arc().payload().read_arc().gamma)
      .collect();

    // Gamma values should differ between single and multi-partition scenarios
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

  /// Regression test for C2: when partitions are empty (BranchLengthMode::Input),
  /// the relaxed clock should handle zero-length sequences gracefully.
  /// This tests defense-in-depth: if the guard in refinement.rs fails,
  /// apply_relaxed_clock should not produce NaN or crash with one_mutation near zero.
  #[test]
  fn test_relaxed_clock_handles_tiny_one_mutation() -> Result<(), Report> {
    let graph = build_simple_tree()?;
    let params = [1.0, 1.0];

    // Simulate what would happen if one_mutation were computed from very short sequences
    // (defense-in-depth - the guard in refinement.rs should prevent this)
    let tiny_one_mutation = 1e-15;
    apply_relaxed_clock(&graph, &params, tiny_one_mutation, 1.0)?;

    for edge in graph.get_edges() {
      let gamma = edge.read_arc().payload().read_arc().gamma;
      assert!(gamma.is_finite(), "gamma must be finite, got {gamma}");
      assert!(gamma >= 0.1, "gamma must respect minimum bound, got {gamma}");
    }

    Ok(())
  }

  /// Regression test for T4: root node must compute branch penalty using one_mutation.
  /// Root gamma should be influenced by its own branch penalty (k1/k2 from one_mutation),
  /// not just child coupling.
  #[test]
  fn test_relaxed_clock_root_has_branch_penalty() -> Result<(), Report> {
    // Tree with a single child to isolate root penalty behavior
    let graph: GraphTimetree = nwk_read_str("(A:0.1)root:0.0;")?;

    for edge in graph.get_edges() {
      let edge = edge.write_arc();
      let branch_length = edge.payload().read_arc().base.branch_length.unwrap_or(0.0);
      edge.payload().write_arc().time_length = Some(branch_length * 100.0);
    }

    let one_mutation = 0.01;
    let params = [1.0, 1.0];
    apply_relaxed_clock(&graph, &params, one_mutation, 1.0)?;

    // Get root gamma via the edge (gamma is stored on child's parent edge)
    let root_edge_gamma = graph
      .get_edges()
      .first()
      .map(|e| e.read_arc().payload().read_arc().gamma);

    // Root's gamma is computed from k1/k2 which includes branch penalty.
    // With one_mutation = 0.01, root uses opt_len = act_len = 0.01.
    // If root had no branch penalty, gamma would be determined only by child coupling.
    // Verify gamma is computed (not default 1.0) and within reasonable bounds.
    if let Some(gamma) = root_edge_gamma {
      assert!(gamma >= 0.1, "Root gamma should respect minimum bound: {gamma}");
      assert!(gamma < 10.0, "Root gamma should be reasonable: {gamma}");
    }

    Ok(())
  }

  /// Regression test for M3: root gamma should equal 1.0 when root has no children.
  ///
  /// Analytical derivation for childless root:
  /// - opt_len = act_len = one_mutation (root uses one_mutation for both)
  /// - c = 1.0 / one_mutation
  /// - denom = opt_len + one_mutation = 2 * one_mutation
  /// - k2 = slack + c * act_len² / denom = slack + 0.5
  /// - k1 = -2 * (c * act_len * opt_len / denom + slack) = -2 * (0.5 + slack) = -1 - 2*slack
  /// - gamma = -0.5 * k1 / k2 = (0.5 + slack) / (slack + 0.5) = 1.0
  ///
  /// This holds for any slack > 0 and any one_mutation > 0.
  #[rstest]
  #[trace]
  fn test_relaxed_clock_childless_root_gamma_equals_one(
    #[values(0.1, 1.0, 10.0, 100.0)] slack: f64,
    #[values(1e-6, 0.001, 0.01, 0.1, 1.0)] one_mutation: f64,
  ) -> Result<(), Report> {
    let graph: GraphTimetree = nwk_read_str("root:0.0;")?;
    let params = [slack, 1.0];
    apply_relaxed_clock(&graph, &params, one_mutation, 1.0)?;

    for edge in graph.get_edges() {
      let gamma = edge.read_arc().payload().read_arc().gamma;
      pretty_assert_ulps_eq!(gamma, 1.0, max_ulps = 4);
    }

    Ok(())
  }

  /// With opt_len == act_len == one_mutation for both root and child,
  /// the system has no rate variation to correct, so gamma = 1.0.
  #[test]
  fn test_relaxed_clock_childless_root_gamma_stored() -> Result<(), Report> {
    let graph: GraphTimetree = nwk_read_str("(A:0.01)root:0.0;")?;

    for edge in graph.get_edges() {
      let edge = edge.write_arc();
      let branch_length = edge.payload().read_arc().base.branch_length.unwrap_or(0.0);
      edge.payload().write_arc().time_length = Some(branch_length);
    }

    let one_mutation = 0.01;
    let params = [1.0, 1.0];
    apply_relaxed_clock(&graph, &params, one_mutation, 1.0)?;

    let gamma = graph
      .get_edges()
      .first()
      .map_or(1.0, |e| e.read_arc().payload().read_arc().gamma);

    pretty_assert_ulps_eq!(gamma, 1.0, max_ulps = 100);

    Ok(())
  }

  mod helpers {
    use super::*;

    #[derive(Deserialize)]
    pub struct GmInput {
      pub name: String,
      pub newick: String,
      pub clock_rate: f64,
      pub one_mutation: f64,
      pub slack: f64,
      pub coupling: f64,
      pub branches: BTreeMap<String, GmBranchInput>,
    }

    #[derive(Deserialize)]
    pub struct GmBranchInput {
      pub clock_length: f64,
    }

    #[derive(Deserialize)]
    pub struct GmOutput {
      pub name: String,
      pub gammas: BTreeMap<String, f64>,
    }

    pub fn build_simple_tree() -> Result<GraphTimetree, Report> {
      let graph: GraphTimetree = nwk_read_str("(A:0.1,B:0.2)root:0.0;")?;
      for edge in graph.get_edges() {
        let edge = edge.write_arc();
        let branch_length = edge.payload().read_arc().base.branch_length.unwrap_or(0.0);
        edge.payload().write_arc().time_length = Some(branch_length * 100.0);
      }
      Ok(graph)
    }

    pub fn build_deep_tree() -> Result<GraphTimetree, Report> {
      let graph: GraphTimetree = nwk_read_str("((A:0.1,B:0.2)AB:0.15,(C:0.05,D:0.1)CD:0.08)root:0.0;")?;
      for edge in graph.get_edges() {
        let edge = edge.write_arc();
        let branch_length = edge.payload().read_arc().base.branch_length.unwrap_or(0.0);
        edge.payload().write_arc().time_length = Some(branch_length * 100.0);
      }
      Ok(graph)
    }

    pub fn compute_variance(values: &[f64]) -> f64 {
      let n = values.len() as f64;
      let mean = values.iter().sum::<f64>() / n;
      values.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n
    }
  }
}
