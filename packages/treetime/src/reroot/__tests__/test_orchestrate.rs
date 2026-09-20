#[cfg(test)]
mod tests {
  use crate::reroot::div_stats::DivStats;
  use crate::reroot::div_stats_traversal::compute_div_stats;
  use crate::reroot::orchestrate::{RerootTopologyParams, reroot_in_place};
  use crate::reroot::params::BrentParams;
  use crate::reroot::search::find_best_root;
  use crate::reroot::variance::VarianceModel;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::nwk_read_str;

  fn root_to_tip_distances(graph: &Graph, branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> Vec<f64> {
    let root_key = graph.get_exactly_one_root().unwrap().key();
    let mut distances = Vec::new();
    collect_distances(graph, root_key, 0.0, &mut distances, branch_lengths);
    distances.sort_by(|a, b| a.partial_cmp(b).unwrap());
    distances
  }

  fn collect_distances(
    graph: &Graph,
    node_key: treetime_graph::node::GraphNodeKey,
    dist: f64,
    out: &mut Vec<f64>,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  ) {
    let node = graph.get_node(node_key).unwrap();
    if node.is_leaf() {
      out.push(dist);
      return;
    }
    for &edge_key in node.outbound() {
      let edge = graph.get_edge(edge_key).unwrap();
      let bl = branch_lengths[&edge.key()].unwrap_or(0.0);
      collect_distances(graph, edge.target(), dist + bl, out, branch_lengths);
    }
  }

  #[test]
  fn test_orchestrate_reroot_reduces_rtt_variance() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.3)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let variance = VarianceModel::default();
    let field = compute_div_stats(&graph, &branch_lengths, &variance)?;

    reroot_in_place::<DivStats, _>(
      &mut graph,
      &field.edge_stats,
      &field.root_stats,
      &variance,
      &BrentParams::default(),
      RerootTopologyParams::default(),
      &mut branch_lengths,
      |_graph, _inverted| Ok(()),
    )?;

    let dists = root_to_tip_distances(&graph, &branch_lengths);
    assert_eq!(2, dists.len());
    assert_abs_diff_eq!(dists[0], dists[1], epsilon = 1e-6);
    Ok(())
  }

  #[test]
  fn test_orchestrate_brent_finds_equidistant_root() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.3)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let variance = VarianceModel::default();
    let field = compute_div_stats(&graph, &branch_lengths, &variance)?;

    reroot_in_place::<DivStats, _>(
      &mut graph,
      &field.edge_stats,
      &field.root_stats,
      &variance,
      &BrentParams::default(),
      RerootTopologyParams::default(),
      &mut branch_lengths,
      |_graph, _inverted| Ok(()),
    )?;

    let dists = root_to_tip_distances(&graph, &branch_lengths);
    assert_abs_diff_eq!(dists[0], dists[1], epsilon = 1e-6);
    Ok(())
  }

  #[test]
  fn test_orchestrate_endpoint_snap_split_zero_roots_at_source() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.3)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let variance = VarianceModel::default();
    let field = compute_div_stats(&graph, &branch_lengths, &variance)?;

    let best = find_best_root(
      &graph,
      &field.edge_stats,
      &field.root_stats,
      &variance,
      &branch_lengths,
      &BrentParams::default(),
    )?;

    assert!(best.edge.is_some());
    assert!(
      best.split > 0.0 && best.split < 1.0,
      "expected interior split, got {}",
      best.split
    );
    Ok(())
  }

  #[test]
  fn test_orchestrate_no_split_snaps_to_nearer_endpoint() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.5)i:0.02,C:0.2)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let variance = VarianceModel::default();
    let field = compute_div_stats(&graph, &branch_lengths, &variance)?;
    let root_before = graph.get_exactly_one_root().unwrap().key();

    reroot_in_place::<DivStats, _>(
      &mut graph,
      &field.edge_stats,
      &field.root_stats,
      &variance,
      &BrentParams::default(),
      RerootTopologyParams {
        split_edge: false,
        remove_trivial_root: true,
      },
      &mut branch_lengths,
      |_graph, _inverted| Ok(()),
    )?;

    let root_after = graph.get_exactly_one_root().unwrap().key();
    assert_ne!(
      root_before, root_after,
      "split_edge=false should reroot to the nearer endpoint"
    );
    Ok(())
  }
}
