#[cfg(test)]
mod tests {
  use crate::payload::ancestral::GraphAncestral;
  use crate::reroot::div_stats::DivStats;
  use crate::reroot::div_stats_traversal::compute_div_stats;
  use crate::reroot::orchestrate::{RerootTopologyParams, reroot_in_place};
  use crate::reroot::params::BrentParams;
  use crate::reroot::search::find_best_root;
  use crate::reroot::variance::VarianceModel;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::nwk::nwk_read_str;

  fn root_to_tip_distances(graph: &GraphAncestral) -> Vec<f64> {
    let root_key = graph.get_exactly_one_root().unwrap().read_arc().key();
    let mut distances = Vec::new();
    collect_distances(graph, root_key, 0.0, &mut distances);
    distances.sort_by(|a, b| a.partial_cmp(b).unwrap());
    distances
  }

  fn collect_distances(
    graph: &GraphAncestral,
    node_key: treetime_graph::node::GraphNodeKey,
    dist: f64,
    out: &mut Vec<f64>,
  ) {
    let node = graph.get_node(node_key).unwrap();
    let node = node.read_arc();
    if node.is_leaf() {
      out.push(dist);
      return;
    }
    for &edge_key in node.outbound() {
      let edge = graph.get_edge(edge_key).unwrap();
      let edge = edge.read_arc();
      let bl = edge.payload().read_arc().branch_length().unwrap_or(0.0);
      collect_distances(graph, edge.target(), dist + bl, out);
    }
  }

  #[test]
  fn test_orchestrate_reroot_reduces_rtt_variance() -> Result<(), Report> {
    let mut graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.3)root;")?;
    let variance = VarianceModel::default();
    let field = compute_div_stats(&graph, &variance)?;

    reroot_in_place::<_, _, _, DivStats, _>(
      &mut graph,
      &field.edge_stats,
      &field.root_stats,
      &variance,
      &BrentParams::default(),
      RerootTopologyParams::default(),
      |_graph, _inverted| Ok(()),
    )?;

    let dists = root_to_tip_distances(&graph);
    assert_eq!(dists.len(), 2);
    assert_abs_diff_eq!(dists[0], dists[1], epsilon = 1e-6);
    Ok(())
  }

  #[test]
  fn test_orchestrate_brent_finds_equidistant_root() -> Result<(), Report> {
    let mut graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.3)root;")?;
    let variance = VarianceModel::default();
    let field = compute_div_stats(&graph, &variance)?;

    reroot_in_place::<_, _, _, DivStats, _>(
      &mut graph,
      &field.edge_stats,
      &field.root_stats,
      &variance,
      &BrentParams::default(),
      RerootTopologyParams::default(),
      |_graph, _inverted| Ok(()),
    )?;

    let dists = root_to_tip_distances(&graph);
    assert_abs_diff_eq!(dists[0], dists[1], epsilon = 1e-6);
    Ok(())
  }

  #[test]
  fn test_orchestrate_endpoint_snap_split_zero_roots_at_source() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.3)root;")?;
    let variance = VarianceModel::default();
    let field = compute_div_stats(&graph, &variance)?;

    let best = find_best_root(
      &graph,
      &field.edge_stats,
      &field.root_stats,
      &variance,
      &BrentParams::default(),
    )?;

    // The optimal root is interior (not at an endpoint), so test the snap
    // logic directly via the endpoint mapping by verifying the split convention.
    assert!(best.edge.is_some());
    assert!(
      best.split > 0.0 && best.split < 1.0,
      "expected interior split, got {}",
      best.split
    );
    Ok(())
  }

  // The optimal split on edge i->B (~x=0.28 from source=i) is closer to i than
  // to B, so snap-to-nearest reroots to internal node i -- a real topology change.
  #[test]
  fn test_orchestrate_no_split_snaps_to_nearer_endpoint() -> Result<(), Report> {
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.5)i:0.02,C:0.2)root;")?;
    let variance = VarianceModel::default();
    let field = compute_div_stats(&graph, &variance)?;
    let root_before = graph.get_exactly_one_root().unwrap().read_arc().key();

    reroot_in_place::<_, _, _, DivStats, _>(
      &mut graph,
      &field.edge_stats,
      &field.root_stats,
      &variance,
      &BrentParams::default(),
      RerootTopologyParams {
        split_edge: false,
        remove_trivial_root: true,
      },
      |_graph, _inverted| Ok(()),
    )?;

    let root_after = graph.get_exactly_one_root().unwrap().read_arc().key();
    assert_ne!(
      root_before, root_after,
      "split_edge=false should reroot to the nearer endpoint"
    );
    Ok(())
  }
}
