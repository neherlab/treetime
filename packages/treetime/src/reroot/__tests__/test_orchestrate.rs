#[cfg(test)]
mod tests {
  use crate::progress::NoopProgress;
  use crate::reroot::orchestrate::{RerootTopologyParams, reroot_min_dev};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::nwk_read;

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
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.3)root;".as_slice())?;
    let names = nwk_parsed.names();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;

    reroot_min_dev(
      &mut graph,
      RerootTopologyParams::default(),
      &mut branch_lengths,
      &names,
      &NoopProgress,
    )?;

    let dists = root_to_tip_distances(&graph, &branch_lengths);
    assert_eq!(2, dists.len());
    assert_abs_diff_eq!(dists[0], 0.2, epsilon = 1e-15);
    assert_abs_diff_eq!(dists[1], 0.2, epsilon = 1e-15);
    Ok(())
  }

  // Root at distance x up the branch to C of the star (A:1, B:2, C:6): 3 Var = 14 - 12 x + 8 x^2 / 3, minimal at x = 9/4
  // with distances 3.25, 4.25, 3.75; on the branches to A and B the variance only grows.
  #[test]
  fn test_orchestrate_reroot_min_dev_places_the_root_at_the_analytic_minimum_of_the_variance() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:1,B:2,C:6)R;".as_slice())?;
    let names = nwk_parsed.names();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;

    reroot_min_dev(
      &mut graph,
      RerootTopologyParams::default(),
      &mut branch_lengths,
      &names,
      &NoopProgress,
    )?;

    let dists = root_to_tip_distances(&graph, &branch_lengths);
    assert_abs_diff_eq!(dists.as_slice(), [3.25, 3.75, 4.25].as_slice(), epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_orchestrate_reroot_min_dev_removes_an_undated_stem_before_the_root_search() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(((A:0.1,B:0.3)X:0.05,C:0.4)R:0.7)STEM;".as_slice())?;
    let names = nwk_parsed.names();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let unstemmed = nwk_read(b"((A:0.1,B:0.3)X:0.05,C:0.4)R;".as_slice())?;
    let unstemmed_names = unstemmed.names();
    let mut unstemmed_graph = unstemmed.graph;
    let mut unstemmed_branch_lengths = unstemmed.branch_lengths;

    let result = reroot_min_dev(
      &mut graph,
      RerootTopologyParams::default(),
      &mut branch_lengths,
      &names,
      &NoopProgress,
    )?;
    reroot_min_dev(
      &mut unstemmed_graph,
      RerootTopologyParams::default(),
      &mut unstemmed_branch_lengths,
      &unstemmed_names,
      &NoopProgress,
    )?;

    assert!(result.stem_removal.is_some());
    assert_eq!(3, graph.get_leaves().count());
    let expected = root_to_tip_distances(&unstemmed_graph, &unstemmed_branch_lengths);
    let actual = root_to_tip_distances(&graph, &branch_lengths);
    assert_abs_diff_eq!(expected.as_slice(), actual.as_slice(), epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_orchestrate_no_split_snaps_to_nearer_endpoint() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:0.1,B:0.5)i:0.02,C:0.2)root;".as_slice())?;
    let names = nwk_parsed.names();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let root_before = graph.get_exactly_one_root().unwrap().key();

    reroot_min_dev(
      &mut graph,
      RerootTopologyParams {
        split_edge: false,
        remove_trivial_root: true,
      },
      &mut branch_lengths,
      &names,
      &NoopProgress,
    )?;

    let root_after = graph.get_exactly_one_root().unwrap().key();
    assert_ne!(
      root_before, root_after,
      "split_edge=false should reroot to the nearer endpoint"
    );
    Ok(())
  }
}
