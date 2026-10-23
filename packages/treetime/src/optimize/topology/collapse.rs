use crate::partition::assembly::AlignmentGroup;
use crate::partition::marginal::reconstruction::SparseReconstruction;
use crate::partition::marginal::sparse::alignment_indels::{composed_edge_pieces, store_edge_pieces};
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;

pub(crate) fn collapse_edge(
  graph: &mut Graph,
  sparse: &mut [SparseReconstruction],
  groups: &[AlignmentGroup],
  edge_key: GraphEdgeKey,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<(), Report> {
  let target_node_key = graph.get_target_node_key(edge_key)?;

  let removed_bl = branch_lengths[&edge_key];
  let (_, _removed_edge, new_edges) = graph.collapse_edge(edge_key)?;

  for new_edge in &new_edges {
    let new_edge_key = *new_edge;

    if let (Some(bl1), Some(bl2)) = (removed_bl, branch_lengths[&new_edge_key]) {
      branch_lengths.insert(new_edge_key, Some(bl1 + bl2));
    }

    for reconstruction in sparse.iter_mut() {
      let obs_edges = &mut reconstruction.partition.obs_edges;
      let removed_edge_data = obs_edges[&edge_key].clone();
      let child_edge = obs_edges.entry(new_edge_key).or_default();
      let merged_subs = removed_edge_data.chain_fitch_subs(child_edge.fitch_subs())?;
      child_edge.set_fitch_subs(merged_subs);
    }
    for group in groups {
      let pieces = composed_edge_pieces(group, sparse, edge_key, new_edge_key)?;
      store_edge_pieces(group, sparse, new_edge_key, pieces);
    }
  }

  for reconstruction in sparse.iter_mut() {
    reconstruction.partition.obs_nodes.remove(&target_node_key);
    reconstruction.partition.obs_edges.remove(&edge_key);
  }

  branch_lengths.remove(&edge_key);

  Ok(())
}
