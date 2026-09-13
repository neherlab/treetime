use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;

/// Collapse a single edge, updating graph topology and partition observations.
///
/// Graph topology: the target node is removed and its children become children of the
/// source node. Edge keys of the former children are preserved (graph reuses them).
///
/// Branch lengths: the collapsed edge's branch length is summed into each former-child
/// edge's branch length. Composed with substitutions, this preserves total evolutionary
/// distance from source to descendants.
///
/// Sparse partitions: substitutions on the collapsed edge are composed with each former
/// child edge's substitutions using the Markov semigroup property. Indels are
/// composed (overlapping/adjacent deletions merged, cancellations applied).
/// See `compose_substitutions()` and `compose_indels()` for details.
///
/// Stale observations for the removed node and removed edge are dropped from every sparse
/// partition so the observation maps stay consistent with the graph. The evolving node
/// states, edge messages, and estimates are reconciled by the caller after the topology
/// batch and rebuilt by the next marginal update.
///
/// # Scientific background
///
/// Edge collapse is a graph contraction operation. When the collapsed edge carries
/// substitutions, they are composed with each child edge's substitutions via the
/// Markov semigroup property rather than taking the set union. Reversions at the
/// same position (e.g. a forward substitution on the collapsed edge followed by the
/// reverse substitution on the child edge) cancel to no net change. See
/// [`SparseEdgeObs::chain_fitch_subs`](crate::partition::storage::sparse::SparseEdgeObs::chain_fitch_subs)
/// for the exact composition semantics.
pub fn collapse_edge(
  graph: &mut Graph,
  sparse: &mut [PartitionMarginalSparse],
  edge_key: GraphEdgeKey,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<(), Report> {
  let target_node_key = graph.get_target_node_key(edge_key)?;

  let removed_bl = branch_lengths[&edge_key];
  let (_, _removed_edge, new_edges) = graph.collapse_edge(edge_key)?;

  for new_edge in &new_edges {
    let new_edge_key = new_edge.read_arc().key();

    // Sum branch lengths: net edge length = collapsed-edge length + child-edge length.
    // Both must be present for a sum; a missing weight (`None`) on either side is preserved.
    if let (Some(bl1), Some(bl2)) = (removed_bl, branch_lengths[&new_edge_key]) {
      branch_lengths.insert(new_edge_key, Some(bl1 + bl2));
    }

    // Compose substitutions and merge indels on each sparse partition's observations.
    // Each graph edge is expected to have a corresponding observation entry (populated
    // by the Fitch pre-pass); strict indexing surfaces any invariant violation.
    for partition in sparse.iter_mut() {
      let obs_edges = &mut partition.obs_edges;
      let removed_edge_data = obs_edges[&edge_key].clone();
      let child_edge = obs_edges.entry(new_edge_key).or_default();
      let merged_subs = removed_edge_data.chain_fitch_subs(child_edge.fitch_subs())?;
      child_edge.set_fitch_subs(merged_subs);
      child_edge.indels = removed_edge_data.chain_fitch_indels(&child_edge.indels);
    }
  }

  // Drop stale observations for the removed node and removed edge in every partition.
  // Downstream passes (e.g. `marginal_update`) recompute any state they need from
  // the remaining entries.
  for partition in sparse.iter_mut() {
    partition.obs_nodes.remove(&target_node_key);
    partition.obs_edges.remove(&edge_key);
  }

  // The collapsed edge no longer exists; drop its stale entry from the branch-length map.
  branch_lengths.remove(&edge_key);

  Ok(())
}
