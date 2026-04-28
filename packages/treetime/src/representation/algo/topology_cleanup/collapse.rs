use crate::representation::partition::marginal_dense::PartitionMarginalDense;
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::representation::payload::ancestral::GraphAncestral;
use eyre::Report;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::edge::{GraphEdgeKey, HasBranchLength};

/// Collapse a single edge, updating graph topology and partition state.
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
/// Stale entries for the removed node and removed edge are dropped from every sparse
/// and dense partition so the partition state stays consistent with the graph.
///
/// # Scientific background
///
/// Edge collapse is a graph contraction operation. When the collapsed edge carries
/// substitutions, they are composed with each child edge's substitutions via the
/// Markov semigroup property rather than taking the set union. Reversions at the
/// same position (e.g. a forward substitution on the collapsed edge followed by the
/// reverse substitution on the child edge) cancel to no net change. See
/// [`compose_substitutions`] for the exact composition semantics.
pub fn collapse_edge(
  graph: &mut GraphAncestral,
  sparse_partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  dense_partitions: &[Arc<RwLock<PartitionMarginalDense>>],
  edge_key: GraphEdgeKey,
) -> Result<(), Report> {
  let target_node_key = graph.get_target_node_key(edge_key)?;

  let (_, removed_edge, new_edges) = graph.collapse_edge(edge_key)?;
  let removed_edge = removed_edge.payload().read_arc();

  for new_edge in &new_edges {
    let new_edge_key = new_edge.read_arc().key();
    let mut new_edge_payload = new_edge.write_arc().payload().write_arc();

    // Sum branch lengths: net edge length = collapsed-edge length + child-edge length
    if let (Some(bl1), Some(bl2)) = (removed_edge.branch_length(), new_edge_payload.branch_length()) {
      new_edge_payload.set_branch_length(Some(bl1 + bl2));
    }

    // Compose substitutions and merge indels on each sparse partition.
    // Each graph edge is expected to have a corresponding partition entry (populated
    // by `compress_sequences`); strict indexing surfaces any invariant violation.
    for partition in sparse_partitions {
      let mut partition = partition.write_arc();
      let removed_edge_data = partition.edges[&edge_key].clone();
      let child_edge = partition.edges.entry(new_edge_key).or_default();
      let merged_subs = removed_edge_data.chain_fitch_subs(child_edge.fitch_subs())?;
      child_edge.set_fitch_subs(merged_subs);
      child_edge.indels = removed_edge_data.chain_fitch_indels(&child_edge.indels);
    }
  }

  // Drop stale entries for the removed node and removed edge in every partition.
  // Downstream passes (e.g. `update_marginal`) recompute any state they need from
  // the remaining entries.
  for partition in sparse_partitions {
    let mut partition = partition.write_arc();
    partition.nodes.remove(&target_node_key);
    partition.edges.remove(&edge_key);
  }
  for partition in dense_partitions {
    let mut partition = partition.write_arc();
    partition.nodes.remove(&target_node_key);
    partition.edges.remove(&edge_key);
  }

  Ok(())
}
