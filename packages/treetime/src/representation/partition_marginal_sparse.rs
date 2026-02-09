use crate::alphabet::alphabet::Alphabet;
use crate::commands::clock::reroot::{EdgeMergeInfo, EdgeSplitInfo};
use crate::commands::timetree::partition_ops::PartitionRerootOps;
use crate::graph::edge::{EdgeOptimizeOps, GraphEdgeKey};
use crate::graph::graph::Graph;
use crate::graph::graph_traverse::{GraphNodeBackward, GraphNodeForward};
use crate::graph::node::{GraphNode, GraphNodeKey, Named};
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::PartitionWithGtrInference;
use crate::io::fasta::FastaRecord;
use crate::make_internal_report;
use crate::representation::graph_ancestral::GraphAncestral;
use crate::representation::graph_sparse::{MarginalSparseSeqDistribution, SparseEdgePartition, SparseNodePartition};
use crate::representation::log_lh::HasLogLh;
use crate::representation::partition_compressed::PartitionCompressed;
use crate::representation::partition_marginal::{PartitionMarginal, PartitionMarginalOps};
use crate::representation::partition_marginal_sparse_passes;
use crate::representation::seq::Seq;
use crate::seq::composition::Composition;
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::Itertools;
use ndarray_stats::QuantileExt;
use std::collections::BTreeMap;
use std::mem;
use treetime_utils::container::get_exactly_one;

#[derive(Clone, Debug)]
pub struct PartitionMarginalSparse {
  pub index: usize,
  pub gtr: GTR,
  pub alphabet: Alphabet,
  pub length: usize,
  pub nodes: BTreeMap<GraphNodeKey, SparseNodePartition>,
  pub edges: BTreeMap<GraphEdgeKey, SparseEdgePartition>,
}

impl PartitionCompressed for PartitionMarginalSparse {
  fn index(&self) -> usize {
    self.index
  }

  fn alphabet(&self) -> &Alphabet {
    &self.alphabet
  }

  fn length(&self) -> usize {
    self.length
  }

  fn nodes(&self) -> &BTreeMap<GraphNodeKey, SparseNodePartition> {
    &self.nodes
  }

  fn edges(&self) -> &BTreeMap<GraphEdgeKey, SparseEdgePartition> {
    &self.edges
  }

  fn nodes_mut(&mut self) -> &mut BTreeMap<GraphNodeKey, SparseNodePartition> {
    &mut self.nodes
  }

  fn edges_mut(&mut self) -> &mut BTreeMap<GraphEdgeKey, SparseEdgePartition> {
    &mut self.edges
  }
}

impl PartitionMarginalSparse {
  #[allow(clippy::same_name_method)]
  pub fn get_sequence_length(&self) -> usize {
    self.length
  }
}

impl HasLogLh for PartitionMarginalSparse {
  fn get_log_lh(&self, node_key: GraphNodeKey) -> f64 {
    self.nodes.get(&node_key).map_or(0.0, |node| node.profile.log_lh)
  }
}

impl PartitionMarginal for PartitionMarginalSparse {}

impl PartitionRerootOps for PartitionMarginalSparse {}

impl<N, E> crate::commands::timetree::partition_ops::PartitionTimetreeOps<N, E> for PartitionMarginalSparse
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  fn create_edge_contribution(
    &self,
    edge_key: GraphEdgeKey,
  ) -> Result<crate::commands::optimize::optimize_unified::OptimizationContribution, Report> {
    crate::commands::optimize::optimize_unified::OptimizationContribution::from_sparse(edge_key, self)
  }

  fn handle_edge_split(&mut self, info: &EdgeSplitInfo) -> Result<(), Report> {
    // Create new node entry with empty placeholder (sequence computed during marginal update)
    self
      .nodes
      .insert(info.new_node_key, SparseNodePartition::empty(&self.alphabet));

    // Move mutations from old edge to child-side edge, parent-side edge is empty
    let old_edge_data = self
      .edges
      .remove(&info.old_edge_key)
      .ok_or_else(|| make_internal_report!("Old edge {:?} must exist for split", info.old_edge_key))?;

    // Child-side edge gets all mutations (they describe parent->child relationship)
    self.edges.insert(info.child_side_edge_key, old_edge_data);

    // Parent-side edge is empty (no mutations between parent and new split node)
    self
      .edges
      .insert(info.parent_side_edge_key, SparseEdgePartition::default());

    Ok(())
  }

  fn handle_edge_merge(&mut self, info: &EdgeMergeInfo) -> Result<(), Report> {
    // Remove the trivial node
    self.nodes.remove(&info.removed_node_key);

    // Get mutations from both edges
    let parent_edge = self
      .edges
      .remove(&info.parent_edge_key)
      .ok_or_else(|| make_internal_report!("Parent edge {:?} must exist for merge", info.parent_edge_key))?;
    let child_edge = self
      .edges
      .remove(&info.child_edge_key)
      .ok_or_else(|| make_internal_report!("Child edge {:?} must exist for merge", info.child_edge_key))?;

    // Compose substitutions: parent then child
    let merged_subs = compose_substitutions(&parent_edge.subs, &child_edge.subs)?;

    // Compose indels: concatenate (parent indels first, then child indels)
    let mut merged_indels = parent_edge.indels;
    merged_indels.extend(child_edge.indels);

    // Create merged edge with composed mutations
    let merged_edge = SparseEdgePartition {
      subs: merged_subs,
      indels: merged_indels,
      ..SparseEdgePartition::default()
    };

    self.edges.insert(info.merged_edge_key, merged_edge);

    Ok(())
  }

  fn update_partition_after_reroot(
    &mut self,
    _old_root_key: GraphNodeKey,
    _new_root_key: GraphNodeKey,
    path_from_old_to_new: &[(GraphNodeKey, Option<GraphEdgeKey>)],
  ) -> Result<(), Report> {
    // Build ordered list of edge keys on the reroot path (old root -> new root direction)
    let reroot_edge_keys = path_from_old_to_new
      .iter()
      .filter_map(|(_, edge_key)| *edge_key)
      .collect_vec();

    // Node sequences are not computed here. The new root node was created by handle_edge_split
    // with an empty placeholder. Sequences are recomputed by the subsequent marginal update pass
    // (process_node_backward + process_node_forward), which overwrites any sequence data.
    //
    // This method only needs to invert edge mutations so they reflect the new parent->child
    // direction after reroot.

    // Invert edge data for each edge on the reroot path
    for edge_key in &reroot_edge_keys {
      let edge_data = self
        .edges
        .get_mut(edge_key)
        .ok_or_else(|| make_internal_report!("Edge {edge_key:?} must exist on reroot path"))?;

      // Invert all substitutions
      for sub in &mut edge_data.subs {
        sub.invert();
      }

      // Invert all indels
      for indel in &mut edge_data.indels {
        indel.invert();
      }

      // Swap directional messages
      mem::swap(&mut edge_data.msg_to_parent, &mut edge_data.msg_to_child);

      // Reset msg_from_child (stale propagated message cache)
      edge_data.msg_from_child = MarginalSparseSeqDistribution::default();
    }

    Ok(())
  }
}

impl<N, E> PartitionMarginalOps<N, E> for PartitionMarginalSparse
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  fn attach_sequences(&mut self, _graph: &Graph<N, E, ()>, _aln: &[FastaRecord]) -> Result<(), Report> {
    // Sparse partitions get sequences attached during compression phase
    Ok(())
  }

  fn process_node_backward(&mut self, node: &GraphNodeBackward<N, E, ()>) -> Result<(), Report> {
    partition_marginal_sparse_passes::process_node_backward(self, node)
  }

  fn process_node_forward(&mut self, graph: &Graph<N, E, ()>, node: &GraphNodeForward<N, E, ()>) -> Result<(), Report> {
    partition_marginal_sparse_passes::process_node_forward(self, graph, node)
  }

  fn extract_ancestral_sequence(&mut self, node_key: GraphNodeKey) -> Seq {
    if let Some(node_data) = self.nodes.get(&node_key) {
      if !node_data.seq.sequence.is_empty() {
        node_data.seq.sequence.clone()
      } else {
        // Return empty seq to be filled later by the reconstruction algorithm
        crate::seq![]
      }
    } else {
      crate::seq![]
    }
  }

  fn reconstruct_node_sequence(&mut self, node: &GraphNodeForward<N, E, ()>, include_leaves: bool) -> Option<Seq> {
    if !include_leaves && node.is_leaf {
      return None;
    }

    let mut node_data = self.nodes.remove(&node.key)?;

    let mut seq = if node.is_root {
      node_data.seq.sequence.clone()
    } else {
      let (parent_key, edge_key) = get_exactly_one(&node.parent_keys).ok()?;
      let parent_data = self.nodes.get(parent_key)?;
      let edge_data = self.edges.get(edge_key)?;

      let mut seq = parent_data.seq.sequence.clone();

      // Implant mutations
      for m in &edge_data.subs {
        seq[m.pos()] = m.qry();
      }

      // Implant indels
      for indel in &edge_data.indels {
        if indel.deletion {
          seq[indel.range.0..indel.range.1].fill(self.alphabet.gap());
        } else {
          seq[indel.range.0..indel.range.1].copy_from_slice(&indel.seq);
        }
      }

      seq
    };

    // At the node itself, mask whatever is unknown in the node.
    let alphabet = &self.alphabet;
    for r in &node_data.seq.unknown {
      let ambig_char = alphabet.unknown();
      seq[r.0..r.1].fill(ambig_char);
    }

    // change variable sites to their most likely state
    for (pos, states) in &node_data.profile.variable {
      seq[*pos] = alphabet.char(states.dis.argmax().unwrap());
    }

    node_data.seq.sequence = seq.clone();
    self.nodes.insert(node.key, node_data);

    Some(seq)
  }

  fn get_sequence_length(&self) -> usize {
    self.length
  }
}

impl PartitionWithGtrInference for PartitionMarginalSparse {
  fn alphabet(&self) -> &Alphabet {
    &self.alphabet
  }

  fn get_seq_composition(&self, node_key: GraphNodeKey) -> &Composition {
    &self.nodes[&node_key].seq.composition
  }

  fn get_edge_substitutions(&self, edge_key: GraphEdgeKey, _graph: &GraphAncestral) -> Vec<Sub> {
    self.edges[&edge_key].subs.clone()
  }
}

/// Compose substitutions from two edges: parent edge then child edge.
/// Applies cancellation rules: A5G + G5T = A5T, A5G + G5A = (no mutation).
fn compose_substitutions(parent_subs: &[Sub], child_subs: &[Sub]) -> Result<Vec<Sub>, Report> {
  // Index child subs by position for efficient lookup
  let child_by_pos: BTreeMap<usize, &Sub> = child_subs.iter().map(|s| (s.pos(), s)).collect();

  let mut result = Vec::with_capacity(parent_subs.len() + child_subs.len());

  // Process parent subs - compose with child if position overlaps
  for parent_sub in parent_subs {
    let pos = parent_sub.pos();
    if let Some(child_sub) = child_by_pos.get(&pos) {
      // Composition: parent reff -> parent qry -> child qry
      // Net effect: parent reff -> child qry
      if parent_sub.reff() != child_sub.qry() {
        // Non-cancelling: A5G + G5T = A5T
        result.push(Sub::new(parent_sub.reff(), pos, child_sub.qry())?);
      }
      // else: Cancelling mutation A5G + G5A = no net change
    } else {
      // No child mutation at this position - keep parent mutation
      result.push(parent_sub.clone());
    }
  }

  // Add child subs at positions not covered by parent
  let parent_by_pos: BTreeMap<usize, &Sub> = parent_subs.iter().map(|s| (s.pos(), s)).collect();
  for child_sub in child_subs {
    if !parent_by_pos.contains_key(&child_sub.pos()) {
      result.push(child_sub.clone());
    }
  }

  Ok(result)
}

#[cfg(test)]
mod tests {
  use super::*;
  use eyre::Result;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  fn sub(reff: char, pos: usize, qry: char) -> Sub {
    Sub::new(reff, pos, qry).unwrap()
  }

  #[test]
  fn test_compose_substitutions_empty_both() -> Result<()> {
    let result = compose_substitutions(&[], &[])?;
    assert_eq!(Vec::<Sub>::new(), result);
    Ok(())
  }

  #[test]
  fn test_compose_substitutions_empty_parent() -> Result<()> {
    let child = vec![sub('A', 5, 'G'), sub('C', 10, 'T')];
    let result = compose_substitutions(&[], &child)?;
    assert_eq!(child, result);
    Ok(())
  }

  #[test]
  fn test_compose_substitutions_empty_child() -> Result<()> {
    let parent = vec![sub('A', 5, 'G'), sub('C', 10, 'T')];
    let result = compose_substitutions(&parent, &[])?;
    assert_eq!(parent, result);
    Ok(())
  }

  #[test]
  fn test_compose_substitutions_no_overlap() -> Result<()> {
    let parent = vec![sub('A', 5, 'G')];
    let child = vec![sub('C', 10, 'T')];
    let result = compose_substitutions(&parent, &child)?;
    assert_eq!(vec![sub('A', 5, 'G'), sub('C', 10, 'T')], result);
    Ok(())
  }

  #[test]
  fn test_compose_substitutions_chain() -> Result<()> {
    // A5G + G5T = A5T (composition at same position)
    let parent = vec![sub('A', 5, 'G')];
    let child = vec![sub('G', 5, 'T')];
    let result = compose_substitutions(&parent, &child)?;
    assert_eq!(vec![sub('A', 5, 'T')], result);
    Ok(())
  }

  #[test]
  fn test_compose_substitutions_cancellation() -> Result<()> {
    // A5G + G5A = no mutation (cancellation)
    let parent = vec![sub('A', 5, 'G')];
    let child = vec![sub('G', 5, 'A')];
    let result = compose_substitutions(&parent, &child)?;
    assert_eq!(Vec::<Sub>::new(), result);
    Ok(())
  }

  #[test]
  fn test_compose_substitutions_mixed() -> Result<()> {
    // Position 5: A->G + G->T = A->T (chain)
    // Position 10: C->T (parent only)
    // Position 15: G->A (child only)
    // Position 20: T->C + C->T = no mutation (cancellation)
    let parent = vec![sub('A', 5, 'G'), sub('C', 10, 'T'), sub('T', 20, 'C')];
    let child = vec![sub('G', 5, 'T'), sub('G', 15, 'A'), sub('C', 20, 'T')];
    let result = compose_substitutions(&parent, &child)?;
    // Expected: pos 5 chained, pos 10 kept, pos 15 added, pos 20 cancelled
    assert_eq!(vec![sub('A', 5, 'T'), sub('C', 10, 'T'), sub('G', 15, 'A')], result);
    Ok(())
  }

  #[rstest]
  #[case('A', 'G', 'T', Some(('A', 'T')))] // chain: A->G->T = A->T
  #[case('A', 'G', 'A', None)] // cancel: A->G->A = none
  #[case('C', 'T', 'G', Some(('C', 'G')))] // chain: C->T->G = C->G
  #[case('G', 'A', 'G', None)] // cancel: G->A->G = none
  fn test_compose_substitutions_single_position(
    #[case] parent_reff: char,
    #[case] intermediate: char,
    #[case] child_qry: char,
    #[case] expected: Option<(char, char)>,
  ) -> Result<()> {
    let parent = vec![sub(parent_reff, 0, intermediate)];
    let child = vec![sub(intermediate, 0, child_qry)];
    let result = compose_substitutions(&parent, &child)?;
    let expected = expected.map_or_else(Vec::new, |(r, q)| vec![sub(r, 0, q)]);
    assert_eq!(expected, result);
    Ok(())
  }
}
