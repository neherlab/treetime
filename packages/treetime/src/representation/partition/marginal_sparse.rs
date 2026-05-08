use crate::alphabet::alphabet::Alphabet;
use crate::commands::optimize::optimize_unified::OptimizationContribution;
use crate::gtr::gtr::GTR;
use crate::representation::algo::topology_cleanup::reroot::RerootChanges;
use crate::representation::partition::marginal_passes;
use crate::representation::partition::traits::{
  BranchTopology, HasGtr, HasLogLh, PartitionBranchOps, PartitionMarginal, PartitionMarginalOps, PartitionOptimizeOps,
  PartitionRerootOps, PartitionTimetreeOps,
};
use crate::representation::payload::sparse::{MarginalSparseSeqDistribution, SparseEdgePartition, SparseNodePartition};
use crate::seq::mutation::Sub;
use crate::{make_error, make_internal_report};
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};
use std::mem;
use treetime_graph::edge::{EdgeOptimizeOps, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::{GraphNodeBackward, GraphNodeForward};
use treetime_graph::node::{GraphNode, GraphNodeKey, Named};
use treetime_io::fasta::FastaRecord;
use treetime_primitives::{Seq, seq};
use treetime_utils::array::ndarray::argmax_first;
use treetime_utils::collections::container::get_exactly_one;
use treetime_utils::interval::range_union::range_union;

#[derive(Clone, Debug)]
pub struct PartitionMarginalSparse {
  pub index: usize,
  pub gtr: GTR,
  pub alphabet: Alphabet,
  pub length: usize,
  pub root_sequence: Seq,
  pub nodes: BTreeMap<GraphNodeKey, SparseNodePartition>,
  pub edges: BTreeMap<GraphEdgeKey, SparseEdgePartition>,
}

impl HasGtr for PartitionMarginalSparse {
  fn gtr(&self) -> &GTR {
    &self.gtr
  }
  fn gtr_mut(&mut self) -> &mut GTR {
    &mut self.gtr
  }
  fn sequence_length(&self) -> usize {
    self.length
  }
}

impl PartitionMarginalSparse {
  #[allow(clippy::same_name_method)]
  pub fn get_sequence_length(&self) -> usize {
    self.length
  }
}

pub(crate) fn reconstruct_map_seq(
  base_seq: &Seq,
  edge: Option<&SparseEdgePartition>,
  node: &SparseNodePartition,
  alphabet: &Alphabet,
) -> Seq {
  let mut seq = if let Some(edge) = edge {
    let mut seq = base_seq.clone();
    for m in edge.fitch_subs() {
      seq[m.pos()] = m.qry();
    }
    for indel in &edge.indels {
      if indel.deletion {
        seq[indel.range.0..indel.range.1].fill(alphabet.gap());
      } else {
        seq[indel.range.0..indel.range.1].copy_from_slice(&indel.seq);
      }
    }
    seq
  } else {
    base_seq.clone()
  };

  for r in &node.seq.unknown {
    seq[r.0..r.1].fill(alphabet.unknown());
  }

  for (pos, states) in &node.profile.variable {
    seq[*pos] = alphabet.char(argmax_first(&states.dis.view()).unwrap_or(0));
  }

  seq
}

impl HasLogLh for PartitionMarginalSparse {
  fn get_log_lh(&self, node_key: GraphNodeKey) -> f64 {
    self.nodes.get(&node_key).map_or(0.0, |node| node.profile.log_lh)
  }
}

impl PartitionMarginal for PartitionMarginalSparse {}

impl PartitionRerootOps for PartitionMarginalSparse {
  fn apply_reroot(&mut self, changes: &RerootChanges) -> Result<(), Report> {
    // Handle edge split: create new node entry, move mutations to child-side edge
    if let Some(info) = &changes.edge_split {
      self
        .nodes
        .insert(info.new_node_key, SparseNodePartition::empty(&self.alphabet));

      let mut old_edge_data = self
        .edges
        .remove(&info.old_edge_key)
        .ok_or_else(|| make_internal_report!("Old edge {:?} must exist for split", info.old_edge_key))?;

      // Child-side edge gets all mutations (they describe parent->child relationship)
      old_edge_data.clear_ml_subs();
      self.edges.insert(info.child_side_edge_key, old_edge_data);

      // Parent-side edge is empty (no mutations between parent and new split node)
      self
        .edges
        .insert(info.parent_side_edge_key, SparseEdgePartition::default());
    }

    // Handle edge merge: compose mutations from two edges.
    // When the removed node is the old root (trivial root removal), update
    // root_sequence using the parent_edge subs before discarding them.
    // The parent_edge goes from the new-root side toward the removed node:
    // sub.reff() = new-root-side state, sub.qry() = removed-node state.
    if let Some(info) = &changes.edge_merge {
      let parent_edge = self
        .edges
        .remove(&info.parent_edge_key)
        .ok_or_else(|| make_internal_report!("Parent edge {:?} must exist for merge", info.parent_edge_key))?;
      let child_edge = self
        .edges
        .remove(&info.child_edge_key)
        .ok_or_else(|| make_internal_report!("Child edge {:?} must exist for merge", info.child_edge_key))?;

      if changes.inverted_edge_keys.is_empty() {
        for sub in parent_edge.fitch_subs() {
          if sub.pos() < self.root_sequence.len() {
            self.root_sequence[sub.pos()] = sub.reff();
          }
        }
        for indel in &parent_edge.indels {
          if indel.range.0 < self.root_sequence.len() && indel.range.1 <= self.root_sequence.len() {
            if indel.deletion {
              self.root_sequence[indel.range.0..indel.range.1].fill(self.alphabet.gap());
            } else {
              self.root_sequence[indel.range.0..indel.range.1].copy_from_slice(&indel.seq);
            }
          }
        }
      }

      self.nodes.remove(&info.removed_node_key);

      let merged_subs = parent_edge.chain_fitch_subs(child_edge.fitch_subs())?;
      let merged_indels = parent_edge.chain_fitch_indels(&child_edge.indels);

      let mut merged_edge = SparseEdgePartition::default();
      merged_edge.set_fitch_subs(merged_subs);
      merged_edge.indels = merged_indels;

      self.edges.insert(info.merged_edge_key, merged_edge);
    }

    // Invert edge data for each edge on the reroot path
    for edge_key in &changes.inverted_edge_keys {
      let edge_data = self
        .edges
        .get_mut(edge_key)
        .ok_or_else(|| make_internal_report!("Edge {edge_key:?} must exist on reroot path"))?;

      // Invert all substitutions
      edge_data.invert_fitch_subs();

      // Clear ML subs (stale after topology change)
      edge_data.clear_ml_subs();

      // Invert all indels
      for indel in &mut edge_data.indels {
        indel.invert();
      }

      // Swap directional messages
      mem::swap(&mut edge_data.msg_to_parent, &mut edge_data.msg_to_child);

      // Reset msg_from_child (stale propagated message cache)
      edge_data.msg_from_child = MarginalSparseSeqDistribution::default();
    }

    // Derive root sequence for the new root from the old root sequence.
    // Walk inverted edges (listed in old_root->new_root order). After sub
    // inversion, each sub's ref() holds the new-root-ward state (was the
    // original child qry). Applying ref() at each position advances the
    // sequence one step toward the new root.
    if !changes.inverted_edge_keys.is_empty() {
      let mut new_root_seq = self.root_sequence.clone();
      for edge_key in &changes.inverted_edge_keys {
        if let Some(edge_data) = self.edges.get(edge_key) {
          for sub in edge_data.fitch_subs() {
            if sub.pos() < new_root_seq.len() {
              new_root_seq[sub.pos()] = sub.reff();
            }
          }
          for indel in &edge_data.indels {
            if indel.range.0 < new_root_seq.len() && indel.range.1 <= new_root_seq.len() {
              // After inversion, deletion flag is flipped. To go in the
              // original (old_root->new_root) direction:
              // - current deletion (was insertion): child had `seq`
              // - current insertion (was deletion): child had gaps
              if indel.deletion {
                new_root_seq[indel.range.0..indel.range.1].copy_from_slice(&indel.seq);
              } else {
                new_root_seq[indel.range.0..indel.range.1].fill(self.alphabet.gap());
              }
            }
          }
        }
      }
      self.root_sequence = new_root_seq;
    }

    Ok(())
  }
}

impl PartitionBranchOps for PartitionMarginalSparse {
  fn sequence_length(&self) -> usize {
    self.length
  }

  fn edge_subs(&self, _graph: &dyn BranchTopology, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    let edge = &self.edges[&edge_key];
    match edge.ml_subs() {
      Some(subs) => Ok(subs.to_vec()),
      None => make_error!("edge_subs() called before marginal inference populated subs_ml for edge {edge_key:?}"),
    }
  }

  fn edge_effective_length(&self, graph: &dyn BranchTopology, edge_key: GraphEdgeKey) -> Result<usize, Report> {
    let (parent_key, child_key) = graph.edge_endpoints(edge_key)?;
    let parent_non_char = &self.nodes[&parent_key].seq.non_char;
    let child_non_char = &self.nodes[&child_key].seq.non_char;

    // non_char covers both gaps and unknowns (positions that do not evolve
    // under the substitution model). For internal nodes, non_char is the
    // intersection of children's non_char (Fitch backward pass), so a
    // position is excluded only when all descendants lack data there.
    let non_char_positions: usize = range_union(&[parent_non_char.clone(), child_non_char.clone()])
      .iter()
      .map(|(start, end)| end - start)
      .sum();

    Ok(self.length.saturating_sub(non_char_positions))
  }
}

impl PartitionOptimizeOps for PartitionMarginalSparse {
  fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report> {
    OptimizationContribution::from_sparse(edge_key, self)
  }

  fn edge_indel_count(&self, edge_key: GraphEdgeKey) -> usize {
    self.edges[&edge_key].indels.len()
  }
}

impl<N, E> PartitionTimetreeOps<N, E> for PartitionMarginalSparse
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  fn reconcile_topology(&mut self, graph: &Graph<N, E, ()>) {
    let graph_node_keys: BTreeSet<GraphNodeKey> = graph.get_nodes().into_iter().map(|n| n.read_arc().key()).collect();
    let graph_edge_keys: BTreeSet<GraphEdgeKey> = graph.get_edges().into_iter().map(|e| e.read_arc().key()).collect();

    // Add missing nodes with empty partition data
    for &key in &graph_node_keys {
      self
        .nodes
        .entry(key)
        .or_insert_with(|| SparseNodePartition::empty(&self.alphabet));
    }

    // Add missing edges with default partition data
    for &key in &graph_edge_keys {
      self.edges.entry(key).or_default();
    }

    // Remove stale entries for nodes/edges no longer in graph
    self.nodes.retain(|k, _| graph_node_keys.contains(k));
    self.edges.retain(|k, _| graph_edge_keys.contains(k));
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
    marginal_passes::process_node_backward(self, node)
  }

  fn process_node_forward(&mut self, graph: &Graph<N, E, ()>, node: &GraphNodeForward<N, E, ()>) -> Result<(), Report> {
    marginal_passes::process_node_forward(self, graph, node)
  }

  fn extract_ancestral_sequence(&self, node_key: GraphNodeKey) -> Seq {
    if let Some(node_data) = self.nodes.get(&node_key) {
      if !node_data.seq.sequence.is_empty() {
        node_data.seq.sequence.clone()
      } else {
        // Return empty seq to be filled later by the reconstruction algorithm
        seq![]
      }
    } else {
      seq![]
    }
  }

  fn reconstruct_node_sequence(&mut self, node: &GraphNodeForward<N, E, ()>, include_leaves: bool) -> Option<Seq> {
    if !include_leaves && node.is_leaf {
      return None;
    }

    let mut node_data = self.nodes.remove(&node.key)?;

    let (base_seq, edge) = if node.is_root {
      (&self.root_sequence, None)
    } else {
      let (parent_key, edge_key) = get_exactly_one(&node.parent_keys).ok()?;
      let parent_data = self.nodes.get(parent_key)?;
      let edge_data = self.edges.get(edge_key)?;
      (&parent_data.seq.sequence, Some(edge_data))
    };

    let seq = reconstruct_map_seq(base_seq, edge, &node_data, &self.alphabet);
    node_data.seq.sequence = seq.clone();
    self.nodes.insert(node.key, node_data);

    Some(seq)
  }

  fn get_sequence_length(&self) -> usize {
    self.length
  }
}
