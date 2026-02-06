use crate::alphabet::alphabet::Alphabet;
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

  fn supports_reroot_edge_split(&self) -> bool {
    false
  }

  fn supports_old_root_removal(&self) -> bool {
    false
  }

  fn reroot_partition_node_only(
    &mut self,
    _graph: &Graph<N, E, ()>,
    old_root_key: GraphNodeKey,
    new_root_key: GraphNodeKey,
    path_from_old_to_new: &[(GraphNodeKey, Option<GraphEdgeKey>)],
  ) -> Result<(), Report> {
    // Build ordered list of edge keys on the reroot path (old root -> new root direction)
    let reroot_edge_keys = path_from_old_to_new
      .iter()
      .filter_map(|(_, edge_key)| *edge_key)
      .collect_vec();

    // Compute new root sequence by applying edge mutations from old root toward new root
    // Path is already in old->new direction, so iterate forward
    let mut new_root_seq = self.nodes[&old_root_key].seq.sequence.clone();
    for edge_key in &reroot_edge_keys {
      let edge_data = &self.edges[edge_key];
      // Apply substitutions
      for sub in &edge_data.subs {
        // Original direction: parent->child means reff->qry
        // We're going child->parent, so we apply reff (the parent state)
        new_root_seq[sub.pos()] = sub.reff();
      }
      // Apply indels (reverse direction)
      for indel in &edge_data.indels {
        if indel.deletion {
          // Original: deletion means child has gap. Reverse: child has sequence
          new_root_seq[indel.range.0..indel.range.1].copy_from_slice(&indel.seq);
        } else {
          // Original: insertion means child has sequence. Reverse: child has gap
          new_root_seq[indel.range.0..indel.range.1].fill(self.alphabet.gap());
        }
      }
    }

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

    // Set new root sequence
    self
      .nodes
      .get_mut(&new_root_key)
      .ok_or_else(|| make_internal_report!("New root node {new_root_key:?} must exist"))?
      .seq
      .sequence = new_root_seq;

    // Clear old root sequence
    self
      .nodes
      .get_mut(&old_root_key)
      .ok_or_else(|| make_internal_report!("Old root node {old_root_key:?} must exist"))?
      .seq
      .sequence = crate::seq![];

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
