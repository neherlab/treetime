use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::sample::SampleMode;
use crate::gtr::gtr::GTR;
use crate::make_error;
use crate::partition::marginal::sparse::pass;
use crate::partition::marginal::sparse::reconstruct::{reconstruct_leaf_sequence, reconstruct_map_seq_sampled};
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::partition::storage::sparse::{SparseEdgePartition, SparseNodePartition};
use crate::partition::traits::{
  BranchTopology, HasGtr, HasLogLh, PartitionBranchOps, PartitionMarginalOps, PartitionMarginalPasses,
  PartitionOptimizeOps, PartitionTimetreeOps,
};
use crate::seq::composition::Composition;
use crate::seq::mutation::Sub;
use eyre::Report;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::{EdgeOptimizeOps, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::GraphNodeForward;
use treetime_graph::node::{GraphNode, GraphNodeKey, Named};
use treetime_io::fasta::FastaRecord;
use treetime_primitives::{AlphabetLike, LogLh, Seq, seq};
use treetime_utils::collections::container::get_exactly_one;
use treetime_utils::interval::range_union::range_union;

#[derive(Clone, Debug, Serialize)]
pub struct PartitionMarginalSparse {
  pub index: usize,
  pub gtr: GTR,
  pub alphabet: Alphabet,
  pub length: usize,
  pub root_sequence: Seq,
  pub nodes: BTreeMap<GraphNodeKey, SparseNodePartition>,
  pub edges: BTreeMap<GraphEdgeKey, SparseEdgePartition>,
}

impl PartitionMarginalSparse {
  #[allow(clippy::same_name_method)]
  pub fn get_sequence_length(&self) -> usize {
    self.length
  }
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

impl HasLogLh for PartitionMarginalSparse {
  fn get_log_lh(&self, node_key: GraphNodeKey) -> LogLh {
    self
      .nodes
      .get(&node_key)
      .map_or(LogLh::ZERO, |node| node.profile.log_lh)
  }

  fn reset_node_log_likelihoods(&mut self) {
    for node_data in self.nodes.values_mut() {
      node_data.profile.log_lh = LogLh::ZERO;
    }
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

  fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<crate::seq::indel::InDel> {
    self.edges[&edge_key].indels.clone()
  }

  fn root_sequence(&self, _graph: &dyn BranchTopology) -> Result<Seq, Report> {
    Ok(self.root_sequence.clone())
  }

  fn node_sequence(&self, node_key: GraphNodeKey) -> Seq {
    self.nodes[&node_key].seq.sequence.clone()
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

impl<N, E> PartitionMarginalPasses<N, E> for PartitionMarginalSparse
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  fn process_backward_pass(&mut self, graph: &Graph<N, E, ()>) -> Result<(), Report> {
    pass::process_backward_indexed(self, graph)
  }

  fn process_forward_pass(&mut self, graph: &Graph<N, E, ()>) -> Result<(), Report> {
    pass::process_forward_indexed(self, graph)
  }

  fn get_sequence_length(&self) -> usize {
    self.length
  }
}

impl<N, E> PartitionMarginalOps<N, E> for PartitionMarginalSparse
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  fn attach_sequences(&mut self, _graph: &Graph<N, E, ()>, _aln: &[FastaRecord]) -> Result<(), Report> {
    Ok(())
  }

  fn extract_ancestral_sequence(&self, node_key: GraphNodeKey) -> Seq {
    if let Some(node_data) = self.nodes.get(&node_key) {
      if !node_data.seq.sequence.is_empty() {
        node_data.seq.sequence.clone()
      } else {
        seq![]
      }
    } else {
      seq![]
    }
  }

  fn reconstruct_node_sequence(
    &mut self,
    node: &GraphNodeForward<N, E>,
    include_leaves: bool,
    impute: bool,
    sample_mode: SampleMode,
    rng: &mut dyn rand::RngCore,
  ) -> Option<Seq> {
    let mut node_data = self.nodes.remove(&node.key)?;

    let (base_seq, edge) = if node.is_root {
      (&self.root_sequence, None)
    } else {
      let (parent_key, edge_key) = get_exactly_one(&node.parent_keys).ok()?;
      let parent_data = self.nodes.get(parent_key)?;
      let edge_data = self.edges.get(edge_key)?;
      (&parent_data.seq.sequence, Some(edge_data))
    };

    // Tips reconstruct from their own observed state; the internal-node and root path is left
    // byte-identical (it feeds optimize/timetree). The leaf branch fixes the tip corruption and,
    // when requested, imputes missing states.
    let seq = if node.is_leaf {
      reconstruct_leaf_sequence(&node_data, edge, base_seq, impute, &self.alphabet)
    } else {
      let sample = sample_mode.samples_node(node.is_root);
      reconstruct_map_seq_sampled(base_seq, edge, &node_data, &self.alphabet, sample, rng)
    };

    node_data.seq.composition = Composition::with_seq(&seq, self.alphabet.chars(), self.alphabet.gap());
    node_data.seq.sequence = seq.clone();
    self.nodes.insert(node.key, node_data);

    // A suppressed tip is still reconstructed above (so the node-data serializer reads the corrected
    // sequence), but is not emitted to the reconstructed-FASTA visitor.
    if !include_leaves && node.is_leaf {
      return None;
    }

    Some(seq)
  }
}
