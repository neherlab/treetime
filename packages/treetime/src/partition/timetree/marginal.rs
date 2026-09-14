use crate::ancestral::sample::SampleMode;
use crate::partition::marginal::shared::update::MarginalPasses;
use crate::partition::timetree::partition::PartitionTimetree;
use crate::seq::indel::InDel;
use crate::seq::mutation::{Mutation, MutationTrack, Sub, combine_edge_mutations};
use eyre::Report;
use rayon::prelude::*;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::GraphNodeForward;
use treetime_graph::node::GraphNodeKey;
use treetime_io::nwk::NwkFastaNodeInput;
use treetime_primitives::{LogLh, Seq, seq};

impl PartitionTimetree {
  /// Sequence length represented by this partition.
  pub fn get_sequence_length(&self) -> usize {
    match self {
      Self::Dense(family) => family.partition.length,
      Self::Sparse(family) => family.partition.length,
    }
  }

  /// Root-relative substitution log likelihood for one node, read from this partition's node states.
  pub fn get_log_lh(&self, node_key: GraphNodeKey) -> LogLh {
    match self {
      Self::Dense(family) => family.partition.get_log_lh(&family.node_states, node_key),
      Self::Sparse(family) => family.partition.get_log_lh(&family.node_states, node_key),
    }
  }

  /// Seed leaf node states for the dense representation from the alignment; a no-op for sparse, whose
  /// node states are seeded from the Fitch handoff at construction.
  pub fn attach_sequences(
    &mut self,
    graph: &Graph,
    node_inputs: &BTreeMap<GraphNodeKey, NwkFastaNodeInput>,
  ) -> Result<(), Report> {
    if let Self::Dense(family) = self {
      family.node_states = family.partition.attach_sequences(graph, node_inputs)?;
    }
    Ok(())
  }

  /// Run a full marginal update, returning the partition at the refreshed node states, messages, and
  /// estimates together with the substitution log likelihood.
  pub fn marginal_update(
    self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  ) -> Result<(Self, LogLh), Report> {
    match self {
      Self::Dense(family) => {
        let (family, log_lh) = family.marginal_update(graph, branch_lengths)?;
        Ok((Self::Dense(family), log_lh))
      },
      Self::Sparse(family) => {
        let (family, log_lh) = family.marginal_update(graph, branch_lengths)?;
        Ok((Self::Sparse(family), log_lh))
      },
    }
  }

  /// The deterministic most-likely-state sequence for one node (convergence reads).
  pub fn extract_ancestral_sequence(&self, node_key: GraphNodeKey) -> Seq {
    match self {
      Self::Dense(family) => family
        .partition
        .extract_ancestral_sequence(&family.node_states, node_key),
      Self::Sparse(family) => family
        .partition
        .extract_ancestral_sequence(&family.node_states, node_key),
    }
  }

  /// Reconstruct one node's output sequence, recording it into the node state, or `None` for a
  /// suppressed tip.
  pub fn reconstruct_node_sequence(
    &mut self,
    node: &GraphNodeForward,
    include_leaves: bool,
    impute: bool,
    sample_mode: SampleMode,
    rng: &mut dyn rand::RngCore,
  ) -> Option<Seq> {
    match self {
      Self::Dense(family) => family.partition.reconstruct_node_sequence(
        &mut family.node_states,
        node,
        include_leaves,
        impute,
        sample_mode,
        rng,
      ),
      Self::Sparse(family) => family.partition.reconstruct_node_sequence(
        &mut family.node_states,
        &family.edges.forward,
        node,
        include_leaves,
        impute,
        sample_mode,
        rng,
      ),
    }
  }

  /// MAP-derived nucleotide substitutions on one edge (parent -> child).
  pub fn edge_subs(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    match self {
      Self::Dense(family) => family.edge_subs(graph, edge_key),
      Self::Sparse(family) => family.edge_subs(edge_key),
    }
  }

  /// Grouped aligned insertions and deletions for one edge.
  pub fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<InDel> {
    match self {
      Self::Dense(family) => family.edge_indels(edge_key),
      Self::Sparse(family) => family.edge_indels(edge_key),
    }
  }

  /// The reconstructed root sequence.
  pub fn root_sequence(&self, graph: &Graph) -> Result<Seq, Report> {
    match self {
      Self::Dense(family) => family.root_sequence(graph),
      Self::Sparse(family) => family.root_sequence(graph),
    }
  }

  /// The reconstructed sequence for one node.
  pub fn node_sequence(&self, node_key: GraphNodeKey) -> Seq {
    match self {
      Self::Dense(family) => family.node_sequence(node_key),
      Self::Sparse(family) => family.node_sequence(node_key),
    }
  }

  /// The substitutions and indels on one edge as one mutation list on the given track.
  pub fn edge_mutations(
    &self,
    graph: &Graph,
    edge_key: GraphEdgeKey,
    track: &MutationTrack,
  ) -> Result<Vec<Mutation>, Report> {
    combine_edge_mutations(self.edge_subs(graph, edge_key)?, &self.edge_indels(edge_key), track)
  }
}

/// Sum of per-partition root log-likelihoods after marginal reconstruction.
pub fn graph_log_lh(graph: &Graph, partitions: &[PartitionTimetree]) -> Result<LogLh, Report> {
  let root_key = graph.get_exactly_one_root()?.key();
  let log_lh = partitions
    .par_iter()
    .map(|partition| partition.get_log_lh(root_key))
    .collect::<Vec<_>>()
    .into_iter()
    .sum();
  Ok(log_lh)
}

/// Run a marginal update over every timetree partition, returning the updated partitions and the
/// summed substitution log likelihood.
pub fn marginal_update_timetree(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  partitions: Vec<PartitionTimetree>,
) -> Result<(Vec<PartitionTimetree>, LogLh), Report> {
  partitions
    .into_iter()
    .try_fold((Vec::new(), LogLh::ZERO), |(mut updated, total), partition| {
      let (partition, log_lh) = partition.marginal_update(graph, branch_lengths)?;
      updated.push(partition);
      Ok((updated, total + log_lh))
    })
}

/// Attach leaf sequences (dense) then run the initial marginal update over every timetree partition.
pub fn initialize_marginal_timetree(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  mut partitions: Vec<PartitionTimetree>,
  node_inputs: &BTreeMap<GraphNodeKey, NwkFastaNodeInput>,
) -> Result<(Vec<PartitionTimetree>, LogLh), Report> {
  for partition in &mut partitions {
    partition.attach_sequences(graph, node_inputs)?;
  }
  marginal_update_timetree(graph, branch_lengths, partitions)
}

/// Walk the graph in preorder, reconstructing every node's sequence from the first partition,
/// emitting each reconstructed sequence to `visitor`, and returning the reconstructed sequences keyed
/// by node id. Mirrors the ancestral command's reconstruction walk for the timetree partitions.
pub fn ancestral_reconstruction_timetree(
  graph: &Graph,
  include_leaves: bool,
  impute: bool,
  partitions: &mut [PartitionTimetree],
  sample_mode: SampleMode,
  rng: &mut dyn rand::RngCore,
  mut visitor: impl FnMut(GraphNodeKey, &Seq) -> Result<(), Report>,
) -> Result<BTreeMap<GraphNodeKey, Seq>, Report> {
  let mut node_sequences = BTreeMap::new();
  graph.iter_depth_first_preorder_forward(|node| {
    if partitions.is_empty() {
      if !include_leaves && node.is_leaf {
        return Ok(());
      }
      let seq = seq![];
      visitor(node.key, &seq)?;
      node_sequences.insert(node.key, seq);
      return Ok(());
    }

    let reconstructed = partitions[0].reconstruct_node_sequence(&node, include_leaves, impute, sample_mode, rng);
    match reconstructed {
      Some(seq) => {
        visitor(node.key, &seq)?;
        node_sequences.insert(node.key, seq);
        Ok(())
      },
      None => Ok(()),
    }
  })?;
  Ok(node_sequences)
}
