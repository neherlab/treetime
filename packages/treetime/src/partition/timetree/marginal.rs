use crate::ancestral::sample::SampleMode;
use crate::partition::timetree::partition::PartitionTimetree;
use eyre::Report;
use rayon::prelude::*;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::GraphNodeForward;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::FastaRecord;
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
    aln: &[FastaRecord],
    names: &BTreeMap<GraphNodeKey, Option<String>>,
  ) -> Result<(), Report> {
    if let Self::Dense(family) = self {
      family.node_states = family.partition.attach_sequences(graph, aln, names)?;
    }
    Ok(())
  }

  /// Run a full marginal update in place, replacing the result maps and returning the substitution
  /// log likelihood.
  pub fn run_marginal_update(
    &mut self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  ) -> Result<LogLh, Report> {
    match self {
      Self::Dense(family) => family.run_marginal_update(graph, branch_lengths),
      Self::Sparse(family) => family.run_marginal_update(graph, branch_lengths),
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
        &family.forward,
        node,
        include_leaves,
        impute,
        sample_mode,
        rng,
      ),
    }
  }
}

/// Sum of per-partition root log-likelihoods after marginal reconstruction.
pub fn graph_log_lh(graph: &Graph, partitions: &[PartitionTimetree]) -> Result<LogLh, Report> {
  let root_key = graph.get_exactly_one_root()?.read_arc().key();
  let log_lh = partitions
    .par_iter()
    .map(|partition| partition.get_log_lh(root_key))
    .collect::<Vec<_>>()
    .into_iter()
    .sum();
  Ok(log_lh)
}

/// Run a marginal update over each timetree partition in place and return the summed substitution log
/// likelihood.
pub fn marginal_update_timetree(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  partitions: &mut [PartitionTimetree],
) -> Result<LogLh, Report> {
  let mut total = LogLh::ZERO;
  for partition in partitions.iter_mut() {
    total += partition.run_marginal_update(graph, branch_lengths)?;
  }
  Ok(total)
}

/// Attach leaf sequences (dense) then run the initial marginal update over every timetree partition.
pub fn initialize_marginal_timetree(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  partitions: &mut [PartitionTimetree],
  aln: &[FastaRecord],
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<LogLh, Report> {
  for partition in partitions.iter_mut() {
    partition.attach_sequences(graph, aln, names)?;
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
