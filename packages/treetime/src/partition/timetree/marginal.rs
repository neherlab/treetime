use crate::ancestral::sample::SampleMode;
use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::timetree::partition::{GraphTimetree, PartitionTimetree};
use crate::partition::traits::{HasLogLh, PartitionMarginalOps, PartitionMarginalPasses};
use crate::payload::timetree::{EdgeTimetree, NodeTimetree};
use eyre::Report;
use treetime_graph::graph_traverse::GraphNodeForward;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::FastaRecord;
use treetime_primitives::{LogLh, Seq};

impl HasLogLh for PartitionTimetree {
  fn get_log_lh(&self, node_key: GraphNodeKey) -> LogLh {
    match self {
      Self::Dense(partition) => partition.get_log_lh(node_key),
      Self::Sparse(partition) => partition.get_log_lh(node_key),
    }
  }

  fn reset_node_log_likelihoods(&mut self) {
    match self {
      Self::Dense(partition) => partition.reset_node_log_likelihoods(),
      Self::Sparse(partition) => partition.reset_node_log_likelihoods(),
    }
  }
}

impl PartitionMarginalPasses<NodeTimetree, EdgeTimetree> for PartitionTimetree {
  fn process_backward_pass(&mut self, graph: &GraphTimetree) -> Result<(), Report> {
    match self {
      Self::Dense(partition) => partition.process_backward_pass(graph),
      Self::Sparse(partition) => partition.process_backward_pass(graph),
    }
  }

  fn process_forward_pass(&mut self, graph: &GraphTimetree) -> Result<(), Report> {
    match self {
      Self::Dense(partition) => partition.process_forward_pass(graph),
      Self::Sparse(partition) => partition.process_forward_pass(graph),
    }
  }

  fn get_sequence_length(&self) -> usize {
    match self {
      Self::Dense(partition) => partition.get_sequence_length(),
      Self::Sparse(partition) => partition.get_sequence_length(),
    }
  }
}

impl PartitionMarginalOps<NodeTimetree, EdgeTimetree> for PartitionTimetree {
  fn attach_sequences(&mut self, graph: &GraphTimetree, aln: &[FastaRecord]) -> Result<(), Report> {
    match self {
      Self::Dense(partition) => partition.attach_sequences(graph, aln),
      Self::Sparse(partition) => partition.attach_sequences(graph, aln),
    }
  }

  fn extract_ancestral_sequence(&self, node_key: GraphNodeKey) -> Seq {
    match self {
      Self::Dense(partition) => {
        <PartitionMarginalDense as PartitionMarginalOps<NodeTimetree, EdgeTimetree>>::extract_ancestral_sequence(
          partition, node_key,
        )
      },
      Self::Sparse(partition) => {
        <PartitionMarginalSparse as PartitionMarginalOps<NodeTimetree, EdgeTimetree>>::extract_ancestral_sequence(
          partition, node_key,
        )
      },
    }
  }

  fn reconstruct_node_sequence(
    &mut self,
    node: &GraphNodeForward<NodeTimetree, EdgeTimetree>,
    include_leaves: bool,
    impute: bool,
    sample_mode: SampleMode,
    rng: &mut dyn rand::RngCore,
  ) -> Option<Seq> {
    match self {
      Self::Dense(partition) => partition.reconstruct_node_sequence(node, include_leaves, impute, sample_mode, rng),
      Self::Sparse(partition) => partition.reconstruct_node_sequence(node, include_leaves, impute, sample_mode, rng),
    }
  }
}
