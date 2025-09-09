use crate::graph::graph::{GraphNodeBackward, GraphNodeForward};
use crate::graph::node::GraphNodeKey;
use crate::io::fasta::FastaRecord;
use crate::representation::graph_ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use crate::representation::seq::Seq;
use eyre::Report;

pub trait PartitionMarginal {}

/// Trait defining the core operations needed for marginal reconstruction
/// These methods encapsulate the differences between sparse and dense implementations
pub trait PartitionMarginalOps: PartitionMarginal + Send + Sync {
  /// Initialize partition with sequence data (for dense implementations)
  fn attach_sequences(&mut self, graph: &GraphAncestral, aln: &[FastaRecord]) -> Result<(), Report>;

  /// Process a node during backward pass (equivalent to run_marginal_*_backward)
  fn process_node_backward(&mut self, node: &GraphNodeBackward<NodeAncestral, EdgeAncestral, ()>)
  -> Result<(), Report>;

  /// Process a node during forward pass (equivalent to run_marginal_*_forward)
  fn process_node_forward(
    &mut self,
    graph: &GraphAncestral,
    node: &GraphNodeForward<NodeAncestral, EdgeAncestral, ()>,
  ) -> Result<(), Report>;

  /// Extract ancestral sequence from node profile
  fn extract_ancestral_sequence(&mut self, node_key: GraphNodeKey) -> Seq;

  /// Reconstruct ancestral sequences for a specific node during tree traversal
  fn reconstruct_node_sequence(
    &mut self,
    node: &GraphNodeForward<NodeAncestral, EdgeAncestral, ()>,
    include_leaves: bool,
  ) -> Option<Seq>;

  /// Get sequence length for branch length normalization
  fn get_length(&self) -> usize;
}
