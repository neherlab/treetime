use treetime_graph::edge::EdgeOptimizeOps;
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::{GraphNodeBackward, GraphNodeForward};
use treetime_graph::node::{GraphNode, GraphNodeKey, Named};
use crate::io::fasta::FastaRecord;
use crate::representation::seq::Seq;
use eyre::Report;

pub trait PartitionMarginal {}

/// Trait defining the core operations needed for marginal reconstruction
/// These methods encapsulate the differences between sparse and dense implementations
pub trait PartitionMarginalOps<N, E>: PartitionMarginal + Send + Sync
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  /// Initialize partition with sequence data (for dense implementations)
  fn attach_sequences(&mut self, graph: &Graph<N, E, ()>, aln: &[FastaRecord]) -> Result<(), Report>;

  /// Process a node during backward pass (equivalent to run_marginal_*_backward)
  fn process_node_backward(&mut self, node: &GraphNodeBackward<N, E, ()>) -> Result<(), Report>;

  /// Process a node during forward pass (equivalent to run_marginal_*_forward)
  fn process_node_forward(&mut self, graph: &Graph<N, E, ()>, node: &GraphNodeForward<N, E, ()>) -> Result<(), Report>;

  /// Extract ancestral sequence from node profile
  fn extract_ancestral_sequence(&mut self, node_key: GraphNodeKey) -> Seq;

  /// Reconstruct ancestral sequences for a specific node during tree traversal
  fn reconstruct_node_sequence(&mut self, node: &GraphNodeForward<N, E, ()>, include_leaves: bool) -> Option<Seq>;

  /// Get sequence length
  fn get_sequence_length(&self) -> usize;
}
