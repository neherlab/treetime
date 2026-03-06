use crate::alphabet::alphabet::Alphabet;
use crate::representation::payload::sparse::{SparseEdgePartition, SparseNodePartition};
use eyre::Report;
use parking_lot::RwLock;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_graph::edge::{EdgeOptimizeOps, GraphEdge, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::{GraphNodeBackward, GraphNodeForward};
use treetime_graph::node::{GraphNode, GraphNodeKey, Named};
use treetime_io::fasta::FastaRecord;
use treetime_primitives::Seq;

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
  fn extract_ancestral_sequence(&self, node_key: GraphNodeKey) -> Seq;

  /// Reconstruct ancestral sequences for a specific node during tree traversal
  fn reconstruct_node_sequence(&mut self, node: &GraphNodeForward<N, E, ()>, include_leaves: bool) -> Option<Seq>;

  /// Get sequence length
  fn get_sequence_length(&self) -> usize;
}

pub trait PartitionCompressed: Sync + Send {
  fn index(&self) -> usize;

  fn alphabet(&self) -> &Alphabet;

  fn length(&self) -> usize;

  fn nodes(&self) -> &BTreeMap<GraphNodeKey, SparseNodePartition>;

  fn edges(&self) -> &BTreeMap<GraphEdgeKey, SparseEdgePartition>;

  fn nodes_mut(&mut self) -> &mut BTreeMap<GraphNodeKey, SparseNodePartition>;

  fn edges_mut(&mut self) -> &mut BTreeMap<GraphEdgeKey, SparseEdgePartition>;

  fn node(&self, key: &GraphNodeKey) -> &SparseNodePartition {
    self.nodes().get(key).expect("Node not found")
  }

  fn edge(&self, key: &GraphEdgeKey) -> &SparseEdgePartition {
    self.edges().get(key).expect("Edge not found")
  }

  fn node_mut(&mut self, key: &GraphNodeKey) -> &mut SparseNodePartition {
    self.nodes_mut().get_mut(key).expect("Node not found")
  }

  fn edge_mut(&mut self, key: &GraphEdgeKey) -> &mut SparseEdgePartition {
    self.edges_mut().get_mut(key).expect("Edge not found")
  }
}

pub trait HasLogLh {
  /// Get the log likelihood contribution for a given node
  fn get_log_lh(&self, node_key: GraphNodeKey) -> f64;
}

/// Calculate the total log likelihood of the graph given the partitions
pub fn graph_log_lh<P, N, E, D>(graph: &Graph<N, E, D>, partitions: &[Arc<RwLock<P>>]) -> Result<f64, Report>
where
  P: HasLogLh + ?Sized,
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send + Default,
{
  let root = graph.get_exactly_one_root()?;
  let root_key = root.read_arc().key();

  let log_lh = partitions
    .iter()
    .map(|partition| partition.read_arc().get_log_lh(root_key))
    .sum();

  Ok(log_lh)
}
