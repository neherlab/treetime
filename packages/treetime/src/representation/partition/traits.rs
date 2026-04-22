use crate::alphabet::alphabet::Alphabet;
use crate::make_internal_error;
use crate::make_internal_report;
use crate::representation::payload::sparse::{SparseEdgePartition, SparseNodePartition};
use crate::seq::mutation::Sub;
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

/// Minimal graph-structure abstraction used by per-branch partition operations.
///
/// Exists so that `PartitionBranchOps` can serve any phylogenetic graph payload
/// (ancestral, timetree, etc.) without binding to a concrete `Graph<N, E, D>`.
/// Only the operations actually needed by branch-level computations are exposed:
/// resolving an edge to its endpoints and walking one step toward the root.
/// Trait-object safe so `&dyn BranchTopology` can flow through dynamic dispatch
/// while callers keep their concrete `&Graph<N, E, D>` thanks to unsized
/// coercion.
pub trait BranchTopology: Send + Sync {
  /// Return `(parent_node_key, child_node_key)` for one edge.
  fn edge_endpoints(&self, edge_key: GraphEdgeKey) -> Result<(GraphNodeKey, GraphNodeKey), Report>;

  /// Return `Some((parent_node_key, parent_edge_key))` for a non-root node,
  /// or `None` when the node is the root. Errors when the node has more than
  /// one parent (the algorithm only supports trees).
  fn node_parent(&self, node_key: GraphNodeKey) -> Result<Option<(GraphNodeKey, GraphEdgeKey)>, Report>;
}

impl<N, E, D> BranchTopology for Graph<N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: Send + Sync,
{
  fn edge_endpoints(&self, edge_key: GraphEdgeKey) -> Result<(GraphNodeKey, GraphNodeKey), Report> {
    let edge = self
      .get_edge(edge_key)
      .ok_or_else(|| make_internal_report!("Edge {edge_key} not found"))?;
    let edge = edge.read_arc();
    Ok((edge.source(), edge.target()))
  }

  fn node_parent(&self, node_key: GraphNodeKey) -> Result<Option<(GraphNodeKey, GraphEdgeKey)>, Report> {
    let node = self
      .get_node(node_key)
      .ok_or_else(|| make_internal_report!("Node {node_key} not found"))?;
    let node = node.read_arc();
    let inbound = node.inbound();
    match inbound.len() {
      0 => Ok(None),
      1 => {
        let parent_edge_key = inbound[0];
        let parent_node_key = self.get_source_node_key(parent_edge_key)?;
        Ok(Some((parent_node_key, parent_edge_key)))
      },
      n => make_internal_error!("Node {node_key} has {n} parents; only trees are supported"),
    }
  }
}

/// Derive branch mutations and effective lengths from the current partition state.
///
/// Both dense and sparse partitions implement this: dense partitions compare
/// MAP states from full marginal posteriors, sparse partitions compare
/// reconstructed states from Fitch data and current variable-site assignments.
///
/// Post-marginal consumers (GTR inference, prune, output) should use this
/// trait instead of reading `SparseEdgePartition.subs` directly, which
/// reflects Fitch parsimony data and becomes stale after marginal inference.
///
/// The graph is accepted as `&dyn BranchTopology` so the same partition trait
/// serves ancestral, timetree, and any future payloads that share the sparse
/// or dense representation.
pub trait PartitionBranchOps: Send + Sync {
  /// Return the sequence length represented by this partition.
  fn sequence_length(&self) -> usize;

  /// Return nucleotide substitutions for one edge derived from current state.
  ///
  /// For sparse partitions, reconstructs parent and child states from Fitch
  /// parsimony data and current variable-site assignments, then compares.
  /// For dense partitions, takes the MAP state (argmax of the node posterior)
  /// at the parent and child endpoints. Non-canonical states and gap positions
  /// are excluded in both cases.
  fn edge_subs(&self, graph: &dyn BranchTopology, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report>;

  /// Return the number of alignment positions where both parent and child
  /// have canonical (non-gap, non-ambiguous) states for one edge.
  fn edge_effective_length(&self, graph: &dyn BranchTopology, edge_key: GraphEdgeKey) -> Result<usize, Report>;
}

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

  /// Post-Fitch finalization hook called at the end of `compress_sequences`.
  ///
  /// Default is a no-op. `PartitionMarginalSparse` overrides this to extract
  /// the root sequence off-tree and clear internal node sequences.
  fn finalize_fitch<N, E>(&mut self, _graph: &Graph<N, E, ()>) -> Result<(), Report>
  where
    N: GraphNode,
    E: GraphEdge,
  {
    Ok(())
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
