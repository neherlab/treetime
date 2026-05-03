use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
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

pub trait HasGtr {
  fn gtr(&self) -> &GTR;
  fn gtr_mut(&mut self) -> &mut GTR;
  fn sequence_length(&self) -> usize;

  fn weighted_rate(&self) -> f64 {
    self.sequence_length() as f64 * self.gtr().mu
  }

  fn normalize_rate(&mut self, scale: f64) {
    self.gtr_mut().mu /= scale;
  }
}

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

/// Derive branch mutations and effective lengths from MAP-derived partition state.
///
/// Both dense and sparse partitions implement this. Dense partitions compare
/// MAP states from full marginal posteriors. Sparse partitions return
/// pre-computed MAP-derived substitutions stored in `subs_ml` during
/// the marginal forward pass.
///
/// Requires marginal inference to have run. Pre-marginal consumers (GTR
/// inference, prune/merge) read `fitch_subs()` directly.
///
/// The graph is accepted as `&dyn BranchTopology` so the same partition trait
/// serves ancestral, timetree, and any future payloads that share the sparse
/// or dense representation.
pub trait PartitionBranchOps: Send + Sync {
  /// Return the sequence length represented by this partition.
  fn sequence_length(&self) -> usize;

  /// Return MAP-derived nucleotide substitutions for one edge.
  ///
  /// For dense partitions, computes MAP state (argmax of node posterior) at
  /// parent and child endpoints. For sparse partitions, returns pre-computed
  /// `subs_ml` (populated during forward pass). Non-canonical states
  /// and gap positions are excluded in both cases.
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
