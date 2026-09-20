use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
use crate::partition::marginal::dense::reroot::reroot_dense;
use crate::partition::marginal::shared::reconcile::{live_node_keys, reconcile_node_states};
use crate::partition::marginal::sparse::reroot::reroot_sparse;
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::partition::storage::dense::DenseNodeState;
use crate::partition::storage::sparse::SparseNodeState;
use crate::partition::timetree::partition::PartitionTimetree;
use eyre::Report;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::reroot::RerootChanges;

#[allow(clippy::multiple_inherent_impl, reason = "PartitionTimetree's inherent impl is split by concern: the type definition in partition.rs, the timetree optimization methods here")]
impl PartitionTimetree {
  /// The sequence length this partition represents.
  pub fn sequence_length(&self) -> usize {
    match self {
      Self::Dense(family) => family.sequence_length(),
      Self::Sparse(family) => family.sequence_length(),
    }
  }

  /// The per-edge branch-length optimization contribution, built from the last update's messages.
  pub fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report> {
    match self {
      Self::Dense(family) => Ok(family.create_edge_contribution(edge_key)),
      Self::Sparse(family) => family.create_edge_contribution(edge_key),
    }
  }

  /// The number of indel events on one edge.
  pub fn edge_indel_count(&self, edge_key: GraphEdgeKey) -> usize {
    match self {
      Self::Dense(family) => family.edge_indel_count(edge_key),
      Self::Sparse(family) => family.edge_indel_count(edge_key),
    }
  }

  /// Apply a reroot at the structural operation, returning the partition over the rerooted topology.
  ///
  /// The partition is consumed: the reroot rewrites the durable observations and carries the node
  /// states across, and the per-edge results of the previous update do not survive it.
  pub fn apply_reroot(self, changes: &RerootChanges) -> Result<Self, Report> {
    Ok(match self {
      Self::Dense(family) => Self::Dense(reroot_dense(family.partition, family.gtr, family.node_states, changes)),
      Self::Sparse(family) => Self::Sparse(reroot_sparse(
        family.partition,
        family.gtr,
        family.node_states,
        changes,
      )?),
    })
  }

  /// Reconcile the partition to the graph after a topology change (polytomy resolution): give the
  /// observations and node states an entry for every current node, drop entries for nodes that are
  /// gone, and carry no per-edge results across. The next marginal update recomputes the values.
  #[must_use]
  pub fn reconcile_topology(self, graph: &Graph) -> Self {
    let live_nodes = live_node_keys(graph);
    match self {
      Self::Dense(family) => Self::Dense(DenseReconstruction::seeded(
        family.partition,
        family.gtr,
        reconcile_node_states(family.node_states, &live_nodes, DenseNodeState::empty),
      )),
      Self::Sparse(family) => {
        let mut partition = family.partition;
        partition.reconcile_topology(graph);
        Self::Sparse(SparseReconstruction::seeded(
          partition,
          family.gtr,
          reconcile_node_states(family.node_states, &live_nodes, SparseNodeState::empty),
        ))
      },
    }
  }
}
