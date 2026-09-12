use crate::partition::optimize::contribution::OptimizationContribution;
use crate::partition::timetree::partition::PartitionTimetree;
use crate::partition::traits::{PartitionOptimizeOps, PartitionRerootOps};
use eyre::Report;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::reroot::RerootChanges;

impl PartitionOptimizeOps for PartitionTimetree {
  fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report> {
    match self {
      Self::Dense(partition) => partition.create_edge_contribution(edge_key),
      Self::Sparse(partition) => partition.create_edge_contribution(edge_key),
    }
  }

  fn edge_indel_count(&self, edge_key: GraphEdgeKey) -> usize {
    match self {
      Self::Dense(partition) => partition.edge_indel_count(edge_key),
      Self::Sparse(partition) => partition.edge_indel_count(edge_key),
    }
  }
}

impl PartitionRerootOps for PartitionTimetree {
  fn apply_reroot(&mut self, changes: &RerootChanges) -> Result<(), Report> {
    match self {
      Self::Dense(partition) => partition.apply_reroot(changes),
      Self::Sparse(partition) => partition.apply_reroot(changes),
    }
  }
}

impl PartitionTimetree {
  /// Ensure the partition has entries for all nodes and edges in the graph.
  ///
  /// After topology changes (polytomy resolution), new nodes/edges may lack partition entries.
  /// This adds empty/default entries for missing elements and removes stale entries for
  /// elements no longer in the graph. The subsequent marginal update pass recomputes values.
  pub fn reconcile_topology(&mut self, graph: &Graph<()>) {
    match self {
      Self::Dense(partition) => partition.reconcile_topology(graph),
      Self::Sparse(partition) => partition.reconcile_topology(graph),
    }
  }
}
