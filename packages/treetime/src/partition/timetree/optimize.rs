use crate::partition::optimize::contribution::OptimizationContribution;
use crate::partition::timetree::partition::PartitionTimetree;
use crate::partition::traits::{PartitionOptimizeOps, PartitionRerootOps, PartitionTimetreeOps};
use eyre::Report;
use treetime_graph::edge::{EdgeOptimizeOps, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
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

impl<N, E> PartitionTimetreeOps<N, E> for PartitionTimetree
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  fn reconcile_topology(&mut self, graph: &Graph<N, E, ()>) {
    match self {
      Self::Dense(partition) => partition.reconcile_topology(graph),
      Self::Sparse(partition) => partition.reconcile_topology(graph),
    }
  }
}
