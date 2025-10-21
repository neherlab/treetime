use crate::graph::edge::{GraphEdge, Weighted};
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, Named};
use crate::representation::edge_timetree::EdgeTimetree;
use crate::representation::log_lh::HasLogLh;
use crate::representation::node_timetree::NodeTimetree;
use crate::representation::partition_marginal::PartitionMarginalOps;
use eyre::Report;
use parking_lot::RwLock;
use std::sync::Arc;

pub type GraphTimetree = Graph<NodeTimetree, EdgeTimetree, ()>;
pub type PartitionTreetimeMarginalVec = Vec<Arc<RwLock<dyn PartitionTreetimeMarginalOps<NodeTimetree, EdgeTimetree>>>>;

pub trait PartitionTreetimeMarginalOps<N, E>:
  PartitionTimetreeOps<N, E> + PartitionMarginalOps<N, E> + HasLogLh
where
  N: GraphNode + Named,
  E: GraphEdge + Weighted,
{
}

/// Blanket implementation: any type implementing all three traits automatically implements PartitionTreetimeMarginalOps
impl<T, N, E> PartitionTreetimeMarginalOps<N, E> for T
where
  T: PartitionTimetreeOps<N, E> + PartitionMarginalOps<N, E> + HasLogLh,
  N: GraphNode + Named,
  E: GraphEdge + Weighted,
{
}

pub trait PartitionTimetree {}

pub trait PartitionTimetreeOps<N, E>: PartitionTimetree + Send + Sync
where
  N: GraphNode + Named,
  E: GraphEdge + Weighted,
{
  fn initialize_nodes(&mut self, graph: &Graph<N, E, ()>) -> Result<(), Report>;

  fn get_sequence_length(&self) -> Option<usize>;
}
