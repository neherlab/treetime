use crate::commands::timetree::data::date_constraints::DateConstraint;
use crate::distribution::distribution::Distribution;
use crate::graph::graph::{GraphNodeBackward, GraphNodeForward};
use crate::graph::node::GraphNodeKey;
use crate::representation::graph_ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use crate::representation::log_lh::HasLogLh;
use crate::representation::partition_marginal::PartitionMarginalOps;
use eyre::Report;
use std::collections::BTreeMap;
use std::sync::Arc;

pub trait PartitionTreetimeMarginalOps: PartitionTimetreeOps + PartitionMarginalOps + HasLogLh {}

/// Blanket implementation: any type implementing all three traits automatically implements PartitionUnified
impl<T> PartitionTreetimeMarginalOps for T where T: PartitionTimetreeOps + PartitionMarginalOps + HasLogLh {}

pub trait PartitionTimetree {}

pub trait PartitionTimetreeOps: PartitionTimetree + Send + Sync {
  fn attach_date_constraints(
    &mut self,
    graph: &GraphAncestral,
    constraints: &BTreeMap<GraphNodeKey, DateConstraint>,
  ) -> Result<(), Report>;

  fn process_node_backward(&mut self, node: &GraphNodeBackward<NodeAncestral, EdgeAncestral, ()>)
  -> Result<(), Report>;

  fn process_node_forward(
    &mut self,
    graph: &GraphAncestral,
    node: &GraphNodeForward<NodeAncestral, EdgeAncestral, ()>,
  ) -> Result<(), Report>;

  fn extract_node_time(&self, node_key: GraphNodeKey) -> Option<f64>;

  fn extract_node_posterior(&self, node_key: GraphNodeKey) -> Option<Arc<Distribution>>;

  fn get_sequence_length(&self) -> Option<usize>;
}
