use crate::distribution::distribution::Distribution;
use crate::graph::edge::GraphEdge;
use crate::graph::node::GraphNode;
use std::sync::Arc;

/// Trait for node types that support timetree inference.
/// Provides access to time distribution and estimated time fields.
pub trait TimetreeNode: GraphNode {
  fn time_distribution(&self) -> &Option<Arc<Distribution>>;
  fn set_time_distribution(&mut self, dist: Option<Arc<Distribution>>);
  fn time(&self) -> Option<f64>;
  fn set_time(&mut self, time: Option<f64>);
}

/// Trait for edge types that support timetree inference.
/// Provides access to branch length distribution and message passing fields.
pub trait TimetreeEdge: GraphEdge {
  fn branch_length_distribution(&self) -> &Option<Arc<Distribution>>;
  fn set_branch_length_distribution(&mut self, dist: Option<Arc<Distribution>>);
  fn msg_to_parent(&self) -> &Option<Arc<Distribution>>;
  fn set_msg_to_parent(&mut self, msg: Option<Arc<Distribution>>);
}
