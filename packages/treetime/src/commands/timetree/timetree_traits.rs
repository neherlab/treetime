use crate::commands::clock::clock_traits::ClockEdge;
use crate::distribution::distribution::Distribution;
use crate::graph::edge::{BranchDistribution, GraphEdge, TimeLength};
use crate::graph::node::{GraphNode, TimeConstraint};
use std::sync::Arc;

/// Trait for node types that support timetree inference.
/// Provides access to time distribution and estimated time fields.
pub trait TimetreeNode: GraphNode + TimeConstraint<Arc<Distribution>> {
  fn time(&self) -> Option<f64>;
  fn set_time(&mut self, time: Option<f64>);
}

/// Trait for edge types that support timetree inference.
/// Combines clock edge capabilities with branch distribution and time length.
pub trait TimetreeEdge: GraphEdge + ClockEdge + BranchDistribution<Arc<Distribution>> + TimeLength {}
