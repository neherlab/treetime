use crate::commands::clock::clock_traits::ClockEdge;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_graph::edge::{BranchDistribution, GraphEdge, TimeLength};
use treetime_graph::node::{GraphNode, TimeConstraint};

/// Trait for node types that support timetree inference.
/// Provides access to time distribution and estimated time fields.
pub trait TimetreeNode: GraphNode + TimeConstraint<Arc<Distribution>> {
  fn time(&self) -> Option<f64>;
  fn set_time(&mut self, time: Option<f64>);
}

/// Trait for edge types that support timetree inference.
/// Combines clock edge capabilities with branch distribution and time length.
pub trait TimetreeEdge: GraphEdge + ClockEdge + BranchDistribution<Arc<Distribution>> + TimeLength {
  /// Per-branch relaxed clock rate multiplier. Default 1.0 (strict clock).
  fn gamma(&self) -> f64;
  fn set_gamma(&mut self, gamma: f64);
}
