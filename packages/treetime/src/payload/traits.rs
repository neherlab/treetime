use crate::payload::clock_set::ClockSet;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_graph::edge::{BranchDistribution, ClockMessages, GraphEdge, HasBranchLength, TimeLength};
use treetime_graph::node::{GraphNode, Named, Outlier, TimeConstraint};

pub trait ClockNode: Outlier + Send + Sync {
  fn likely_time(&self) -> Option<f64>;
  fn div(&self) -> f64;
  fn set_div(&mut self, div: f64);
  fn clock_set(&self) -> &ClockSet;
  fn clock_set_mut(&mut self) -> &mut ClockSet;
  fn set_clock_set(&mut self, clock_set: ClockSet) {
    *self.clock_set_mut() = clock_set;
  }
}

pub trait ClockEdge: ClockMessages<ClockSet> + HasBranchLength + TimeLength + Send + Sync {
  /// Per-branch relaxed clock rate multiplier. Default 1.0 (strict clock).
  fn gamma(&self) -> f64 {
    1.0
  }
}

pub trait DateConstraintNode: GraphNode + Named + TimeConstraint<Arc<Distribution>> {}

/// Trait for node types that support timetree inference.
/// Provides access to time distribution and estimated time fields.
pub trait TimetreeNode: GraphNode + TimeConstraint<Arc<Distribution>> {
  fn time(&self) -> Option<f64>;
  fn set_time(&mut self, time: Option<f64>);
}

/// Trait for edge types that support timetree inference.
/// Combines clock edge capabilities with branch distribution and time length.
pub trait TimetreeEdge: GraphEdge + ClockEdge + BranchDistribution<Arc<Distribution>> + TimeLength {
  fn set_gamma(&mut self, gamma: f64);
}
