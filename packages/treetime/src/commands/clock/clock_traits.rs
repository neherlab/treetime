use crate::commands::clock::clock_set::ClockSet;
use treetime_graph::edge::{ClockMessages, HasBranchLength};
use treetime_graph::node::Outlier;

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

pub trait ClockEdge: ClockMessages<ClockSet> + HasBranchLength + Send + Sync {}
