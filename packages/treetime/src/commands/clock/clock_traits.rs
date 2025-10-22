use crate::commands::clock::clock_set::ClockSet;

pub trait ClockNode: Send + Sync {
  fn likely_time(&self) -> Option<f64>;
  fn div(&self) -> f64;
  fn is_outlier(&self) -> bool;
  fn clock_set(&self) -> &ClockSet;
  fn clock_set_mut(&mut self) -> &mut ClockSet;
  fn set_clock_set(&mut self, clock_set: ClockSet) {
    *self.clock_set_mut() = clock_set;
  }
}

pub trait ClockEdge: Send + Sync {
  fn branch_length(&self) -> Option<f64>;
  fn set_branch_length(&mut self, length: Option<f64>);
  fn to_parent(&self) -> &ClockSet;
  fn to_parent_mut(&mut self) -> &mut ClockSet;
  fn to_child(&self) -> &ClockSet;
  fn to_child_mut(&mut self) -> &mut ClockSet;
  fn from_child(&self) -> &ClockSet;
  fn from_child_mut(&mut self) -> &mut ClockSet;
}
