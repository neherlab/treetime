use std::ops::{Add, Sub};

pub trait RootStats: Clone + Default + Add<Output = Self> + Sub<Output = Self> + Send + Sync {
  fn leaf(time: Option<f64>, branch_length: f64, variance: f64) -> Self;

  #[must_use]
  fn propagate(&self, branch_length: f64, variance: f64) -> Self;

  fn score(&self) -> f64;
}
