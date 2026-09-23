use crate::reroot::traits::RootStats;
use argmin::core::{CostFunction, Error};

impl<S: RootStats> CostFunction for &EdgeCostFn<S> {
  type Param = f64;
  type Output = f64;

  fn cost(&self, x: &Self::Param) -> Result<Self::Output, Error> {
    if *x < 0.0 || *x > 1.0 {
      return Ok(f64::INFINITY);
    }
    Ok(self.evaluate(*x).score())
  }
}

pub struct EdgeCostFn<S: RootStats> {
  pub to_parent: S,
  pub to_child: S,
  pub branch_length: f64,
  pub branch_variance: f64,
  pub is_leaf: bool,
  pub leaf_time: Option<f64>,
  pub variance_offset_leaf: f64,
}

impl<S: RootStats> EdgeCostFn<S> {
  pub(crate) fn evaluate(&self, x: f64) -> S {
    let child = if self.is_leaf {
      S::leaf(
        self.leaf_time,
        self.branch_length * (1.0 - x),
        self.branch_variance * (1.0 - x) + self.variance_offset_leaf,
      )
    } else {
      self
        .to_parent
        .propagate(self.branch_length * (1.0 - x), self.branch_variance * (1.0 - x))
    };

    let parent = self
      .to_child
      .propagate(self.branch_length * x, self.branch_variance * x);

    parent + child
  }
}
