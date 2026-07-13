use crate::reroot::traits::RootStats;
use argmin::core::{CostFunction, Error};

/// Cost function for optimizing the root position along a single edge.
///
/// Holds the directional statistics on the edge -- `to_parent` (the child
/// subtree as seen at the child node) and `to_child` (the rest of the tree as
/// seen at the parent node) -- plus the branch geometry. Evaluating at split
/// fraction `x` places the root at distance `x * branch_length` from the source
/// (parent) end: `x = 0` roots at the target node, `x = 1` at the source node.
///
/// The variance is split linearly across the branch (the full-branch variance
/// scaled by the split fraction), matching the v0 branch-point objective. Leaf
/// branches add the leaf variance offset on the child side after scaling.
pub struct EdgeCostFn<S: RootStats> {
  pub to_parent: S,
  pub to_child: S,
  pub branch_length: f64,
  /// Full-branch variance excluding the leaf offset (`VarianceModel::branch`).
  pub branch_variance: f64,
  pub is_leaf: bool,
  /// Tip date for the target node when it is a leaf; `None` for date-free objectives.
  pub leaf_time: Option<f64>,
  pub variance_offset_leaf: f64,
}

impl<S: RootStats> EdgeCostFn<S> {
  /// Combine the child-side and parent-side statistics at split fraction `x`.
  pub fn evaluate(&self, x: f64) -> S {
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
