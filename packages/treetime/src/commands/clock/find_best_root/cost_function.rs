use crate::commands::clock::clock_graph::ClockEdgePayload;
use crate::commands::clock::clock_regression::ClockOptions;
use crate::commands::clock::clock_set::ClockSet;
use argmin::core::{CostFunction, Error};
use eyre::Report;

/// Cost function for branch point optimization using various optimization methods
pub struct BranchPointCostFunction<'a> {
  pub edge_payload: &'a ClockEdgePayload,
  pub branch_length: f64,
  pub branch_variance: f64,
  pub is_leaf: bool,
  pub node_date: Option<f64>,
  pub options: &'a ClockOptions,
}

impl BranchPointCostFunction<'_> {
  /// Evaluate the clock set at a given position (used to get the final result)
  pub fn evaluate_clock_set(&self, x: f64) -> Result<ClockSet, Report> {
    // determine contribution of child/target first -- terminal nodes need special handling
    let child_contribution = if self.is_leaf {
      ClockSet::leaf_contribution_to_parent(
        self.node_date,
        self.branch_length * (1.0 - x),
        self.branch_variance * (1.0 - x) + self.options.variance_offset_leaf,
      )
    } else {
      self
        .edge_payload
        .to_parent
        .propagate_averages(self.branch_length * (1.0 - x), self.branch_variance * (1.0 - x))
    };

    let clock_total = self
      .edge_payload
      .to_child
      .propagate_averages(self.branch_length * x, self.branch_variance * x)
      + child_contribution;

    Ok(clock_total)
  }
}



impl CostFunction for &BranchPointCostFunction<'_> {
  type Param = f64;
  type Output = f64;

  fn cost(&self, x: &Self::Param) -> Result<Self::Output, Error> {
    // Ensure x is within bounds
    if *x < 0.0 || *x > 1.0 {
      return Ok(f64::INFINITY);
    }

    // Evaluate the clock set and return chi-squared
    Ok(
      self
        .evaluate_clock_set(*x)
        .map_or(f64::INFINITY, |clock_total| clock_total.chisq()),
    )
  }
}
