use crate::clock::clock_regression::ClockVarianceParams;
use crate::clock::clock_set::ClockSet;
use crate::clock::clock_state::{ClockInputs, ClockState};
use crate::clock::find_best_root::params::RootObjective;
use argmin::core::{CostFunction, Error};
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_utils::make_report;

/// Cost function for branch point optimization using various optimization methods
pub struct BranchPointCostFunction<'a> {
  pub to_parent: ClockSet,
  pub to_child: ClockSet,
  pub branch_length: f64,
  pub branch_variance: f64,
  pub is_leaf: bool,
  pub node_time: Option<f64>,
  pub options: &'a ClockVarianceParams,
  pub objective: RootObjective,
}

impl<'a> BranchPointCostFunction<'a> {
  pub fn new(
    graph: &Graph,
    inputs: &ClockInputs,
    state: &ClockState,
    edge: GraphEdgeKey,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    options: &'a ClockVarianceParams,
    objective: RootObjective,
  ) -> Result<BranchPointCostFunction<'a>, Report> {
    let edge_obj = graph
      .get_edge(edge)
      .ok_or_else(|| make_report!("Edge not found: {edge}"))?;
    let target_key = edge_obj.target();
    let target_node = graph
      .get_node(target_key)
      .ok_or_else(|| make_report!("Target node not found for edge: {edge}"))?;
    let is_leaf = target_node.is_leaf();
    let node_time = inputs.likely_time(target_key);
    let branch_length = branch_lengths[&edge].ok_or_else(|| make_report!("Edge {edge} has no weight"))?;
    let branch_variance = options.variance_factor * branch_length + options.variance_offset;

    let edge_state = state.edge(edge);
    Ok(BranchPointCostFunction {
      to_parent: edge_state.clock_to_parent.clone(),
      to_child: edge_state.clock_to_child.clone(),
      branch_length,
      branch_variance,
      is_leaf,
      node_time,
      options,
      objective,
    })
  }

  /// Evaluate the clock set at a given position (used to get the final result)
  pub fn evaluate_clock_set(&self, x: f64) -> Result<ClockSet, Report> {
    // determine contribution of child/target first -- terminal nodes need special handling
    let child_contribution = if self.is_leaf {
      ClockSet::leaf_contribution_to_parent(
        self.node_time,
        self.branch_length * (1.0 - x),
        self.branch_variance * (1.0 - x) + self.options.variance_offset_leaf,
      )
    } else {
      self
        .to_parent
        .propagate_averages(self.branch_length * (1.0 - x), self.branch_variance * (1.0 - x))
    };

    let clock_set = self
      .to_child
      .propagate_averages(self.branch_length * x, self.branch_variance * x)
      + child_contribution;

    Ok(clock_set)
  }

  pub fn score_clock_set(&self, clock_set: &ClockSet) -> f64 {
    self.objective.score(clock_set)
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

    // Evaluate the clock set and return the configured objective value.
    let result = self
      .evaluate_clock_set(*x)
      .map_or(f64::INFINITY, |clock_set| self.score_clock_set(&clock_set));

    Ok(result)
  }
}
