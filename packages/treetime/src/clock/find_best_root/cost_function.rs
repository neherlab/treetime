use crate::clock::clock_regression::ClockParams;
use crate::clock::clock_state::ClockState;
use crate::clock::find_best_root::params::RootObjective;
use crate::payload::clock_set::ClockSet;
use argmin::core::{CostFunction, Error};
use eyre::Report;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_utils::make_report;

/// Cost function for branch point optimization using various optimization methods
pub struct BranchPointCostFunction<'a> {
  pub to_parent: ClockSet,
  pub to_child: ClockSet,
  pub branch_length: f64,
  pub branch_variance: f64,
  pub is_leaf: bool,
  pub node_time: Option<f64>,
  pub options: &'a ClockParams,
  pub objective: RootObjective,
}

impl<'a> BranchPointCostFunction<'a> {
  pub fn new<N, E, D>(
    graph: &Graph<N, E, D>,
    state: &ClockState,
    edge: GraphEdgeKey,
    options: &'a ClockParams,
    objective: RootObjective,
  ) -> Result<BranchPointCostFunction<'a>, Report>
  where
    N: GraphNode,
    E: GraphEdge + HasBranchLength,
    D: Send + Sync,
  {
    let edge_obj = graph
      .get_edge(edge)
      .ok_or_else(|| make_report!("Edge not found: {edge}"))?;
    let target_key = edge_obj.read_arc().target();
    let target_node = graph
      .get_node(target_key)
      .ok_or_else(|| make_report!("Target node not found for edge: {edge}"))?;
    let is_leaf = target_node.read_arc().is_leaf();
    let node_time = state.node(target_key).likely_time();
    let branch_length = edge_obj
      .read_arc()
      .payload()
      .read_arc()
      .branch_length()
      .ok_or_else(|| make_report!("Edge {edge} has no weight"))?;
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
