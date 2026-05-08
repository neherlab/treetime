use crate::commands::clock::clock_regression::ClockParams;
use crate::representation::payload::clock_set::ClockSet;
use crate::representation::payload::traits::{ClockEdge, ClockNode};
use crate::commands::clock::find_best_root::cost_function::BranchPointCostFunction;
use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::commands::clock::find_best_root::{method_brent, method_golden_section, method_grid_search};
use eyre::Report;
use serde::{Deserialize, Serialize};
use treetime_graph::edge::{GraphEdge, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;

#[derive(Debug, Serialize, Deserialize)]
pub struct FindRootResult {
  pub edge: Option<GraphEdgeKey>,

  /// The best root might be somewhere half-way on an existing edge.
  /// The position of this best root is the split. If this is 0 or 1, then we reroot on either the parent (source) or
  /// the child (target) of the edge. If this is 0 < x < 1, we put a new node at that point, and reroot on that new node.
  pub split: f64,

  pub clock_set: ClockSet,

  pub chisq: f64,
}

/// Find the best split point along an edge using the specified optimization method
pub fn find_best_split<N, E, D>(
  graph: &Graph<N, E, D>,
  edge: GraphEdgeKey,
  options: &ClockParams,
  params: &BranchPointOptimizationParams,
) -> Result<FindRootResult, Report>
where
  N: GraphNode + ClockNode,
  E: GraphEdge + ClockEdge,
  D: Send + Sync,
{
  // Create cost function once
  let cost_fn = BranchPointCostFunction::new(graph, edge, options)?;

  match params {
    BranchPointOptimizationParams::Grid(params) => method_grid_search::optimize_grid_search(edge, &cost_fn, params),
    BranchPointOptimizationParams::Brent(params) => method_brent::optimize_brent(edge, &cost_fn, params),
    BranchPointOptimizationParams::GoldenSection(params) => {
      method_golden_section::optimize_golden_section(edge, &cost_fn, params)
    },
  }
}
