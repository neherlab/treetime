use crate::commands::clock::clock_graph::ClockGraph;
use crate::commands::clock::clock_regression::ClockOptions;
use crate::commands::clock::find_best_root::cost_function::BranchPointCostFunction;
use crate::commands::clock::find_best_root::find_best_split::FindRootResult;
use crate::commands::clock::find_best_root::params::GridSearchParams;
use crate::graph::edge::{GraphEdgeKey, Weighted};
use eyre::Report;
use ndarray::Array1;

/// Grid search optimization for finding the best split point along an edge
pub fn optimize_grid_search(
  graph: &ClockGraph,
  edge: GraphEdgeKey,
  options: &ClockOptions,
  params: &GridSearchParams,
) -> Result<FindRootResult, Report> {
  let edge = graph.get_edge(edge).expect("Edge not found");
  let edge_payload = edge.read_arc().payload().read_arc();
  let target_node = graph.get_node(edge.read_arc().target()).unwrap();
  let target_node_payload = target_node.read_arc().payload().read_arc();
  let is_leaf = target_node.read_arc().is_leaf();
  let node_date = target_node_payload.date;
  let branch_length = edge_payload.weight().expect("Encountered an edge without a weight");
  let branch_variance = options.variance_factor * branch_length + options.variance_offset;

  // Create cost function
  let cost_fn = BranchPointCostFunction {
    edge_payload: edge_payload.clone(),
    branch_length,
    branch_variance,
    is_leaf,
    node_date,
    options,
  };

  // Grid search - interrogate different positions along the branch
  let mut best_chisq = f64::INFINITY;
  let mut best_split = f64::NAN;
  let mut best_total = cost_fn
    .evaluate_clock_set(0.0)
    .expect("Failed to evaluate initial position");

  for x in Array1::linspace(0.0, 1.0, params.n_points) {
    let clock_total = cost_fn.evaluate_clock_set(x).expect("Failed to evaluate clock set");
    let chisq = clock_total.chisq();
    if chisq < best_chisq {
      best_split = x;
      best_total = clock_total;
      best_chisq = chisq;
    }
  }

  Ok(FindRootResult {
    edge: Some(edge.read_arc().key()),
    split: best_split,
    total: best_total,
    chisq: best_chisq,
  })
}
