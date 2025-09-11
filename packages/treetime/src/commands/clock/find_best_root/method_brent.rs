use crate::commands::clock::clock_graph::ClockGraph;
use crate::commands::clock::clock_regression::ClockOptions;
use crate::commands::clock::find_best_root::cost_function::BranchPointCostFunction;
use crate::commands::clock::find_best_root::find_best_split::FindRootResult;
use crate::commands::clock::find_best_root::params::BrentParams;
use crate::graph::edge::{GraphEdgeKey, Weighted};
use crate::make_report;
use argmin::core::Executor;
use argmin::solver::brent::BrentOpt;
use eyre::Report;

/// Brent's method optimization for finding the best split point along an edge
pub fn optimize_brent(
  graph: &ClockGraph,
  edge: GraphEdgeKey,
  options: &ClockOptions,
  params: &BrentParams,
) -> Result<FindRootResult, Report> {
  let edge_obj = graph.get_edge(edge).expect("Edge not found");
  let edge_payload = edge_obj.read_arc().payload().read_arc();
  let target_node = graph.get_node(edge_obj.read_arc().target()).unwrap();
  let target_node_payload = target_node.read_arc().payload().read_arc();
  let is_leaf = target_node.read_arc().is_leaf();
  let node_date = target_node_payload.date;
  let branch_length = edge_payload.weight().expect("Encountered an edge without a weight");
  let branch_variance = options.variance_factor * branch_length + options.variance_offset;

  let cost_fn = BranchPointCostFunction {
    edge_payload: edge_payload.clone(),
    branch_length,
    branch_variance,
    is_leaf,
    node_date,
    options,
  };

  // Set up Brent solver with bounds [0.0, 1.0]
  // 0.0 means placing the root at the target node, 1.0 means placing it at the source node.
  let solver = BrentOpt::new(0.0, 1.0);

  // Run optimization
  let result = Executor::new(cost_fn, solver)
    .configure(|cfg| cfg.max_iters(params.max_iters as u64).target_cost(params.tolerance))
    .run()
    .map_err(|e| make_report!("Brent optimization failed: {}", e))?;

  let best_split = result.state.best_param.unwrap();
  let best_chisq = result.state.best_cost;

  // Evaluate the cost function one more time to get the ClockSet data
  let cost_fn = BranchPointCostFunction {
    edge_payload: edge_payload.clone(),
    branch_length,
    branch_variance,
    is_leaf,
    node_date,
    options,
  };

  let best_total = cost_fn.evaluate_clock_set(best_split)?;

  Ok(FindRootResult {
    edge: Some(edge),
    split: best_split,
    total: best_total,
    chisq: best_chisq,
  })
}
