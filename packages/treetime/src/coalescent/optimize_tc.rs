use crate::coalescent::edge_data::{CoalescentEdgeData, collect_coalescent_edges, sum_coalescent_cost};
use crate::coalescent::events::collect_tree_events;
use crate::coalescent::integration::compute_integral_merger_rate;
use crate::coalescent::lineage_dynamics::compute_lineage_count_distribution;
use crate::make_report;
use crate::optimize::observer::OptimizationObserver;
use crate::payload::traits::TimetreeNode;
use argmin::core::observers::ObserverMode;
use argmin::core::{CostFunction, Error, Executor};
use argmin::solver::brent::BrentOpt;
use eyre::Report;
use log::{debug, info};
use treetime_distribution::Distribution;
use treetime_graph::edge::{GraphEdge, TimeLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;

/// Result of Tc optimization.
pub struct OptimizeTcResult {
  /// Optimized coalescence time scale.
  pub tc: f64,
  /// Total coalescent likelihood at optimized Tc.
  pub likelihood: f64,
  /// Whether optimization succeeded.
  pub success: bool,
}

/// Optimizes the coalescence time scale Tc to maximize coalescent likelihood.
///
/// Uses Brent's method in log space to find the Tc that maximizes the total
/// coalescent likelihood of the tree. The bracket [-20.0, 2.0] in log space
/// corresponds to roughly [2e-9, 7.4] in linear space.
///
/// # Algorithm
///
/// The total coalescent likelihood is:
///   LH = -Σ cost(t_node, branch_length)
///
/// where cost is the negative log probability of the coalescent process
/// for each branch in the tree.
///
/// # Returns
///
/// Returns the optimized Tc value, or falls back to the initial Tc on failure.
pub fn optimize_tc<N, E, D>(graph: &Graph<N, E, D>, initial_tc: f64) -> Result<OptimizeTcResult, Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + TimeLength,
  D: Sync + Send,
{
  info!("Optimizing coalescent time scale Tc (initial Tc = {initial_tc:.6e})");

  let cost_fn = TcCostFunction::new(graph)?;

  let initial_log_tc = initial_tc.ln();
  debug!("Initial log(Tc) = {initial_log_tc:.4}");

  // Brent's method with bracket in log space
  // Python v0 uses bracket=[-20.0, 2.0]
  let solver = BrentOpt::new(-20.0, 2.0);

  let result = Executor::new(&cost_fn, solver)
    .configure(|cfg| cfg.max_iters(100).target_cost(1e-10))
    .add_observer(
      OptimizationObserver {
        label: "Tc optimization",
        early_threshold: 3,
      },
      ObserverMode::Always,
    )
    .run()
    .map_err(|e| make_report!("Tc optimization failed: {e}"))?;

  let best_log_tc = result.state.best_param.unwrap_or(initial_log_tc);
  let best_cost = result.state.best_cost;
  let optimized_tc = best_log_tc.exp();

  // Check if we got a valid result
  let success = optimized_tc.is_finite() && best_cost.is_finite();

  info!(
    "Tc optimization completed after {} iterations: Tc = {optimized_tc:.6e} (log(Tc) = {best_log_tc:.4}), LH = {:.4}",
    result.state.iter, -best_cost
  );

  Ok(OptimizeTcResult {
    tc: optimized_tc,
    likelihood: -best_cost,
    success,
  })
}

/// Cost function for Tc optimization.
///
/// Precomputes lineage counts and per-edge data from the tree once,
/// then evaluates the total coalescent likelihood for each candidate Tc
/// during Brent's method iterations.
struct TcCostFunction {
  lineage_counts: PiecewiseConstantFn,
  edges: Vec<CoalescentEdgeData>,
}

impl TcCostFunction {
  fn new<N, E, D>(graph: &Graph<N, E, D>) -> Result<Self, Report>
  where
    N: GraphNode + TimetreeNode,
    E: GraphEdge + TimeLength,
    D: Sync + Send,
  {
    let (present_time, events_calendar) = collect_tree_events(graph)?;
    let events_tbp: Vec<_> = events_calendar
      .iter()
      .map(|(t, delta)| (t.to_tbp(present_time), *delta))
      .collect();

    let lineage_counts = compute_lineage_count_distribution(&events_tbp)?;
    let edges = collect_coalescent_edges(graph, present_time);

    Ok(Self { lineage_counts, edges })
  }

  fn compute_total_lh(&self, tc: f64) -> Result<f64, Report> {
    let tc_dist = Distribution::constant(tc);
    let integral_merger_rate = compute_integral_merger_rate(&tc_dist, &self.lineage_counts)?;
    sum_coalescent_cost(&self.edges, &integral_merger_rate, &self.lineage_counts, &tc_dist)
  }
}

impl CostFunction for &TcCostFunction {
  type Param = f64;
  type Output = f64;

  fn cost(&self, log_tc: &Self::Param) -> Result<Self::Output, Error> {
    let tc = log_tc.exp();

    // Return negative likelihood (minimization problem)
    let lh = self
      .compute_total_lh(tc)
      .map_err(|e| Error::msg(format!("Failed to compute likelihood: {e}")))?;

    Ok(-lh)
  }
}
