use crate::commands::timetree::coalescent::events::collect_tree_events;
use crate::commands::timetree::coalescent::integration::compute_integral_merger_rate;
use crate::commands::timetree::coalescent::lineage_dynamics::compute_lineage_count_distribution;
use crate::commands::timetree::coalescent::piecewise_constant_fn::PiecewiseConstantFn;
use crate::commands::timetree::coalescent::time_coordinate::{CalendarTime, Tbp};
use crate::commands::timetree::timetree_traits::TimetreeNode;
use crate::make_report;
use argmin::core::observers::{Observe, ObserverMode};
use argmin::core::{CostFunction, Error, Executor, State};
use argmin::solver::brent::BrentOpt;
use eyre::Report;
use log::{debug, info, warn};
use treetime_distribution::Distribution;
use treetime_graph::edge::{GraphEdge, TimeLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;

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
    .add_observer(TcOptimizationObserver, ObserverMode::Always)
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
/// Computes negative total coalescent likelihood for a given log(Tc) value.
struct TcCostFunction {
  /// Lineage count distribution k(t), precomputed from tree events.
  lineage_counts: PiecewiseConstantFn,
  /// Node branch data for likelihood computation.
  node_branches: Vec<NodeBranchData>,
}

struct NodeBranchData {
  /// Time before present of the node.
  t_node: Tbp,
  /// Branch length to parent.
  branch_length: f64,
  /// Number of children (2 for binary, higher for polytomies).
  multiplicity: f64,
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

    // Precompute lineage count distribution - depends only on tree structure, not Tc
    let lineage_counts = compute_lineage_count_distribution(&events_tbp)?;

    // Collect node branch data for likelihood computation
    let mut node_branches = Vec::new();

    graph.iter_breadth_first_forward(|node| {
      // Only process nodes with parent edges (skip root)
      if node.parent_keys.is_empty() {
        return;
      }

      let Some(time_dist) = node.payload.time_distribution() else {
        return;
      };

      let Some(t_calendar) = time_dist.likely_time() else {
        return;
      };

      let t_node = CalendarTime::new(t_calendar).to_tbp(present_time);

      // Get branch length from parent edge
      // parent_keys is Vec<(GraphNodeKey, GraphEdgeKey)>
      let (parent_node_key, parent_edge_key) = node.parent_keys[0];

      // Try to get time length from edge, or compute from parent-child time difference
      let branch_length = graph
        .get_edge(parent_edge_key)
        .and_then(|e| e.read_arc().payload().read_arc().time_length())
        .or_else(|| {
          // Fall back to computing from parent-child time difference
          graph.get_node(parent_node_key).and_then(|parent| {
            parent
              .read_arc()
              .payload()
              .read_arc()
              .time_distribution()
              .as_ref()
              .and_then(|d| d.likely_time())
              .map(|parent_t| {
                let parent_tbp = CalendarTime::new(parent_t).to_tbp(present_time);
                parent_tbp - t_node
              })
          })
        });

      let Some(branch_length) = branch_length else {
        warn!(
          "Tc optimization: skipping node (key={:?}) with undetermined branch length",
          node.key
        );
        return;
      };

      // Multiplicity: number of children for internal nodes, 2 for leaves
      let multiplicity = if node.child_edges.is_empty() {
        2.0 // Leaves contribute as binary merger
      } else {
        node.child_edges.len() as f64
      };

      node_branches.push(NodeBranchData {
        t_node,
        branch_length,
        multiplicity,
      });
    });

    Ok(Self {
      lineage_counts,
      node_branches,
    })
  }

  /// Compute total coalescent likelihood for a given Tc value.
  ///
  /// This follows Python v0's total_LH() method:
  ///   LH = -Σ cost(t_node, branch_length)
  ///
  /// where cost = I(t_merger) - I(t_node) - log(λ(t_merger)) * (m-1)/m
  fn compute_total_lh(&self, tc: f64) -> Result<f64, Report> {
    let tc_dist = Distribution::constant(tc);
    let integral_merger_rate = compute_integral_merger_rate(&tc_dist, &self.lineage_counts)?;

    let mut total_lh = 0.0;

    for node_data in &self.node_branches {
      let t_node = node_data.t_node;
      let branch_length = node_data.branch_length;
      let multiplicity = node_data.multiplicity;

      // Merger time is when the branch joins its parent
      let t_merger = t_node + branch_length.max(0.0);

      // Compute cost following Python v0's cost() method
      let i_merger = integral_merger_rate.eval(t_merger.value());
      let i_node = integral_merger_rate.eval(t_node.value());

      // Total merger rate λ(t) = k(k-1)/(2*Tc)
      let k = self.lineage_counts.eval(t_merger.value());
      let k_clamped = f64::max(0.5, k - 1.0);
      let lambda = 0.5 * k_clamped * (k_clamped + 1.0) / tc;

      let cost = (i_merger - i_node) - lambda.ln() * (multiplicity - 1.0) / multiplicity;

      total_lh -= cost;
    }

    Ok(total_lh)
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

/// Observer for Tc optimization progress.
struct TcOptimizationObserver;

impl<I> Observe<I> for TcOptimizationObserver
where
  I: State,
  <I as State>::Param: std::fmt::Debug,
  <I as State>::Float: std::fmt::LowerExp,
{
  fn observe_iter(&mut self, state: &I, _kv: &argmin::core::KV) -> Result<(), Error> {
    if state.get_iter().is_multiple_of(10) || state.get_iter() <= 3 {
      debug!(
        "Tc optimization iteration {}: log(Tc) = {:?}, cost = {:.6e}",
        state.get_iter(),
        state.get_best_param(),
        state.get_best_cost()
      );
    }
    Ok(())
  }
}
