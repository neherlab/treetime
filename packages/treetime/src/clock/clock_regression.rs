use crate::clock::clock_model::{ClockModel, ClockRegression};
use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RootObjective};
use crate::clock::reroot::{RerootParams, reroot_in_place};
use crate::payload::clock_set::ClockSet;
use crate::payload::traits::{ClockEdge, ClockNode};
use eyre::Report;
use log::{debug, info};
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use treetime_graph::breadth_first::GraphTraversalContinuation;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_graph::reroot::RerootResult;

#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct ClockParams {
  /// Variance scaling factor proportional to branch length
  #[cfg_attr(feature = "clap", clap(long, default_value_t = ClockParams::default().variance_factor))]
  #[default = 0.0]
  pub variance_factor: f64,

  /// Constant variance offset for all branches
  #[cfg_attr(feature = "clap", clap(long, default_value_t = ClockParams::default().variance_offset))]
  #[default = 0.0]
  pub variance_offset: f64,

  /// Additional variance offset for leaf (terminal) nodes
  #[cfg_attr(feature = "clap", clap(long, default_value_t = ClockParams::default().variance_offset_leaf))]
  #[default = 1.0]
  pub variance_offset_leaf: f64,
}

/// Result of clock estimation with optional rerooting.
///
/// Contains either a raw regression result (estimated rate, any sign) or a
/// validated `ClockModel` (fixed rate, positive). Callers that need a validated
/// model call `into_clock_model()` at their boundary.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ClockRerootResult {
  regression: Option<ClockRegression>,
  clock_model: Option<ClockModel>,
  reroot_result: Option<RerootResult>,
}

impl ClockRerootResult {
  pub fn into_clock_model(self) -> Result<ClockModel, Report> {
    if let Some(model) = self.clock_model {
      return Ok(model);
    }
    let regression = self
      .regression
      .expect("ClockRerootResult has neither regression nor clock_model");
    ClockModel::from_regression(&regression)
  }

  pub fn regression(&self) -> &ClockRegression {
    self
      .regression
      .as_ref()
      .expect("regression() called on fixed-rate result")
  }

  pub fn reroot_result(&self) -> Option<&RerootResult> {
    self.reroot_result.as_ref()
  }
}

/// Runs backward clock regression pass.
///
/// `prev_clock_rate`: when `Some(rate)`, uses solver-updated `time_length * rate * gamma`
/// as divergence (re-estimation mode). When `None`, uses input `branch_length()` (initial estimation).
pub fn clock_regression_backward<N, E, D>(
  graph: &Graph<N, E, D>,
  options: &ClockParams,
  prev_clock_rate: Option<f64>,
) -> Result<(), Report>
where
  N: GraphNode + ClockNode,
  E: GraphEdge + ClockEdge,
  D: Send + Sync,
{
  graph.par_iter_breadth_first_backward(|mut n| {
    let date = n.payload.likely_time();
    let q_to_parent = if n.is_leaf {
      if n.payload.is_outlier() {
        ClockSet::outlier_contribution()
      } else {
        ClockSet::leaf_contribution(date)
      }
    } else {
      let mut q_dest = ClockSet::default();
      for (_, e) in &n.children {
        q_dest += e.read_arc().from_child();
      }
      q_dest
    };

    let is_leaf = n.is_leaf;
    let is_root = n.is_root;

    if is_root {
      *n.payload.clock_set_mut() = q_to_parent;
    } else {
      let edge_to_parent = n.get_exactly_one_parent_edge()?;
      *edge_to_parent.to_parent_mut() = q_to_parent;
      let edge_len = edge_divergence(
        edge_to_parent.branch_length(),
        edge_to_parent.time_length(),
        edge_to_parent.gamma(),
        prev_clock_rate,
      );
      let mut branch_variance = options.variance_factor * edge_len + options.variance_offset;
      *edge_to_parent.from_child_mut() = if is_leaf {
        branch_variance += options.variance_offset_leaf;
        ClockSet::leaf_contribution_to_parent(date, edge_len, branch_variance)
      } else {
        edge_to_parent.to_parent().propagate_averages(edge_len, branch_variance)
      };
    }
    Ok(GraphTraversalContinuation::Continue)
  })
}

/// Runs forward clock regression pass.
///
/// `prev_clock_rate`: when `Some(rate)`, uses solver-updated `time_length * rate * gamma`
/// as divergence (re-estimation mode). When `None`, uses input `branch_length()` (initial estimation).
pub fn clock_regression_forward<N, E, D>(
  graph: &Graph<N, E, D>,
  options: &ClockParams,
  prev_clock_rate: Option<f64>,
) -> Result<(), Report>
where
  N: GraphNode + ClockNode,
  E: GraphEdge + ClockEdge,
  D: Sync + Send,
{
  graph.par_iter_breadth_first_forward(|mut n| {
    if !n.is_root {
      let (_, edge) = n.get_exactly_one_parent()?;
      let edge = edge.read_arc();
      let edge_len = edge_divergence(edge.branch_length(), edge.time_length(), edge.gamma(), prev_clock_rate);
      let branch_variance = options.variance_factor * edge_len + options.variance_offset;

      let mut q_dest = edge.to_parent().clone();
      q_dest += edge.to_child().propagate_averages(edge_len, branch_variance);
      *n.payload.clock_set_mut() = q_dest;
    }

    for mut child_edge in n.child_edges {
      let mut q = n.payload.clock_set().clone();
      q -= child_edge.from_child();
      *child_edge.to_child_mut() = q;
    }
    Ok(GraphTraversalContinuation::Continue)
  })
}

/// Estimates clock model with optional rerooting using default policy.
///
/// `prev_clock_rate`: when `Some(rate)`, regression uses solver-updated time lengths
/// converted to divergence (re-estimation mode). When `None`, uses input branch lengths.
pub fn estimate_clock_model_with_reroot<N, E, D>(
  graph: &mut Graph<N, E, D>,
  options: &ClockParams,
  clock_rate: Option<f64>,
  keep_root: bool,
  optimization_params: &BranchPointOptimizationParams,
  prev_clock_rate: Option<f64>,
) -> Result<ClockModel, Report>
where
  N: GraphNode + ClockNode + Named + Default,
  E: GraphEdge + ClockEdge + Default,
  D: Send + Sync,
{
  let reroot_params = RerootParams::default();
  let result = estimate_clock_model_with_reroot_policy(
    graph,
    options,
    clock_rate,
    keep_root,
    optimization_params,
    &reroot_params,
    prev_clock_rate,
  )?;
  result.into_clock_model()
}

/// Estimates clock model with optional rerooting using explicit policy.
///
/// `prev_clock_rate`: when `Some(rate)`, regression uses solver-updated time lengths
/// converted to divergence (re-estimation mode). When `None`, uses input branch lengths.
pub fn estimate_clock_model_with_reroot_policy<N, E, D>(
  graph: &mut Graph<N, E, D>,
  options: &ClockParams,
  clock_rate: Option<f64>,
  keep_root: bool,
  optimization_params: &BranchPointOptimizationParams,
  reroot_params: &RerootParams,
  prev_clock_rate: Option<f64>,
) -> Result<ClockRerootResult, Report>
where
  N: GraphNode + ClockNode + Named + Default,
  E: GraphEdge + ClockEdge + Default,
  D: Send + Sync,
{
  if let Some(rate) = clock_rate {
    info!("## Estimating clock model with fixed rate {rate:.6e} (keep_root={keep_root})");
  } else {
    info!("## Estimating clock model (keep_root={keep_root})");
  }

  info!("### Running backward regression");
  clock_regression_backward(graph, options, prev_clock_rate)?;
  debug!("Backward regression completed");

  let reroot_result = if !keep_root {
    info!("### Running forward regression to find optimal root");
    clock_regression_forward(graph, options, prev_clock_rate)?;
    debug!("Forward regression completed");

    info!("### Finding best root and rerooting tree");
    let reroot_params = clock_rate.map_or_else(
      || reroot_params.clone(),
      |rate| reroot_params.with_objective(RootObjective::FixedRate(rate)),
    );
    let reroot_result = reroot_in_place(graph, options, optimization_params, &reroot_params)?;
    info!("Rerooted to node {}", reroot_result.new_root_key.0);
    debug!("Rerooting completed");
    Some(reroot_result)
  } else {
    info!("### Keeping original root (--keep-root enabled)");
    None
  };

  info!("### Extracting clock model from root");
  let root = graph.get_exactly_one_root()?;
  let root = root.read_arc().payload().read_arc();

  let (regression, clock_model) = if let Some(rate) = clock_rate {
    info!("### Using fixed clock rate: {rate:.6e}");
    (None, Some(ClockModel::with_fixed_rate(root.clock_set(), rate)?))
  } else {
    info!("### Using estimated clock rate");
    let regression = ClockRegression::from_clock_set(root.clock_set())?;
    (Some(regression), None)
  };

  let rate = clock_model
    .as_ref()
    .map_or_else(|| regression.as_ref().unwrap().clock_rate(), |m| m.clock_rate());
  let intercept = clock_model
    .as_ref()
    .map_or_else(|| regression.as_ref().unwrap().intercept(), |m| m.intercept());
  info!("**Clock rate:** {rate:.6e}");
  info!("**Intercept:** {intercept:.4}");
  if let Some(reg) = &regression {
    info!("**R²:** {:.4}", reg.r_val() * reg.r_val());
    info!("**χ²:** {:.4}", reg.chisq());
    info!("**Hessian:**\n{}", reg.hessian());
  }

  Ok(ClockRerootResult {
    regression,
    clock_model,
    reroot_result,
  })
}

/// Compute divergence (substitutions/site) for an edge.
///
/// In re-estimation mode (`prev_clock_rate` is `Some`), converts solver-updated time length
/// back to divergence: `time_length * rate * gamma`. Falls back to input `branch_length`
/// when `time_length` is not yet populated (initial estimation or pre-solver edges).
fn edge_divergence(
  branch_length: Option<f64>,
  time_length: Option<f64>,
  gamma: f64,
  prev_clock_rate: Option<f64>,
) -> f64 {
  if let Some(rate) = prev_clock_rate {
    if let Some(tl) = time_length {
      return tl * rate * gamma;
    }
  }
  branch_length.expect("Encountered an edge without a weight")
}
