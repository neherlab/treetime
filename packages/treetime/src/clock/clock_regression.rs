use crate::clock::clock_model::{ClockModel, ClockRegression};
use crate::clock::clock_set::ClockSet;
use crate::clock::clock_state::{ClockEdgeState, ClockInputs, ClockNodeState, ClockState};
use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RootObjective};
use crate::clock::reroot::{RerootParams, reroot_in_place};
use eyre::Report;
use log::{debug, info};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::fmt::Debug;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pass::{GraphPassBackwardContext, GraphPassNodeOutput};
use treetime_graph::reroot::RerootResult;

#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct ClockParams {
  /// Variance scaling factor proportional to branch length
  #[default = 0.0]
  pub variance_factor: f64,

  /// Constant variance offset for all branches
  #[default = 0.0]
  pub variance_offset: f64,

  /// Additional variance offset for leaf (terminal) nodes
  #[default = 1.0]
  pub variance_offset_leaf: f64,
}

/// Result of clock estimation with optional rerooting.
///
/// Contains either a raw regression result (estimated rate, any sign) or a
/// validated `ClockModel` (fixed rate, positive). Callers that require a positive
/// rate for time inference (timetree) call `into_clock_model()`, which errors on a
/// non-positive rate. The clock command, which only reports the regression, calls
/// `into_clock_model_allow_negative()`, which warns and continues.
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

  /// Like `into_clock_model`, but permits a non-positive estimated rate (warning
  /// instead of error). Used by the clock command; see `from_regression_allow_negative`.
  pub fn into_clock_model_allow_negative(self) -> ClockModel {
    if let Some(model) = self.clock_model {
      return model;
    }
    let regression = self
      .regression
      .expect("ClockRerootResult has neither regression nor clock_model");
    ClockModel::from_regression_allow_negative(&regression)
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
pub fn clock_regression_backward(
  graph: &Graph,
  inputs: &ClockInputs,
  state: &mut ClockState,
  options: &ClockParams,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  prev_clock_rate: Option<f64>,
) -> Result<(), Report> {
  state.map_backward(graph, |context| {
    clock_regression_backward_node(options, prev_clock_rate, branch_lengths, inputs, &context)
  })
}

fn clock_regression_backward_node(
  options: &ClockParams,
  prev_clock_rate: Option<f64>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  inputs: &ClockInputs,
  context: &GraphPassBackwardContext<'_, ClockNodeState, ClockEdgeState, ClockNodeState, ClockEdgeState>,
) -> Result<GraphPassNodeOutput<ClockNodeState, ClockEdgeState>, Report> {
  let mut node = context.input.clone();
  let is_leaf = context.is_leaf;
  let date = inputs.likely_time(context.key);
  let q_to_parent = if is_leaf {
    if node.is_outlier {
      ClockSet::outlier_contribution()
    } else {
      ClockSet::leaf_contribution(date)
    }
  } else {
    // Children arrive in the graph's canonical `children_of` order, so the moment sums fold in that
    // order, keeping the floating-point result byte-for-byte identical.
    context.children.iter().fold(ClockSet::default(), |mut total, child| {
      let edge = child.edge.expect("Non-root indexed node must own its parent edge");
      total += &edge.clock_from_child;
      total
    })
  };

  let parent_message = if let Some((edge_key, edge)) = context.parent_edge {
    let mut edge = edge.clone();
    edge.clock_to_parent = q_to_parent;
    let edge_input = inputs.edge(edge_key);
    let edge_len = edge_divergence(
      branch_lengths[&edge_key],
      edge_input.time_length,
      edge_input.gamma,
      prev_clock_rate,
    );
    let mut branch_variance = options.variance_factor * edge_len + options.variance_offset;
    edge.clock_from_child = if is_leaf {
      branch_variance += options.variance_offset_leaf;
      ClockSet::leaf_contribution_to_parent(date, edge_len, branch_variance)
    } else {
      edge.clock_to_parent.propagate_averages(edge_len, branch_variance)
    };
    Some(edge)
  } else {
    node.clock_set = q_to_parent;
    None
  };

  Ok(GraphPassNodeOutput { node, parent_message })
}

/// Runs forward clock regression pass.
///
/// `prev_clock_rate`: when `Some(rate)`, uses solver-updated `time_length * rate * gamma`
/// as divergence (re-estimation mode). When `None`, uses input `branch_length()` (initial estimation).
pub fn clock_regression_forward(
  graph: &Graph,
  inputs: &ClockInputs,
  state: &mut ClockState,
  options: &ClockParams,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  prev_clock_rate: Option<f64>,
) -> Result<(), Report> {
  state.map_forward(graph, |context| {
    let mut node = context.input.clone();
    let parent_message = if let Some((edge_key, edge)) = context.parent_edge {
      let mut edge = edge.clone();
      let parent = context.parent.expect("Non-root node must have a parent");
      let mut q_to_child = parent.clock_set.clone();
      q_to_child -= &edge.clock_from_child;
      edge.clock_to_child = q_to_child;

      let edge_input = inputs.edge(edge_key);
      let edge_len = edge_divergence(
        branch_lengths[&edge_key],
        edge_input.time_length,
        edge_input.gamma,
        prev_clock_rate,
      );
      let branch_variance = options.variance_factor * edge_len + options.variance_offset;
      let mut q_dest = edge.clock_to_parent.clone();
      q_dest += edge.clock_to_child.propagate_averages(edge_len, branch_variance);
      node.clock_set = q_dest;
      Some(edge)
    } else {
      None
    };
    Ok(GraphPassNodeOutput { node, parent_message })
  })
}

/// Estimates clock model with optional rerooting using explicit policy.
///
/// `prev_clock_rate`: when `Some(rate)`, regression uses solver-updated time lengths
/// converted to divergence (re-estimation mode). When `None`, uses input branch lengths.
pub fn estimate_clock_model_with_reroot_policy(
  graph: &mut Graph,
  inputs: &mut ClockInputs,
  mut state: ClockState,
  options: &ClockParams,
  clock_rate: Option<f64>,
  keep_root: bool,
  optimization_params: &BranchPointOptimizationParams,
  reroot_params: &RerootParams,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  prev_clock_rate: Option<f64>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<(ClockState, ClockRerootResult), Report> {
  if let Some(rate) = clock_rate {
    info!("## Estimating clock model with fixed rate {rate:.6e} (keep_root={keep_root})");
  } else {
    info!("## Estimating clock model (keep_root={keep_root})");
  }

  info!("### Running backward regression");
  clock_regression_backward(graph, inputs, &mut state, options, branch_lengths, prev_clock_rate)?;
  debug!("Backward regression completed");

  let reroot_result = if !keep_root {
    info!("### Running forward regression to find optimal root");
    clock_regression_forward(graph, inputs, &mut state, options, branch_lengths, prev_clock_rate)?;
    debug!("Forward regression completed");

    info!("### Finding best root and rerooting tree");
    let reroot_params = clock_rate.map_or_else(
      || reroot_params.clone(),
      |rate| reroot_params.with_objective(RootObjective::FixedRate(rate)),
    );
    let (new_state, reroot_result) = reroot_in_place(
      graph,
      inputs,
      state,
      options,
      optimization_params,
      &reroot_params,
      branch_lengths,
      names,
    )?;
    state = new_state;
    info!("Rerooted to node {}", reroot_result.new_root_key.0);
    debug!("Rerooting completed");
    Some(reroot_result)
  } else {
    info!("### Keeping original root (--keep-root enabled)");
    None
  };

  info!("### Extracting clock model from root");
  let root_key = graph.get_exactly_one_root()?.key();
  let root_clock_set = state.node(root_key).clock_set.clone();

  let (regression, clock_model) = if let Some(rate) = clock_rate {
    info!("### Using fixed clock rate: {rate:.6e}");
    (None, Some(ClockModel::with_fixed_rate(&root_clock_set, rate)?))
  } else {
    info!("### Using estimated clock rate");
    let regression = ClockRegression::from_clock_set(&root_clock_set)?;
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

  Ok((
    state,
    ClockRerootResult {
      regression,
      clock_model,
      reroot_result,
    },
  ))
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
