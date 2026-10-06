use crate::clock::clock_model::{ClockModel, ClockRegression};
use crate::clock::clock_set::ClockSet;
use crate::clock::clock_state::{ClockEdgeState, ClockInputs, ClockState};
use crate::clock::divergence::root_to_node_divergences;
use crate::clock::find_best_root::find_best_split::FindRootResult;
use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RootObjective};
use crate::clock::reroot::{RerootParams, remove_undated_stem, reroot_clock_tree, select_root};
use crate::node_label::node_label;
use crate::progress::LogSink;
use crate::progress_info;
use crate::reroot::placement::root_moves;
use eyre::Report;
use log::debug;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Debug;
use std::mem;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pass::{GraphPassBackwardContext, GraphPassNodeOutput};
use treetime_graph::reroot::{RerootResult, StemRemovalInfo};

#[allow(
  clippy::unwrap_used,
  reason = "unwrap on a value an upstream invariant guarantees is present"
)]
pub(crate) fn estimate_clock_model_with_reroot_policy(
  tree: ClockTree,
  outliers: &BTreeSet<GraphNodeKey>,
  options: &ClockVarianceParams,
  clock_rate: Option<f64>,
  keep_root: bool,
  optimization_params: &BranchPointOptimizationParams,
  reroot_params: &RerootParams,
  prev_clock_rate: Option<f64>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> Result<(ClockTree, ClockRerootResult), Report> {
  if let Some(rate) = clock_rate {
    progress_info!(
      log,
      "## Estimating clock model with fixed rate {rate:.6e} (keep_root={keep_root})"
    );
  } else {
    progress_info!(log, "## Estimating clock model (keep_root={keep_root})");
  }

  progress_info!(log, "### Running backward regression");
  let state = clock_regression_backward(
    &tree.graph,
    &tree.inputs,
    outliers,
    options,
    &tree.branch_lengths,
    prev_clock_rate,
  )?;
  debug!("Backward regression completed");

  let (tree, state, reroot_result) = if keep_root {
    progress_info!(log, "### Keeping original root (--keep-root enabled)");
    (tree, state, None)
  } else {
    let reroot_params = clock_rate.map_or_else(
      || reroot_params.clone(),
      |rate| reroot_params.with_objective(RootObjective::FixedRate(rate)),
    );
    let reroot = RootSearch {
      outliers,
      options,
      optimization_params,
      reroot_params: &reroot_params,
      prev_clock_rate,
    };
    let (tree, state, reroot_result) = reroot_at_best_root(tree, state, &reroot, names, log)?;
    (tree, state, Some(reroot_result))
  };

  let points = clock_regression_points(
    &tree.graph,
    &tree.inputs,
    outliers,
    &tree.branch_lengths,
    prev_clock_rate,
  )?;

  progress_info!(log, "### Extracting clock model from root");
  let root_key = tree.graph.get_exactly_one_root()?.key();
  let root_clock_set = state.node(root_key).clone();

  let (regression, clock_model) = if let Some(rate) = clock_rate {
    progress_info!(log, "### Using fixed clock rate: {rate:.6e}");
    (None, Some(ClockModel::with_fixed_rate(&root_clock_set, rate)?))
  } else {
    progress_info!(log, "### Using estimated clock rate");
    let regression = ClockRegression::try_from(&root_clock_set)?;
    (Some(regression), None)
  };

  let rate = clock_model
    .as_ref()
    .map_or_else(|| regression.as_ref().unwrap().clock_rate(), |m| m.clock_rate());
  let intercept = clock_model
    .as_ref()
    .map_or_else(|| regression.as_ref().unwrap().intercept(), |m| m.intercept());
  progress_info!(log, "**Clock rate:** {rate:.6e}");
  progress_info!(log, "**Intercept:** {intercept:.4}");
  if let Some(reg) = &regression {
    progress_info!(log, "**R²:** {:.4}", reg.r_squared());
    progress_info!(log, "**χ²:** {:.4}", reg.chisq());
    progress_info!(log, "**Hessian:**\n{}", reg.hessian());
  }

  Ok((
    tree,
    ClockRerootResult {
      regression,
      clock_model,
      reroot_result,
      points,
    },
  ))
}

struct RootSearch<'a> {
  outliers: &'a BTreeSet<GraphNodeKey>,
  options: &'a ClockVarianceParams,
  optimization_params: &'a BranchPointOptimizationParams,
  reroot_params: &'a RerootParams,
  prev_clock_rate: Option<f64>,
}

fn reroot_at_best_root(
  tree: ClockTree,
  state: ClockState,
  search: &RootSearch<'_>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> Result<(ClockTree, ClockState, RerootResult), Report> {
  progress_info!(log, "### Running forward regression to find optimal root");
  let state = clock_regression_forward(
    &tree.graph,
    &tree.inputs,
    state,
    search.options,
    &tree.branch_lengths,
    search.prev_clock_rate,
  )?;
  debug!("Forward regression completed");

  progress_info!(log, "### Finding best root and rerooting tree");
  let best_root = search_root(&tree, &state, search, names, log)?;
  let (tree, state, best_root, stem_removal) = if root_moves(
    &tree.graph,
    best_root.edge,
    best_root.split,
    search.reroot_params.split_edge,
  )? {
    without_undated_stem(tree, state, best_root, search, names, log)?
  } else {
    (tree, state, best_root, None)
  };
  let (tree, state, reroot_result) = reroot_clock_tree(
    tree,
    state,
    best_root,
    search.options,
    search.reroot_params,
    stem_removal,
    names,
  )?;
  progress_info!(log, "Rerooted to {}", node_label(names, reroot_result.new_root_key));
  debug!("Rerooting completed");
  Ok((tree, state, reroot_result))
}

fn without_undated_stem(
  mut tree: ClockTree,
  state: ClockState,
  best_root: FindRootResult,
  search: &RootSearch<'_>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> Result<(ClockTree, ClockState, FindRootResult, Option<StemRemovalInfo>), Report> {
  let Some(stem) = remove_undated_stem(&mut tree)? else {
    return Ok((tree, state, best_root, None));
  };
  let state = clock_regression_backward(
    &tree.graph,
    &tree.inputs,
    search.outliers,
    search.options,
    &tree.branch_lengths,
    search.prev_clock_rate,
  )?;
  let state = clock_regression_forward(
    &tree.graph,
    &tree.inputs,
    state,
    search.options,
    &tree.branch_lengths,
    search.prev_clock_rate,
  )?;
  let best_root = search_root(&tree, &state, search, names, log)?;
  Ok((tree, state, best_root, Some(stem)))
}

fn search_root(
  tree: &ClockTree,
  state: &ClockState,
  search: &RootSearch<'_>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> Result<FindRootResult, Report> {
  select_root(
    &tree.graph,
    &tree.inputs,
    state,
    search.options,
    search.optimization_params,
    search.reroot_params,
    &tree.branch_lengths,
    names,
    log,
  )
}

#[derive(Clone, Debug, Serialize, Deserialize, deser::Serialize, deser::Deserialize)]
pub struct ClockRerootResult {
  regression: Option<ClockRegression>,
  clock_model: Option<ClockModel>,
  reroot_result: Option<RerootResult>,
  #[serde(skip)]
  #[deser(skip)]
  points: Vec<ClockRegressionPoint>,
}

pub(crate) struct ClockTree {
  pub graph: Graph,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  pub inputs: ClockInputs,
}

#[derive(Clone, Debug)]
pub(crate) struct ClockFit {
  pub model: ClockModel,
  pub points: Vec<ClockRegressionPoint>,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ClockRegressionPoint {
  pub key: GraphNodeKey,
  pub date: Option<f64>,
  pub div: f64,
  pub is_outlier: bool,
}

impl ClockRerootResult {
  pub(crate) fn into_clock_fit(mut self) -> Result<ClockFit, Report> {
    let points = mem::take(&mut self.points);
    Ok(ClockFit {
      model: self.into_clock_model()?,
      points,
    })
  }

  #[allow(
    clippy::expect_used,
    reason = "expect on a value an upstream invariant guarantees is present"
  )]
  fn into_clock_model(self) -> Result<ClockModel, Report> {
    if let Some(model) = self.clock_model {
      return Ok(model);
    }
    let regression = self
      .regression
      .expect("ClockRerootResult has neither regression nor clock_model");
    ClockModel::from_regression(&regression)
  }

  #[allow(
    clippy::expect_used,
    reason = "expect on a value an upstream invariant guarantees is present"
  )]
  pub(crate) fn into_clock_model_allow_negative(self, log: &dyn LogSink) -> ClockModel {
    if let Some(model) = self.clock_model {
      return model;
    }
    let regression = self
      .regression
      .expect("ClockRerootResult has neither regression nor clock_model");
    ClockModel::from_regression_allow_negative(&regression, log)
  }

  #[allow(
    clippy::expect_used,
    reason = "expect on a value an upstream invariant guarantees is present"
  )]
  pub(crate) fn regression(&self) -> &ClockRegression {
    self
      .regression
      .as_ref()
      .expect("regression() called on fixed-rate result")
  }

  pub(crate) fn reroot_result(&self) -> Option<&RerootResult> {
    self.reroot_result.as_ref()
  }
}

fn clock_regression_backward(
  graph: &Graph,
  inputs: &ClockInputs,
  outliers: &BTreeSet<GraphNodeKey>,
  options: &ClockVarianceParams,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  prev_clock_rate: Option<f64>,
) -> Result<ClockState, Report> {
  let mut state = ClockState::new(graph);
  state.map_backward(graph, |context| {
    clock_regression_backward_node(options, prev_clock_rate, branch_lengths, inputs, outliers, &context)
  })?;
  Ok(state)
}

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
fn clock_regression_backward_node(
  options: &ClockVarianceParams,
  prev_clock_rate: Option<f64>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  inputs: &ClockInputs,
  outliers: &BTreeSet<GraphNodeKey>,
  context: &GraphPassBackwardContext<'_, &ClockSet, ClockEdgeState, ClockSet, ClockEdgeState>,
) -> Result<GraphPassNodeOutput<ClockSet, ClockEdgeState>, Report> {
  let mut node = context.input.clone();
  let is_leaf = context.is_leaf;
  let date = inputs.likely_time(context.key);
  let q_to_parent = if is_leaf {
    if outliers.contains(&context.key) {
      ClockSet::outlier_contribution()
    } else {
      ClockSet::leaf_contribution(date)
    }
  } else {
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
    node = q_to_parent;
    None
  };

  Ok(GraphPassNodeOutput { node, parent_message })
}

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
fn clock_regression_forward(
  graph: &Graph,
  inputs: &ClockInputs,
  mut state: ClockState,
  options: &ClockVarianceParams,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  prev_clock_rate: Option<f64>,
) -> Result<ClockState, Report> {
  state.map_forward(graph, |context| {
    let mut node = context.input.clone();
    let parent_message = if let Some((edge_key, edge)) = context.parent_edge {
      let mut edge = edge.clone();
      let parent = context.parent.expect("Non-root node must have a parent");
      let mut q_to_child = parent.clone();
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
      node = q_dest;
      Some(edge)
    } else {
      None
    };
    Ok(GraphPassNodeOutput { node, parent_message })
  })?;
  Ok(state)
}

fn clock_regression_points(
  graph: &Graph,
  inputs: &ClockInputs,
  outliers: &BTreeSet<GraphNodeKey>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  prev_clock_rate: Option<f64>,
) -> Result<Vec<ClockRegressionPoint>, Report> {
  let divergences = root_to_node_divergences(graph, |edge_key| {
    let edge_input = inputs.edge(edge_key);
    edge_divergence(
      branch_lengths[&edge_key],
      edge_input.time_length,
      edge_input.gamma,
      prev_clock_rate,
    )
  })?;
  let mut points = vec![];
  graph.iter_depth_first_preorder_forward(|node| {
    if node.is_leaf {
      points.push(ClockRegressionPoint {
        key: node.key,
        date: inputs.likely_time(node.key),
        div: divergences[&node.key],
        is_outlier: outliers.contains(&node.key),
      });
    }
    Ok(())
  })?;
  Ok(points)
}

#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault, JsonSchema, deser::Serialize, deser::Deserialize)]
#[serde(default, deny_unknown_fields)]
#[deser(default, deny_unknown_fields)]
pub struct ClockVarianceParams {
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

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
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
