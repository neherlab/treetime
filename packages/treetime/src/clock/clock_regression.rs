use crate::branch_lengths::branch_length_or_zero;
use crate::clock::clock_model::{ClockModel, ClockRegression};
use crate::clock::clock_set::ClockSet;
use crate::clock::clock_state::{ClockEdgeState, ClockInputs, ClockState};
use crate::clock::divergence::root_to_node_divergences;
use crate::clock::find_best_root::params::RootObjective;
use crate::clock::reroot::{RerootParams, reroot_at_best_root};
use crate::progress::LogSink;
use crate::progress_info;
use eyre::Report;
use log::debug;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Debug;
use std::mem;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pass::{GraphPassBackwardContext, GraphPassNodeOutput};
use treetime_graph::reroot::RerootResult;

pub(crate) fn estimate_clock_model_with_reroot_policy(
  tree: ClockTree,
  params: &ClockFitParams<'_>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> Result<(ClockTree, ClockRerootResult), Report> {
  let ClockFitParams {
    outliers,
    clock_rate,
    keep_root,
    prev_clock_rate,
    ..
  } = *params;
  if let Some(rate) = clock_rate {
    progress_info!(
      log,
      "## Estimating clock model with fixed rate {rate:.6e} (keep_root={keep_root})"
    );
  } else {
    progress_info!(log, "## Estimating clock model (keep_root={keep_root})");
  }

  let mut tree = tree;
  let reroot_result = if keep_root {
    progress_info!(log, "### Keeping original root (--keep-root enabled)");
    None
  } else {
    progress_info!(log, "### Finding best root and rerooting tree");
    let reroot_params = clock_rate.map_or_else(
      || params.reroot_params.clone(),
      |rate| params.reroot_params.with_objective(RootObjective::FixedRate(rate)),
    );
    let ClockTree {
      graph,
      branch_lengths,
      inputs,
    } = &mut tree;
    Some(reroot_at_best_root(
      graph,
      branch_lengths,
      inputs,
      outliers,
      &reroot_params,
      names,
      log,
    )?)
  };

  progress_info!(log, "### Running backward regression");
  let state = clock_regression_backward(
    &tree.graph,
    &tree.inputs,
    outliers,
    &tree.branch_lengths,
    prev_clock_rate,
  )?;
  debug!("Backward regression completed");

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

  let estimate = if let Some(rate) = clock_rate {
    progress_info!(log, "### Using fixed clock rate: {rate:.6e}");
    ClockEstimate::Fixed(ClockModel::with_fixed_rate(&root_clock_set, rate)?)
  } else {
    progress_info!(log, "### Using estimated clock rate");
    ClockEstimate::Estimated(ClockRegression::try_from(&root_clock_set)?)
  };

  let (rate, intercept) = match &estimate {
    ClockEstimate::Estimated(regression) => (regression.clock_rate(), regression.intercept()),
    ClockEstimate::Fixed(model) => (model.clock_rate(), model.intercept()),
  };
  progress_info!(log, "**Clock rate:** {rate:.6e}");
  progress_info!(log, "**Intercept:** {intercept:.4}");
  if let ClockEstimate::Estimated(reg) = &estimate {
    progress_info!(log, "**R²:** {:.4}", reg.r_squared());
    progress_info!(log, "**χ²:** {:.4}", reg.chisq());
  }

  Ok((
    tree,
    ClockRerootResult {
      estimate,
      reroot_result,
      points,
      root_clock_set,
    },
  ))
}

#[derive(Clone, Copy)]
pub(crate) struct ClockFitParams<'a> {
  pub outliers: &'a BTreeSet<GraphNodeKey>,
  pub clock_rate: Option<f64>,
  pub keep_root: bool,
  pub reroot_params: &'a RerootParams,
  pub prev_clock_rate: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct ClockRerootResult {
  estimate: ClockEstimate,
  reroot_result: Option<RerootResult>,
  points: Vec<ClockRegressionPoint>,
  root_clock_set: ClockSet,
}

impl ClockRerootResult {
  pub(crate) fn into_clock_fit(mut self) -> Result<ClockFit, Report> {
    let points = mem::take(&mut self.points);
    Ok(ClockFit {
      model: self.into_clock_model()?,
      points,
    })
  }

  pub(crate) fn into_clock_fit_or(self, fallback: impl FnOnce() -> Result<f64, Report>) -> Result<ClockFit, Report> {
    match &self.estimate {
      ClockEstimate::Estimated(regression) if regression.clock_rate() <= 0.0 => {
        let model = ClockModel::with_fixed_rate(&self.root_clock_set, fallback()?)?;
        Ok(ClockFit {
          model,
          points: self.points,
        })
      },
      _ => self.into_clock_fit(),
    }
  }

  fn into_clock_model(self) -> Result<ClockModel, Report> {
    match self.estimate {
      ClockEstimate::Estimated(regression) => ClockModel::from_regression(&regression),
      ClockEstimate::Fixed(model) => Ok(model),
    }
  }

  pub(crate) fn into_clock_model_allow_negative(self, log: &dyn LogSink) -> ClockModel {
    match self.estimate {
      ClockEstimate::Estimated(regression) => ClockModel::from_regression_allow_negative(&regression, log),
      ClockEstimate::Fixed(model) => model,
    }
  }

  pub(crate) fn regression(&self) -> Option<&ClockRegression> {
    match &self.estimate {
      ClockEstimate::Estimated(regression) => Some(regression),
      ClockEstimate::Fixed(_) => None,
    }
  }

  pub(crate) fn reroot_result(&self) -> Option<&RerootResult> {
    self.reroot_result.as_ref()
  }
}

#[derive(Clone, Debug)]
enum ClockEstimate {
  Estimated(ClockRegression),
  Fixed(ClockModel),
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

pub(crate) fn clock_regression(
  graph: &Graph,
  inputs: &ClockInputs,
  outliers: &BTreeSet<GraphNodeKey>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  prev_clock_rate: Option<f64>,
) -> Result<ClockState, Report> {
  let state = clock_regression_backward(graph, inputs, outliers, branch_lengths, prev_clock_rate)?;
  clock_regression_forward(graph, inputs, state, branch_lengths, prev_clock_rate)
}

fn clock_regression_backward(
  graph: &Graph,
  inputs: &ClockInputs,
  outliers: &BTreeSet<GraphNodeKey>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  prev_clock_rate: Option<f64>,
) -> Result<ClockState, Report> {
  let mut state = ClockState::new(graph);
  state.map_backward(graph, |context| {
    clock_regression_backward_node(prev_clock_rate, branch_lengths, inputs, outliers, &context)
  })?;
  Ok(state)
}

fn clock_regression_backward_node(
  prev_clock_rate: Option<f64>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  inputs: &ClockInputs,
  outliers: &BTreeSet<GraphNodeKey>,
  context: &GraphPassBackwardContext<'_, &ClockSet, ClockEdgeState, ClockSet, ClockEdgeState>,
) -> Result<GraphPassNodeOutput<ClockSet, ClockEdgeState>, Report> {
  let q_to_parent = if !context.is_leaf {
    context
      .children
      .iter()
      .try_fold(ClockSet::default(), |mut total, child| {
        total += &child.edge_message()?.clock_from_child;
        Ok::<_, Report>(total)
      })?
  } else if outliers.contains(&context.key) {
    ClockSet::outlier_contribution()
  } else {
    ClockSet::leaf_contribution(inputs.likely_time(context.key))
  };

  let Some((edge_key, edge)) = context.parent_edge else {
    return Ok(GraphPassNodeOutput {
      node: q_to_parent,
      parent_message: None,
    });
  };
  let mut edge = edge.clone();
  edge.clock_from_child = q_to_parent.shift(edge_length(inputs, branch_lengths, edge_key, prev_clock_rate));
  edge.clock_to_parent = q_to_parent;
  Ok(GraphPassNodeOutput {
    node: context.input.clone(),
    parent_message: Some(edge),
  })
}

fn clock_regression_forward(
  graph: &Graph,
  inputs: &ClockInputs,
  mut state: ClockState,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  prev_clock_rate: Option<f64>,
) -> Result<ClockState, Report> {
  state.map_forward(graph, |context| {
    let (Some((edge_key, edge)), Some(parent)) = (context.parent_edge, context.parent) else {
      return Ok(GraphPassNodeOutput {
        node: context.input.clone(),
        parent_message: None,
      });
    };
    let mut edge = edge.clone();
    edge.clock_to_child = parent - &edge.clock_from_child;
    let node = &edge.clock_to_parent
      + edge
        .clock_to_child
        .shift(edge_length(inputs, branch_lengths, edge_key, prev_clock_rate));
    Ok(GraphPassNodeOutput {
      node,
      parent_message: Some(edge),
    })
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
    edge_length(inputs, branch_lengths, edge_key, prev_clock_rate)
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

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ClockRegressionPoint {
  pub key: GraphNodeKey,
  pub date: Option<f64>,
  pub div: f64,
  pub is_outlier: bool,
}

pub(crate) fn edge_length(
  inputs: &ClockInputs,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  edge_key: GraphEdgeKey,
  prev_clock_rate: Option<f64>,
) -> f64 {
  let edge_input = inputs.edge(edge_key);
  match (prev_clock_rate, edge_input.time_length) {
    (Some(rate), Some(time_length)) => time_length * rate * edge_input.gamma,
    _ => branch_length_or_zero(branch_lengths, edge_key),
  }
}
