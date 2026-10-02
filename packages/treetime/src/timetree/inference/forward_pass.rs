use crate::clock::date_constraints::DateConstraints;
use crate::node_label::node_label;
use crate::progress::ProgressSink;
use crate::progress_warn;
use crate::timetree::inference::runner::{EPS, GRID_POINTS};
use crate::timetree::inference::time_inference::{BranchLikelihood, NodePosterior, TimeBackward, TimeMessage};
use eyre::{Report, WrapErr};
use log::{Level, debug, log_enabled};
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_distribution::NegLog;
use treetime_distribution::convolve_across_edge;
use treetime_distribution::distribution_division;
use treetime_distribution::distribution_multiplication;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pass::{GraphMapOutputs, GraphPass, GraphPassForwardContext, GraphPassNodeOutput};
use treetime_grid::Side;
use treetime_utils::make_internal_report;

pub(crate) fn propagate_distributions_forward(
  graph: &Graph,
  constraints: &DateConstraints,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  branches: &BTreeMap<GraphEdgeKey, BranchLikelihood>,
  backward: &TimeBackward,
  progress: &dyn ProgressSink,
) -> Result<BTreeMap<GraphNodeKey, NodePosterior>, Report> {
  let pass = GraphPass::new(graph)?;
  let GraphMapOutputs { nodes: posterior, .. } = pass.map_forward(
    &backward.subtree,
    branches,
    |key| Err(make_internal_report!("Backward pass output is missing node {key}")),
    |context| propagate_distributions_forward_node(constraints, names, &backward.messages, &context, progress),
  )?;

  let contradicted = posterior.values().filter(|node| node.contradicted).count();
  if contradicted > 0 {
    progress_warn!(
      progress,
      "Timetree forward pass: {contradicted} node(s) carry a date that the rest of the tree gives \
       no probability at all, so their posterior came out empty and each kept the date it was \
       given, unrefined. The usual cause is a sequence whose divergence implies a date far from \
       the one it is stamped with, which the clock filter reports separately. Run with \
       `--verbosity=debug` to see which nodes and where the tree puts each of them."
    );
  }

  Ok(posterior)
}

fn propagate_distributions_forward_node(
  constraints: &DateConstraints,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  messages: &BTreeMap<GraphEdgeKey, TimeMessage>,
  context: &GraphPassForwardContext<'_, TimeMessage, BranchLikelihood, NodePosterior>,
  progress: &dyn ProgressSink,
) -> Result<GraphPassNodeOutput<NodePosterior, ()>, Report> {
  let date_constraint = constraints.date_constraint(context.key);
  let subtree = context.input;
  let parent_edge = context.parent_edge.map(|(edge_key, branch)| ParentEdge {
    branch,
    msg_to_parent: messages[&edge_key].as_deref(),
  });

  let refinement = refine_distribution_from_parent(
    names,
    context.key,
    date_constraint.as_ref(),
    context.parent,
    parent_edge,
    subtree.as_deref(),
  )?;
  let (distribution, contradicted) = match refinement {
    ForwardRefinement::Refined(distribution) => (Some(Arc::new(distribution)), false),
    ForwardRefinement::Unrefined => (subtree.clone(), false),
    ForwardRefinement::ContradictedGivenDate => (subtree.clone(), true),
  };

  let time = commit_node_time(
    names,
    context.key,
    date_constraint.as_ref(),
    context.parent,
    context.is_leaf,
    distribution.as_deref(),
    progress,
  )?;

  Ok(GraphPassNodeOutput {
    node: NodePosterior {
      distribution,
      time,
      contradicted,
    },
    parent_message: None,
  })
}

#[derive(Clone, Copy)]
struct ParentEdge<'a> {
  branch: &'a BranchLikelihood,
  msg_to_parent: Option<&'a Distribution<NegLog>>,
}

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
fn refine_distribution_from_parent(
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  key: GraphNodeKey,
  date_constraint: Option<&Arc<Distribution<NegLog>>>,
  parent: Option<&NodePosterior>,
  edge: Option<ParentEdge<'_>>,
  subtree: Option<&Distribution<NegLog>>,
) -> Result<ForwardRefinement, Report> {
  let Some(parent) = parent else {
    return Ok(ForwardRefinement::Unrefined);
  };
  if has_exact_date(date_constraint) {
    return Ok(ForwardRefinement::Unrefined);
  }

  let edge = edge.expect("Non-root indexed node must own its parent edge");

  let (Some(parent_time_dist), Some(branch_dist)) = (&parent.distribution, &edge.branch.distribution) else {
    return Ok(ForwardRefinement::Unrefined);
  };

  let Some(subtree_dist) = subtree else {
    let dist_from_parent = convolve_across_edge(parent_time_dist, branch_dist, Side::Right, EPS, GRID_POINTS)?;
    log_refinement(names, key, parent_time_dist, &dist_from_parent);
    return Ok(ForwardRefinement::Refined(dist_from_parent));
  };

  let parent_except_subtree = match edge.msg_to_parent {
    Some(msg_to_parent) => distribution_division(parent_time_dist, msg_to_parent)?,
    None => parent_time_dist.as_ref().clone(),
  };
  let dist_from_parent = convolve_across_edge(&parent_except_subtree, branch_dist, Side::Right, EPS, GRID_POINTS)?;

  let combined = distribution_multiplication(&dist_from_parent, subtree_dist)?
    .normalize()
    .wrap_err_with(|| format!("When normalizing the time distribution of node {key}"))?;
  log_refinement(names, key, parent_time_dist, &combined);

  if combined.likely_time()?.is_none() && date_constraint.is_some() {
    log_kept_given_date(names, key, date_constraint, &dist_from_parent);
    return Ok(ForwardRefinement::ContradictedGivenDate);
  }
  Ok(ForwardRefinement::Refined(combined))
}

#[derive(Debug, Clone, PartialEq)]
enum ForwardRefinement {
  Unrefined,
  Refined(Distribution<NegLog>),
  ContradictedGivenDate,
}

fn commit_node_time(
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  key: GraphNodeKey,
  date_constraint: Option<&Arc<Distribution<NegLog>>>,
  parent: Option<&NodePosterior>,
  is_leaf: bool,
  distribution: Option<&Distribution<NegLog>>,
  progress: &dyn ProgressSink,
) -> Result<Option<f64>, Report> {
  let parent_time = if has_exact_date(date_constraint) {
    None
  } else {
    parent.and_then(|parent| parent.time)
  };

  let time = committed_time(distribution, parent_time)?;
  let is_dateable = !is_leaf || date_constraint.is_some();
  if time.is_none() && is_dateable {
    let name = node_label(names, key);
    progress_warn!(
      progress,
      "Timetree forward pass: node '{name}' has an empty time distribution; no date was assigned. \
       The messages meeting at this node leave no time with any probability: the dates below it \
       and the times the rest of the tree implies have disjoint support."
    );
  }
  Ok(time)
}

fn has_exact_date(date_constraint: Option<&Arc<Distribution<NegLog>>>) -> bool {
  date_constraint.is_some_and(|dist| dist.is_point())
}

pub(super) fn committed_time(
  distribution: Option<&Distribution<NegLog>>,
  parent_time: Option<f64>,
) -> Result<Option<f64>, Report> {
  let Some(distribution) = distribution else {
    return Ok(None);
  };
  let Some(time) = distribution.likely_time()? else {
    return Ok(None);
  };
  Ok(Some(parent_time.map_or(time, |parent_time| time.max(parent_time))))
}

fn log_refinement(
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  key: GraphNodeKey,
  parent: &Distribution<NegLog>,
  refined: &Distribution<NegLog>,
) {
  if !log_enabled!(Level::Debug) {
    return;
  }
  let name = node_label(names, key);
  debug!(
    "Timetree forward pass: node '{name}': parent {} -> refined {}",
    describe_grid(parent),
    describe_grid(refined)
  );
}

fn log_kept_given_date(
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  key: GraphNodeKey,
  date_constraint: Option<&Arc<Distribution<NegLog>>>,
  dist_from_parent: &Distribution<NegLog>,
) {
  if !log_enabled!(Level::Debug) {
    return;
  }
  let name = node_label(names, key);
  let given = date_constraint.map_or_else(|| "none".to_owned(), |constraint| describe_grid(constraint.as_ref()));
  debug!(
    "Timetree forward pass: node '{name}' keeps the date it was given, {given}: the rest of the \
     tree puts it at {}, which leaves no probability on that date",
    describe_grid(dist_from_parent)
  );
}

fn describe_grid(dist: &Distribution<NegLog>) -> String {
  let n_points = dist.t().len();
  match dist.time_bounds() {
    Some((t_min, t_max)) => format!("n={n_points}, [{t_min:.6}, {t_max:.6}]"),
    None => "empty".to_owned(),
  }
}
