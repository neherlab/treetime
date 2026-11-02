use crate::clock::date_constraints::DateConstraints;
use crate::coalescent::mean_field::CoalescentPrior;
use crate::timetree::inference::bad_branches::applied_date_constraint;
use crate::timetree::inference::coalescent_factor::{apply_coalescent_prior, apply_leaf_coalescent_factor};
use crate::timetree::inference::grid_request::{
  GridOperand, GridOperation, GridOrigin, GridRequest, OperandRole, WrapGridRequest,
};
use crate::timetree::inference::result::{BranchLikelihood, NodeBackward, TimeDistribution};
use crate::timetree::inference::runner::{EPS, GRID_POINTS};
use eyre::{Report, WrapErr};
use itertools::{Itertools, izip};
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_distribution::NegLog;
use treetime_distribution::convolve_across_edge;
use treetime_distribution::distribution_multiplication;
use treetime_distribution::distribution_product;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::node::GraphNodeKey;
use treetime_graph::schedule::GraphSchedule;
use treetime_graph::schedule_scan::PerEdge;
use treetime_grid::{MaxGridPoints, Side};

pub(crate) fn propagate_distributions_backward(
  graph: &Graph,
  constraints: &DateConstraints,
  coalescent_model: Option<&CoalescentPrior>,
  bad_branches: &NodeMap<bool>,
  branches: &EdgeMap<BranchLikelihood>,
  max_points: MaxGridPoints,
) -> Result<NodeMap<NodeBackward>, Report> {
  let schedule = GraphSchedule::new(graph)?;
  schedule.require_tree()?;
  let pass = BackwardPass {
    schedule: &schedule,
    constraints,
    coalescent_model,
    branches,
    max_points,
  };
  let outputs = schedule.scan_up(
    bad_branches.iter().map(|(key, &bad_branch)| (key, bad_branch)),
    |key, bad_branch, messages| pass.node(key, bad_branch, messages),
  )?;
  Ok(outputs.into_map())
}

struct BackwardPass<'a> {
  schedule: &'a GraphSchedule<'a>,
  constraints: &'a DateConstraints,
  coalescent_model: Option<&'a CoalescentPrior>,
  branches: &'a EdgeMap<BranchLikelihood>,
  max_points: MaxGridPoints,
}

impl BackwardPass<'_> {
  fn node(
    &self,
    key: GraphNodeKey,
    bad_branch: bool,
    messages: PerEdge<TimeMessage>,
  ) -> Result<(NodeBackward, TimeMessage), Report> {
    let BackwardPass {
      schedule,
      constraints,
      coalescent_model,
      branches,
      max_points,
    } = *self;
    let date_constraint = applied_date_constraint(constraints, key, bad_branch);
    let children = schedule.children(key);
    let messages = izip!(children, messages)
      .filter_map(|(&(child, _), message)| {
        let distribution = message.distribution?;
        let operand = GridOperand::of(OperandRole::ChildMessage(child), &distribution, message.origin);
        Some((distribution, operand))
      })
      .collect_vec();
    let (distribution, origin) = combine_child_messages(key, &messages, max_points)?;
    let distribution = apply_coalescent_prior(coalescent_model, key, children.len(), distribution)?;
    let (distribution, subtree_origin) =
      apply_date_constraint(key, date_constraint.as_ref(), distribution, origin, max_points)?;
    let distribution =
      apply_leaf_coalescent_factor(coalescent_model.filter(|_| children.is_empty()), key, distribution)?;

    let subtree = if matches!(distribution, Distribution::Empty) {
      None
    } else {
      let distribution = distribution
        .repeak()
        .wrap_err_with(|| format!("When re-peaking the time distribution of node {key}"))?;
      Some(Arc::new(distribution))
    };

    let to_parent = schedule
      .parents(key)
      .first()
      .map(|&(_, edge_key)| {
        send_backward_message(
          key,
          bad_branch,
          subtree.as_deref(),
          subtree_origin,
          &branches[&edge_key],
          max_points,
        )
      })
      .transpose()?
      .unwrap_or_default();
    let output = NodeBackward {
      subtree,
      to_parent: to_parent.distribution.clone(),
      subtree_origin,
      to_parent_origin: to_parent.origin,
    };
    Ok((output, to_parent))
  }
}

fn combine_child_messages(
  key: GraphNodeKey,
  messages: &[(Arc<Distribution<NegLog>>, GridOperand)],
  max_points: MaxGridPoints,
) -> Result<(Distribution<NegLog>, GridOrigin), Report> {
  let factors = messages
    .iter()
    .filter(|(message, _)| !matches!(message.as_ref(), Distribution::Empty))
    .collect_vec();
  if factors.is_empty() {
    return Ok((Distribution::Empty, GridOrigin::default()));
  }
  let operands = factors.iter().map(|(_, operand)| *operand);
  let distributions = factors.iter().map(|(message, _)| message.as_ref()).collect_vec();
  let product = distribution_product(&distributions, max_points)
    .wrap_grid_request(|| GridRequest::new(GridOperation::ChildMessages(key), operands.clone()))?;
  Ok((product, GridOperation::ChildMessages(key).origin(operands)))
}

fn apply_date_constraint(
  key: GraphNodeKey,
  date_constraint: Option<&Arc<Distribution<NegLog>>>,
  distribution: Distribution<NegLog>,
  origin: GridOrigin,
  max_points: MaxGridPoints,
) -> Result<(Distribution<NegLog>, GridOrigin), Report> {
  let Some(constraint) = date_constraint else {
    return Ok((distribution, origin));
  };
  if matches!(distribution, Distribution::Empty) {
    return Ok((constraint.as_ref().clone(), GridOrigin::date_constraint(key)));
  }
  let operands = [
    GridOperand::of(OperandRole::ChildProduct(key), &distribution, origin),
    GridOperand::of(
      OperandRole::DateConstraint(key),
      constraint,
      GridOrigin::date_constraint(key),
    ),
  ];
  let product = distribution_multiplication(&distribution, constraint, max_points)
    .wrap_grid_request(|| GridRequest::new(GridOperation::DateConstraint(key), operands))?;
  Ok((product, GridOperation::DateConstraint(key).origin(operands.into_iter())))
}

fn send_backward_message(
  key: GraphNodeKey,
  bad_branch: bool,
  subtree: Option<&Distribution<NegLog>>,
  subtree_origin: GridOrigin,
  branch: &BranchLikelihood,
  max_points: MaxGridPoints,
) -> Result<TimeMessage, Report> {
  if bad_branch {
    return Ok(TimeMessage::default());
  }
  let (Some(distribution), Some(branch_length_distribution)) = (subtree, &branch.distribution) else {
    return Ok(TimeMessage::default());
  };

  let negated_branch = branch_length_distribution.negate()?;
  let operands = [
    GridOperand::of(OperandRole::Subtree(key), distribution, subtree_origin),
    GridOperand::of(OperandRole::Branch(key), &negated_branch, GridOrigin::branch(key)),
  ];
  let message = convolve_across_edge(distribution, &negated_branch, Side::Left, EPS, GRID_POINTS, max_points)
    .wrap_grid_request(|| GridRequest::new(GridOperation::BackwardMessage(key), operands))?;

  Ok(TimeMessage {
    distribution: Some(Arc::new(message)),
    origin: GridOperation::BackwardMessage(key).origin(operands.into_iter()),
  })
}

#[derive(Clone, Debug, Default)]
struct TimeMessage {
  distribution: TimeDistribution,
  origin: GridOrigin,
}
