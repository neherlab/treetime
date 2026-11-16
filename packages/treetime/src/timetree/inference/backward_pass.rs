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
use treetime_ops::FftPlans;

pub(crate) fn propagate_distributions_backward(
  graph: &Graph,
  constraints: &DateConstraints,
  coalescent_model: Option<&CoalescentPrior>,
  bad_branches: &NodeMap<bool>,
  branches: &EdgeMap<BranchLikelihood>,
  max_points: MaxGridPoints,
  fft_plans: &FftPlans,
) -> Result<NodeMap<NodeBackward>, Report> {
  let schedule = GraphSchedule::new(graph)?;
  schedule.require_tree()?;
  let from_parents = EdgeMap::new();
  let pass = BackwardPass {
    schedule: &schedule,
    constraints,
    coalescent_model,
    branches,
    max_points,
    fft_plans,
    from_parents: &from_parents,
    parent_messages: ParentMessages::default(),
  };
  pass.run(bad_branches)
}

pub(crate) fn propagate_distributions_backward_on_graph(
  graph: &Graph,
  inputs: &BackwardInputs<'_>,
  from_parents: &EdgeMap<Arc<Distribution<NegLog>>>,
  parent_messages: ParentMessages,
) -> Result<NodeMap<NodeBackward>, Report> {
  let schedule = GraphSchedule::new(graph)?;
  let pass = BackwardPass {
    schedule: &schedule,
    constraints: inputs.constraints,
    coalescent_model: inputs.coalescent_model,
    branches: inputs.branches,
    max_points: inputs.max_points,
    fft_plans: inputs.fft_plans,
    from_parents,
    parent_messages,
  };
  pass.run(inputs.bad_branches)
}

pub(crate) struct BackwardInputs<'a> {
  pub constraints: &'a DateConstraints,
  pub coalescent_model: Option<&'a CoalescentPrior>,
  pub bad_branches: &'a NodeMap<bool>,
  pub branches: &'a EdgeMap<BranchLikelihood>,
  pub max_points: MaxGridPoints,
  pub fft_plans: &'a FftPlans,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum ParentMessages {
  AllEdges,
  #[default]
  ChildrenOnly,
}

struct BackwardPass<'a> {
  schedule: &'a GraphSchedule<'a>,
  constraints: &'a DateConstraints,
  coalescent_model: Option<&'a CoalescentPrior>,
  branches: &'a EdgeMap<BranchLikelihood>,
  max_points: MaxGridPoints,
  fft_plans: &'a FftPlans,
  from_parents: &'a EdgeMap<Arc<Distribution<NegLog>>>,
  parent_messages: ParentMessages,
}

type ParentSlotMessages = Arc<PerEdge<TimeMessage>>;

impl BackwardPass<'_> {
  fn run(&self, bad_branches: &NodeMap<bool>) -> Result<NodeMap<NodeBackward>, Report> {
    let outputs = self.schedule.scan_up(
      bad_branches.iter().map(|(key, &bad_branch)| (key, bad_branch)),
      |key, bad_branch, messages| self.node(key, bad_branch, &messages),
    )?;
    Ok(outputs.into_map())
  }

  fn node(
    &self,
    key: GraphNodeKey,
    bad_branch: bool,
    messages: &[ParentSlotMessages],
  ) -> Result<(NodeBackward, ParentSlotMessages), Report> {
    let BackwardPass {
      schedule,
      constraints,
      coalescent_model,
      branches,
      max_points,
      fft_plans,
      ..
    } = *self;
    let date_constraint = applied_date_constraint(constraints, key, bad_branch);
    let children = schedule.children(key);
    let messages = izip!(children, messages)
      .filter_map(|(&(child, edge_key), messages)| {
        let message = messages.get(schedule.parent_slot(edge_key))?;
        let distribution = message.distribution.clone()?;
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

    let parents = schedule.parents(key);
    let to_parents = parents
      .iter()
      .enumerate()
      .map(|(slot, &(_, edge_key))| {
        let base = self.base_for_parent(key, slot, subtree.as_ref())?;
        send_backward_message(
          key,
          bad_branch,
          base.as_deref(),
          subtree_origin,
          &branches[&edge_key],
          max_points,
          fft_plans,
        )
      })
      .collect::<Result<PerEdge<_>, Report>>()?;
    let output = NodeBackward {
      subtree,
      to_parents: to_parents.iter().map(|message| message.distribution.clone()).collect(),
      subtree_origin,
      to_parent_origins: to_parents.iter().map(|message| message.origin).collect(),
    };
    Ok((output, Arc::new(to_parents)))
  }

  fn base_for_parent(
    &self,
    key: GraphNodeKey,
    slot: usize,
    subtree: Option<&Arc<Distribution<NegLog>>>,
  ) -> Result<Option<Arc<Distribution<NegLog>>>, Report> {
    let parents = self.schedule.parents(key);
    if parents.len() < 2 || self.parent_messages == ParentMessages::ChildrenOnly {
      return Ok(subtree.cloned());
    }
    let others = parents
      .iter()
      .enumerate()
      .filter(|&(other, _)| other != slot)
      .filter_map(|(_, (_, edge_key))| self.from_parents.get(edge_key))
      .collect_vec();
    if others.is_empty() {
      return Ok(subtree.cloned());
    }
    let factors = subtree
      .into_iter()
      .chain(others)
      .map(Arc::as_ref)
      .filter(|factor| !matches!(factor, Distribution::Empty))
      .collect_vec();
    let product = distribution_product(&factors, self.max_points)
      .wrap_err_with(|| format!("When combining the messages to the parents of node {key}"))?;
    if matches!(product, Distribution::Empty) {
      return Ok(None);
    }
    Ok(Some(Arc::new(product.repeak()?)))
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
  fft_plans: &FftPlans,
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
  let message = convolve_across_edge(
    distribution,
    &negated_branch,
    Side::Left,
    EPS,
    GRID_POINTS,
    max_points,
    fft_plans,
  )
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
