use crate::clock::date_constraints::DateConstraints;
use crate::coalescent::mean_field::CoalescentPrior;
use crate::progress::LogSink;
use crate::progress_warn;
use crate::timetree::inference::bad_branches::applied_date_constraint;
use crate::timetree::inference::coalescent_factor::apply_leaf_coalescent_factor;
use crate::timetree::inference::grid_request::{
  GridOperand, GridOperation, GridOrigin, GridRequest, OperandRole, WrapGridRequest,
};
use crate::timetree::inference::result::{
  BranchLikelihood, NodeBackward, PosteriorView, TimeDistribution, likely_time_of_node,
};
use crate::timetree::inference::runner::{EPS, GRID_POINTS, TRIM_EPS};
use crate::timetree::max_posterior_region::posterior_median;
use eyre::{Report, WrapErr};
use log::{Level, debug, log_enabled};
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_distribution::ForwardPosterior;
use treetime_distribution::MessageSpacing;
use treetime_distribution::NegLog;
use treetime_distribution::convolve_forward_across_edge;
use treetime_distribution::distribution_division;
use treetime_distribution::distribution_multiplication;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_graph::schedule::GraphSchedule;
use treetime_graph::schedule_scan::PerEdge;
use treetime_grid::MaxGridPoints;

pub(crate) fn propagate_distributions_forward<S: Send + Sync>(
  graph: &Graph,
  inputs: &ForwardPassInputs<'_>,
  backward: NodeMap<NodeBackward>,
  summarize: impl Fn(GraphNodeKey, PosteriorView<'_>) -> Result<S, Report> + Sync,
  log: &dyn LogSink,
) -> Result<NodeMap<S>, Report> {
  let schedule = GraphSchedule::new(graph)?;
  schedule.require_tree()?;
  let pass = ForwardPass {
    schedule: &schedule,
    inputs,
    summarize,
  };
  let outputs = schedule
    .scan_down(backward, |key, backward, parents: PerEdge<Posterior>| {
      pass.node(key, backward, parents.first())
    })?
    .into_map();

  for &key in schedule.preorder().iter().filter(|&&key| outputs[key].undated) {
    let name = &inputs.names[key];
    progress_warn!(
      log,
      "Timetree forward pass: node '{name}' has an empty time distribution; no date was assigned. \
       The messages meeting at this node leave no time with any probability: the dates below it \
       and the times the rest of the tree implies have disjoint support."
    );
  }
  let contradicted = outputs.values().filter(|output| output.contradicted).count();
  if contradicted > 0 {
    progress_warn!(
      log,
      "Timetree forward pass: {contradicted} node(s) carry a date that the rest of the tree gives \
       no probability at all, so their posterior came out empty and each kept the date it was \
       given, unrefined. The usual cause is a sequence whose divergence implies a date far from \
       the one it is stamped with, which the clock filter reports separately. Run with \
       `--verbosity=debug` to see which nodes and where the tree puts each of them."
    );
  }

  Ok(outputs.into_iter().map(|(key, output)| (key, output.summary)).collect())
}

struct ForwardPass<'a, F> {
  schedule: &'a GraphSchedule<'a>,
  inputs: &'a ForwardPassInputs<'a>,
  summarize: F,
}

impl<S, F> ForwardPass<'_, F>
where
  F: Fn(GraphNodeKey, PosteriorView<'_>) -> Result<S, Report>,
{
  fn node(
    &self,
    key: GraphNodeKey,
    backward: NodeBackward,
    parent: Option<&Posterior>,
  ) -> Result<(ForwardOutput<S>, Posterior), Report> {
    let ForwardPassInputs {
      constraints,
      coalescent,
      bad_branches,
      names,
      branches,
      max_points,
    } = *self.inputs;
    let date_constraint = applied_date_constraint(constraints, key, bad_branches[&key]);
    let NodeBackward {
      subtree,
      to_parent,
      subtree_origin,
      to_parent_origin,
    } = backward;
    let parent_edge = self
      .schedule
      .parents(key)
      .first()
      .map(|&(parent_key, edge_key)| ParentEdge {
        parent_key,
        branch: &branches[&edge_key],
        msg_to_parent: to_parent.as_deref(),
        msg_origin: to_parent_origin,
      });

    let is_leaf = self.schedule.children(key).is_empty();
    let refinement = refine_distribution_from_parent(ForwardStep {
      names,
      key,
      leaf_prior: coalescent.filter(|_| is_leaf),
      date_constraint: date_constraint.as_ref(),
      parent,
      edge: parent_edge,
      subtree: subtree.as_deref(),
      subtree_origin,
      max_points,
    })?;
    let normalized_subtree = || -> Result<_, Report> {
      subtree
        .as_deref()
        .map(|subtree| subtree.normalize().map(Arc::new))
        .transpose()
    };
    let (distribution, origin, contradicted) = match refinement {
      ForwardRefinement::Refined { distribution, origin } => (Some(Arc::new(distribution)), origin, false),
      ForwardRefinement::Unrefined => (normalized_subtree()?, subtree_origin, false),
      ForwardRefinement::ContradictedGivenDate => (normalized_subtree()?, subtree_origin, true),
    };

    let is_root = self.schedule.parents(key).is_empty();
    let likely_time = if is_root && !is_leaf {
      distribution
        .as_deref()
        .map(posterior_median)
        .transpose()
        .wrap_err_with(|| format!("When finding the median time of root {key}"))?
        .flatten()
    } else {
      likely_time_of_node(key, distribution.as_deref())?
    };
    let time = commit_node_time(date_constraint.as_ref(), parent, likely_time);
    let undated = time.is_none() && (!is_leaf || date_constraint.is_some());

    let view = PosteriorView {
      distribution: distribution.as_ref(),
      likely_time,
      time,
      contradicted,
      subtree_support: subtree.as_deref().and_then(Distribution::support_bounds),
    };
    let output = ForwardOutput {
      summary: (self.summarize)(key, view)?,
      undated,
      contradicted,
    };
    let posterior = Posterior {
      distribution,
      time,
      origin,
    };
    Ok((output, posterior))
  }
}

pub(crate) struct ForwardPassInputs<'a> {
  pub constraints: &'a DateConstraints,
  pub coalescent: Option<&'a CoalescentPrior>,
  pub bad_branches: &'a NodeMap<bool>,
  pub names: &'a NodeNames,
  pub branches: &'a EdgeMap<BranchLikelihood>,
  pub max_points: MaxGridPoints,
}

struct ForwardOutput<S> {
  summary: S,
  undated: bool,
  contradicted: bool,
}

fn refine_distribution_from_parent(step: ForwardStep<'_>) -> Result<ForwardRefinement, Report> {
  let ForwardStep {
    names,
    key,
    leaf_prior,
    date_constraint,
    parent,
    edge,
    subtree,
    subtree_origin,
    max_points,
  } = step;
  let (Some(parent), Some(edge)) = (parent, edge) else {
    return Ok(ForwardRefinement::Unrefined);
  };
  if has_exact_date(date_constraint) {
    return Ok(ForwardRefinement::Unrefined);
  }

  let (Some(parent_time_dist), Some(branch_dist)) = (&parent.distribution, &edge.branch.distribution) else {
    return Ok(ForwardRefinement::Unrefined);
  };
  let branch_operand = GridOperand::of(OperandRole::Branch(key), branch_dist, GridOrigin::branch(key));

  let Some(subtree_dist) = subtree else {
    let operands = [
      GridOperand::of(OperandRole::Posterior(edge.parent_key), parent_time_dist, parent.origin),
      branch_operand,
    ];
    let ForwardPosterior { posterior, spacing, .. } =
      convolve_forward_across_edge(parent_time_dist, branch_dist, EPS, GRID_POINTS, max_points)
        .and_then(|forward| {
          forward.resolve_posterior(TRIM_EPS, |message| {
            apply_leaf_coalescent_factor(leaf_prior, key, message.clone())?
              .normalize()
              .wrap_err_with(|| format!("When normalizing the time distribution of node {key}"))
          })
        })
        .wrap_grid_request(|| GridRequest::new(GridOperation::ForwardMessage(key), operands))?
        .into_untrimmed();
    log_refinement(names, key, parent_time_dist, &posterior, spacing, &posterior);
    return Ok(ForwardRefinement::Refined {
      distribution: posterior,
      origin: GridOperation::ForwardMessage(key).origin(operands.into_iter()),
    });
  };

  let parent_operand = GridOperand::of(OperandRole::Posterior(edge.parent_key), parent_time_dist, parent.origin);
  let (parent_except_subtree, parent_operand) =
    parent_without_child(key, parent_time_dist, parent_operand, edge, max_points)?;
  let operands = [parent_operand, branch_operand];
  let message_origin = GridOperation::ForwardMessage(key).origin(operands.into_iter());
  let subtree_operand = |subtree_dist| GridOperand::of(OperandRole::Subtree(key), subtree_dist, subtree_origin);
  let ForwardPosterior {
    message: dist_from_parent,
    posterior: combined,
    spacing,
  } = convolve_forward_across_edge(&parent_except_subtree, branch_dist, EPS, GRID_POINTS, max_points)
    .and_then(|forward| {
      forward.resolve_posterior(TRIM_EPS, |message| {
        let message_operand = GridOperand::of(OperandRole::ForwardMessage(key), message, message_origin);
        normalized_product(
          key,
          message,
          subtree_dist,
          [message_operand, subtree_operand(subtree_dist)],
          max_points,
        )
      })
    })
    .wrap_grid_request(|| GridRequest::new(GridOperation::ForwardMessage(key), operands))?
    .trim_posterior_to_mass(TRIM_EPS)
    .wrap_err_with(|| format!("When trimming the time distribution of node {key} to its mass"))?;
  log_refinement(names, key, parent_time_dist, &dist_from_parent, spacing, &combined);

  if combined.likely_time()?.is_none() && date_constraint.is_some() {
    log_kept_given_date(names, key, date_constraint, &dist_from_parent);
    return Ok(ForwardRefinement::ContradictedGivenDate);
  }
  let product_operands = [
    GridOperand::of(OperandRole::ForwardMessage(key), &dist_from_parent, message_origin),
    subtree_operand(subtree_dist),
  ];
  Ok(ForwardRefinement::Refined {
    distribution: combined,
    origin: GridOperation::Posterior(key).origin(product_operands.into_iter()),
  })
}

fn parent_without_child(
  key: GraphNodeKey,
  parent_time_dist: &Distribution<NegLog>,
  parent_operand: GridOperand,
  edge: ParentEdge<'_>,
  max_points: MaxGridPoints,
) -> Result<(Distribution<NegLog>, GridOperand), Report> {
  let Some(msg_to_parent) = edge.msg_to_parent else {
    return Ok((parent_time_dist.clone(), parent_operand));
  };
  let operands = [
    parent_operand,
    GridOperand::of(OperandRole::ChildMessage(key), msg_to_parent, edge.msg_origin),
  ];
  let divided = distribution_division(parent_time_dist, msg_to_parent, max_points)
    .wrap_grid_request(|| GridRequest::new(GridOperation::ParentWithoutChild(key), operands))?;
  let origin = GridOperation::ParentWithoutChild(key).origin(operands.into_iter());
  let operand = GridOperand::of(OperandRole::ParentWithoutChild(key), &divided, origin);
  Ok((divided, operand))
}

#[derive(Clone, Copy)]
struct ForwardStep<'a> {
  names: &'a NodeNames,
  key: GraphNodeKey,
  leaf_prior: Option<&'a CoalescentPrior>,
  date_constraint: Option<&'a Arc<Distribution<NegLog>>>,
  parent: Option<&'a Posterior>,
  edge: Option<ParentEdge<'a>>,
  subtree: Option<&'a Distribution<NegLog>>,
  subtree_origin: GridOrigin,
  max_points: MaxGridPoints,
}

fn normalized_product(
  key: GraphNodeKey,
  dist_from_parent: &Distribution<NegLog>,
  subtree_dist: &Distribution<NegLog>,
  operands: [GridOperand; 2],
  max_points: MaxGridPoints,
) -> Result<Distribution<NegLog>, Report> {
  distribution_multiplication(dist_from_parent, subtree_dist, max_points)
    .wrap_grid_request(|| GridRequest::new(GridOperation::Posterior(key), operands))?
    .normalize()
    .wrap_err_with(|| format!("When normalizing the time distribution of node {key}"))
}

#[derive(Clone, Copy)]
struct ParentEdge<'a> {
  parent_key: GraphNodeKey,
  branch: &'a BranchLikelihood,
  msg_to_parent: Option<&'a Distribution<NegLog>>,
  msg_origin: GridOrigin,
}

#[derive(Debug, Clone, PartialEq)]
enum ForwardRefinement {
  Unrefined,
  Refined {
    distribution: Distribution<NegLog>,
    origin: GridOrigin,
  },
  ContradictedGivenDate,
}

fn commit_node_time(
  date_constraint: Option<&Arc<Distribution<NegLog>>>,
  parent: Option<&Posterior>,
  likely_time: Option<f64>,
) -> Option<f64> {
  let parent_time = if has_exact_date(date_constraint) {
    None
  } else {
    parent.and_then(|parent| parent.time)
  };
  likely_time.map(|likely_time| committed_time(likely_time, parent_time))
}

#[derive(Clone)]
struct Posterior {
  distribution: TimeDistribution,
  time: Option<f64>,
  origin: GridOrigin,
}

fn has_exact_date(date_constraint: Option<&Arc<Distribution<NegLog>>>) -> bool {
  date_constraint.is_some_and(|dist| dist.is_point())
}

pub(super) fn committed_time(likely_time: f64, parent_time: Option<f64>) -> f64 {
  parent_time.map_or(likely_time, |parent_time| likely_time.max(parent_time))
}

fn log_refinement(
  names: &NodeNames,
  key: GraphNodeKey,
  parent: &Distribution<NegLog>,
  convolved: &Distribution<NegLog>,
  spacing: MessageSpacing,
  refined: &Distribution<NegLog>,
) {
  if !log_enabled!(Level::Debug) {
    return;
  }
  let name = &names[key];
  debug!(
    "Timetree forward pass: node '{name}': parent {} -> convolved {} at {spacing:?} spacing -> refined {}",
    describe_grid(parent),
    describe_grid(convolved),
    describe_grid(refined)
  );
}

fn log_kept_given_date(
  names: &NodeNames,
  key: GraphNodeKey,
  date_constraint: Option<&Arc<Distribution<NegLog>>>,
  dist_from_parent: &Distribution<NegLog>,
) {
  if !log_enabled!(Level::Debug) {
    return;
  }
  let name = &names[key];
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
