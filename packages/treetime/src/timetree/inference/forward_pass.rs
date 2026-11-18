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
use itertools::{Itertools, izip};
use log::{Level, debug, log_enabled};
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_distribution::ForwardPosterior;
use treetime_distribution::MessageSpacing;
use treetime_distribution::NegLog;
use treetime_distribution::convolve_forward_across_edge;
use treetime_distribution::distribution_division;
use treetime_distribution::distribution_multiplication;
use treetime_distribution::distribution_product;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_graph::schedule::GraphSchedule;
use treetime_graph::schedule_scan::PerEdge;
use treetime_grid::MaxGridPoints;
use treetime_ops::FftPlans;

#[cfg(test)]
pub(crate) fn propagate_distributions_forward<S: Send + Sync>(
  graph: &Graph,
  inputs: &ForwardPassInputs<'_>,
  backward: NodeMap<NodeBackward>,
  summarize: impl Fn(GraphNodeKey, PosteriorView<'_>) -> Result<S, Report> + Sync,
  log: &dyn LogSink,
) -> Result<NodeMap<S>, Report> {
  Ok(propagate_distributions_forward_on_graph(graph, inputs, backward, summarize, log)?.0)
}

pub(crate) fn propagate_distributions_forward_on_graph<S: Send + Sync>(
  graph: &Graph,
  inputs: &ForwardPassInputs<'_>,
  backward: NodeMap<NodeBackward>,
  summarize: impl Fn(GraphNodeKey, PosteriorView<'_>) -> Result<S, Report> + Sync,
  log: &dyn LogSink,
) -> Result<(NodeMap<S>, EdgeMap<Arc<Distribution<NegLog>>>), Report> {
  let schedule = GraphSchedule::new(graph)?;
  let pass = ForwardPass {
    schedule: &schedule,
    inputs,
    summarize,
  };
  let outputs = schedule
    .scan_down(backward, |key, backward, parents: PerEdge<Posterior>| {
      pass.node(key, backward, &parents)
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
  let disjoint = schedule
    .preorder()
    .iter()
    .filter(|&&key| outputs[key].parents_disjoint)
    .map(|&key| &inputs.names[key])
    .collect_vec();
  if !disjoint.is_empty() {
    progress_warn!(
      log,
      "Timetree forward pass: {} node(s) with several parents received parent messages that leave no \
       time with any probability together with their subtree, so each was dated from its subtree alone. \
       The parents of such a node place it at incompatible times, for example a reassortant whose \
       segment trees disagree on its age: {}",
      disjoint.len(),
      disjoint.iter().join(", ")
    );
  }

  let mut summaries = NodeMap::new();
  let mut from_parents = EdgeMap::new();
  for (key, output) in outputs {
    from_parents.extend(output.from_parents);
    summaries.insert(key, output.summary);
  }
  Ok((summaries, from_parents))
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
    parents: &[Posterior],
  ) -> Result<(ForwardOutput<S>, Posterior), Report> {
    if self.schedule.parents(key).len() > 1 {
      return self.node_with_parents(key, backward, parents);
    }
    self.node_with_parent(key, backward, parents.first())
  }

  fn node_with_parent(
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
      fft_plans,
    } = *self.inputs;
    let date_constraint = applied_date_constraint(constraints, key, bad_branches[&key]);
    let NodeBackward {
      subtree,
      to_parents,
      subtree_origin,
      to_parent_origins,
    } = backward;
    let to_parent = to_parents.into_iter().next().flatten();
    let to_parent_origin = to_parent_origins.into_iter().next().unwrap_or_default();
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
      fft_plans,
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
      parents_disjoint: false,
      from_parents: PerEdge::new(),
    };
    let posterior = Posterior {
      distribution,
      time,
      origin,
    };
    Ok((output, posterior))
  }

  fn node_with_parents(
    &self,
    key: GraphNodeKey,
    backward: NodeBackward,
    parents: &[Posterior],
  ) -> Result<(ForwardOutput<S>, Posterior), Report> {
    let ForwardPassInputs {
      constraints,
      bad_branches,
      branches,
      max_points,
      fft_plans,
      ..
    } = *self.inputs;
    let date_constraint = applied_date_constraint(constraints, key, bad_branches[&key]);
    let NodeBackward {
      subtree,
      to_parents,
      subtree_origin,
      to_parent_origins,
    } = backward;
    let mut from_parents = PerEdge::new();
    for (&(parent_key, edge_key), parent, to_parent, &origin) in
      izip!(self.schedule.parents(key), parents, &to_parents, &to_parent_origins)
    {
      let edge = ParentEdge {
        parent_key,
        branch: &branches[&edge_key],
        msg_to_parent: to_parent.as_deref(),
        msg_origin: origin,
      };
      if let Some(message) = message_from_parent(key, parent, edge, max_points, fft_plans)? {
        from_parents.push((edge_key, Arc::new(message)));
      }
    }
    let exact = has_exact_date(date_constraint.as_ref());
    let factors = subtree
      .as_deref()
      .into_iter()
      .chain(
        from_parents
          .iter()
          .filter(|_| !exact)
          .map(|(_, message)| message.as_ref()),
      )
      .collect_vec();
    let product = if factors.is_empty() {
      Distribution::Empty
    } else {
      distribution_product(&factors, max_points)
        .wrap_err_with(|| format!("When combining the messages at node {key}"))?
        .normalize()
        .wrap_err_with(|| format!("When normalizing the time distribution of node {key}"))?
    };
    let empty = matches!(product, Distribution::Empty);
    let contradicted = empty && date_constraint.is_some();
    let parents_disjoint = empty && !contradicted && !exact && !from_parents.is_empty();
    let distribution = if empty {
      subtree
        .as_deref()
        .map(Distribution::normalize)
        .transpose()?
        .map(Arc::new)
    } else {
      Some(Arc::new(product))
    };
    let likely_time = likely_time_of_node(key, distribution.as_deref())?;
    let parent_time = if exact {
      None
    } else {
      parents.iter().filter_map(|parent| parent.time).reduce(f64::max)
    };
    let time = likely_time.map(|likely_time| committed_time(likely_time, parent_time));
    let is_leaf = self.schedule.children(key).is_empty();
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
      parents_disjoint,
      from_parents,
    };
    let posterior = Posterior {
      distribution,
      time,
      origin: subtree_origin,
    };
    Ok((output, posterior))
  }
}

fn message_from_parent(
  key: GraphNodeKey,
  parent: &Posterior,
  edge: ParentEdge<'_>,
  max_points: MaxGridPoints,
  fft_plans: &FftPlans,
) -> Result<Option<Distribution<NegLog>>, Report> {
  let (Some(parent_time_dist), Some(branch_dist)) = (&parent.distribution, &edge.branch.distribution) else {
    return Ok(None);
  };
  let parent_operand = GridOperand::of(OperandRole::Posterior(edge.parent_key), parent_time_dist, parent.origin);
  let (parent_except_child, parent_operand) =
    parent_without_child(key, parent_time_dist, parent_operand, edge, max_points)?;
  let branch_operand = GridOperand::of(OperandRole::Branch(key), branch_dist, GridOrigin::branch(key));
  let operands = [parent_operand, branch_operand];
  let ForwardPosterior { message, .. } = convolve_forward_across_edge(
    &parent_except_child,
    branch_dist,
    EPS,
    GRID_POINTS,
    max_points,
    fft_plans,
  )
  .and_then(|forward| forward.resolve_posterior(TRIM_EPS, Distribution::normalize))
  .wrap_grid_request(|| GridRequest::new(GridOperation::ForwardMessage(key), operands))?
  .into_untrimmed();
  let message = message
    .normalize()
    .wrap_err_with(|| format!("When normalizing the message from parent to node {key}"))?;
  Ok((!matches!(message, Distribution::Empty)).then_some(message))
}

pub(crate) struct ForwardPassInputs<'a> {
  pub constraints: &'a DateConstraints,
  pub coalescent: Option<&'a CoalescentPrior>,
  pub bad_branches: &'a NodeMap<bool>,
  pub names: &'a NodeNames,
  pub branches: &'a EdgeMap<BranchLikelihood>,
  pub max_points: MaxGridPoints,
  pub fft_plans: &'a FftPlans,
}

struct ForwardOutput<S> {
  summary: S,
  undated: bool,
  contradicted: bool,
  parents_disjoint: bool,
  from_parents: PerEdge<(GraphEdgeKey, Arc<Distribution<NegLog>>)>,
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
    fft_plans,
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
      convolve_forward_across_edge(parent_time_dist, branch_dist, EPS, GRID_POINTS, max_points, fft_plans)
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
  } = convolve_forward_across_edge(
    &parent_except_subtree,
    branch_dist,
    EPS,
    GRID_POINTS,
    max_points,
    fft_plans,
  )
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
  fft_plans: &'a FftPlans,
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
