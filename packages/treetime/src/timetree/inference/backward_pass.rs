use crate::clock::date_constraints::DateConstraints;
use crate::coalescent::mean_field::CoalescentPrior;
use crate::timetree::inference::result::{BranchLikelihood, NodeBackward, TimeDistribution};
use crate::timetree::inference::runner::{EPS, GRID_POINTS};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use ndarray::Array1;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::NegLog;
use treetime_distribution::convolve_across_edge;
use treetime_distribution::distribution_multiplication;
use treetime_distribution::distribution_multiply_by_fn;
use treetime_distribution::distribution_product;
use treetime_distribution::{BoundaryBehavior, Distribution, DistributionFunction};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::tree_view::TreeView;
use treetime_grid::{MaxGridPoints, Side};
use treetime_utils::cast::{ToF64Exact, f64_to_usize_saturating};

pub(crate) fn propagate_distributions_backward(
  graph: &Graph,
  constraints: &DateConstraints,
  coalescent_model: Option<&CoalescentPrior>,
  bad_branches: &BTreeMap<GraphNodeKey, bool>,
  branches: &BTreeMap<GraphEdgeKey, BranchLikelihood>,
  max_points: MaxGridPoints,
) -> Result<BTreeMap<GraphNodeKey, NodeBackward>, Report> {
  let tree = TreeView::new(graph)?;
  let pass = BackwardPass {
    tree: &tree,
    constraints,
    coalescent_model,
    branches,
    max_points,
  };
  let outputs = tree.scan_up(
    bad_branches.iter().map(|(&key, &bad_branch)| (key, bad_branch)),
    |key, bad_branch, messages| pass.node(key, bad_branch, messages),
  )?;
  Ok(outputs.into_map())
}

struct BackwardPass<'a> {
  tree: &'a TreeView<'a>,
  constraints: &'a DateConstraints,
  coalescent_model: Option<&'a CoalescentPrior>,
  branches: &'a BTreeMap<GraphEdgeKey, BranchLikelihood>,
  max_points: MaxGridPoints,
}

impl BackwardPass<'_> {
  fn node(
    &self,
    key: GraphNodeKey,
    bad_branch: bool,
    messages: Vec<TimeDistribution>,
  ) -> Result<(NodeBackward, TimeDistribution), Report> {
    let BackwardPass {
      tree,
      constraints,
      coalescent_model,
      branches,
      max_points,
    } = *self;
    let date_constraint = constraints.date_constraint(key);
    let n_children = tree.children(key).len();
    let messages = messages.into_iter().flatten().collect_vec();
    let distribution = combine_child_messages(&messages, max_points)?;
    let distribution = apply_coalescent_prior(coalescent_model, key, n_children, distribution)?;
    let distribution = apply_date_constraint(date_constraint.as_ref(), distribution, max_points)?;
    let distribution = apply_leaf_coalescent_factor(coalescent_model.filter(|_| n_children == 0), key, distribution)?;

    let subtree = if matches!(distribution, Distribution::Empty) {
      None
    } else {
      let distribution = distribution
        .repeak()
        .wrap_err_with(|| format!("When re-peaking the time distribution of node {key}"))?;
      Some(Arc::new(distribution))
    };

    let to_parent = tree
      .parent(key)
      .map(|(_, edge_key)| {
        send_backward_message(
          bad_branch,
          subtree.as_deref(),
          edge_key,
          &branches[&edge_key],
          max_points,
        )
      })
      .transpose()?
      .flatten();
    let output = NodeBackward {
      subtree,
      to_parent: to_parent.clone(),
    };
    Ok((output, to_parent))
  }
}

fn combine_child_messages(
  messages: &[Arc<Distribution<NegLog>>],
  max_points: MaxGridPoints,
) -> Result<Distribution<NegLog>, Report> {
  let factors: Vec<&Distribution<NegLog>> = messages
    .iter()
    .map(Arc::as_ref)
    .filter(|message| !matches!(message, Distribution::Empty))
    .collect();
  if factors.is_empty() {
    return Ok(Distribution::Empty);
  }
  distribution_product(&factors, max_points)
}

fn apply_coalescent_prior(
  coalescent_model: Option<&CoalescentPrior>,
  key: GraphNodeKey,
  n_children: usize,
  distribution: Distribution<NegLog>,
) -> Result<Distribution<NegLog>, Report> {
  let Some(model) = coalescent_model else {
    return Ok(distribution);
  };
  if n_children == 0 || matches!(distribution, Distribution::Empty) {
    return Ok(distribution);
  }
  distribution_multiply_by_fn(
    &distribution,
    |time| model.node_contribution(key, n_children, time),
    model.node_contribution_exterior_slopes(key, n_children),
  )
}

fn apply_date_constraint(
  date_constraint: Option<&Arc<Distribution<NegLog>>>,
  distribution: Distribution<NegLog>,
  max_points: MaxGridPoints,
) -> Result<Distribution<NegLog>, Report> {
  let Some(constraint) = date_constraint else {
    return Ok(distribution);
  };
  if matches!(distribution, Distribution::Empty) {
    return Ok(constraint.as_ref().clone());
  }
  distribution_multiplication(&distribution, constraint, max_points)
}

fn apply_leaf_coalescent_factor(
  coalescent_model: Option<&CoalescentPrior>,
  key: GraphNodeKey,
  distribution: Distribution<NegLog>,
) -> Result<Distribution<NegLog>, Report> {
  let Some(model) = coalescent_model else {
    return Ok(distribution);
  };
  let distribution = match distribution.time_bounds() {
    Some((start, end)) if matches!(distribution, Distribution::Range(_)) => {
      let n_points =
        (f64_to_usize_saturating(((end - start) / model.time_resolution()).ceil()) + 1).clamp(2, GRID_POINTS);
      let dx = (end - start) / (n_points - 1).to_f64_exact();
      let values = Array1::from_elem(n_points, distribution.eval(start)?);
      Distribution::Function(
        DistributionFunction::from_start_dx_values(start, dx, values)?
          .with_left_extrap(BoundaryBehavior::Hard)?
          .with_right_extrap(BoundaryBehavior::Hard)?,
      )
    },
    _ => distribution,
  };
  distribution_multiply_by_fn(
    &distribution,
    |time| model.node_contribution(key, 0, time),
    model.node_contribution_exterior_slopes(key, 0),
  )
}

fn send_backward_message(
  bad_branch: bool,
  subtree: Option<&Distribution<NegLog>>,
  edge_key: GraphEdgeKey,
  branch: &BranchLikelihood,
  max_points: MaxGridPoints,
) -> Result<TimeDistribution, Report> {
  if bad_branch {
    return Ok(None);
  }
  let (Some(distribution), Some(branch_length_distribution)) = (subtree, &branch.distribution) else {
    return Ok(None);
  };

  let negated_branch = branch_length_distribution.negate()?;
  let message = convolve_across_edge(distribution, &negated_branch, Side::Left, EPS, GRID_POINTS, max_points)
    .wrap_err_with(|| format!("When sending the time message backward along edge {edge_key}"))?;

  Ok(Some(Arc::new(message)))
}
