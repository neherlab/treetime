use crate::clock::date_constraints::DateConstraints;
use crate::coalescent::coalescent::CoalescentModel;
use crate::timetree::inference::runner::{EPS, GRID_POINTS};
use crate::timetree::inference::time_inference::{BranchLikelihood, TimeBackward};
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_distribution::NegLog;
use treetime_distribution::convolve_across_edge;
use treetime_distribution::distribution_multiplication;
use treetime_distribution::distribution_multiply_by_fn;
use treetime_distribution::distribution_product;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pass::{
  GraphMapOutputs, GraphPass, GraphPassBackwardContext, GraphPassChildBackward, GraphPassNodeOutput,
};
use treetime_grid::Side;
use treetime_utils::make_internal_report;

pub(crate) fn propagate_distributions_backward(
  graph: &Graph,
  constraints: &DateConstraints,
  coalescent_model: Option<&CoalescentModel>,
  bad_branches: &BTreeMap<GraphNodeKey, bool>,
  branches: &BTreeMap<GraphEdgeKey, BranchLikelihood>,
) -> Result<TimeBackward, Report> {
  let pass = GraphPass::new(graph)?;
  let GraphMapOutputs { nodes, edges } = pass.map_backward(
    bad_branches,
    branches,
    |key| Err(make_internal_report!("Bad-branch flags are missing node {key}")),
    |context| propagate_distributions_backward_node(constraints, coalescent_model, &context),
  )?;
  Ok(TimeBackward {
    subtree: nodes,
    messages: edges,
  })
}

type TimeMessage = Option<Arc<Distribution<NegLog>>>;

fn propagate_distributions_backward_node(
  constraints: &DateConstraints,
  coalescent_model: Option<&CoalescentModel>,
  context: &GraphPassBackwardContext<'_, bool, BranchLikelihood, TimeMessage, TimeMessage>,
) -> Result<GraphPassNodeOutput<TimeMessage, TimeMessage>, Report> {
  let date_constraint = constraints.date_constraint(context.key);
  let messages = gather_child_messages(context.children);
  let distribution = combine_child_messages(&messages)?;
  let distribution = apply_coalescent_prior(coalescent_model, context.is_root, context.children.len(), distribution)?;
  let distribution = apply_date_constraint(date_constraint.as_ref(), distribution)?;

  let subtree = if matches!(distribution, Distribution::Empty) {
    None
  } else {
    let distribution = distribution
      .normalize()
      .wrap_err_with(|| format!("When normalizing the time distribution of node {}", context.key))?;
    Some(Arc::new(distribution))
  };

  let bad_branch = *context.input;
  let parent_message = send_backward_message(
    coalescent_model,
    context.is_leaf,
    bad_branch,
    subtree.as_deref(),
    context.parent_edge,
  )?;
  Ok(GraphPassNodeOutput {
    node: subtree,
    parent_message,
  })
}

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
fn gather_child_messages(
  children: &[GraphPassChildBackward<'_, TimeMessage, TimeMessage>],
) -> Vec<Arc<Distribution<NegLog>>> {
  children
    .iter()
    .filter_map(|child| {
      child
        .edge
        .expect("Non-root indexed node must own its parent edge")
        .clone()
    })
    .collect()
}

fn combine_child_messages(messages: &[Arc<Distribution<NegLog>>]) -> Result<Distribution<NegLog>, Report> {
  let factors: Vec<&Distribution<NegLog>> = messages
    .iter()
    .map(Arc::as_ref)
    .filter(|message| !matches!(message, Distribution::Empty))
    .collect();
  if factors.is_empty() {
    return Ok(Distribution::Empty);
  }
  distribution_product(&factors)
}

fn apply_coalescent_prior(
  coalescent_model: Option<&CoalescentModel>,
  is_root: bool,
  n_children: usize,
  distribution: Distribution<NegLog>,
) -> Result<Distribution<NegLog>, Report> {
  let Some(model) = coalescent_model else {
    return Ok(distribution);
  };
  if matches!(distribution, Distribution::Empty) {
    return Ok(distribution);
  }
  distribution_multiply_by_fn(&distribution, |time| {
    if is_root {
      model.root_contribution(time, n_children)
    } else {
      model.internal_contribution(time, n_children)
    }
  })
}

fn apply_date_constraint(
  date_constraint: Option<&Arc<Distribution<NegLog>>>,
  distribution: Distribution<NegLog>,
) -> Result<Distribution<NegLog>, Report> {
  let Some(constraint) = date_constraint else {
    return Ok(distribution);
  };
  if matches!(distribution, Distribution::Empty) {
    return Ok(constraint.as_ref().clone());
  }
  distribution_multiplication(&distribution, constraint)
}

fn send_backward_message(
  coalescent_model: Option<&CoalescentModel>,
  is_leaf: bool,
  bad_branch: bool,
  subtree: Option<&Distribution<NegLog>>,
  parent_edge: Option<(GraphEdgeKey, &BranchLikelihood)>,
) -> Result<Option<TimeMessage>, Report> {
  let Some((edge_key, branch)) = parent_edge else {
    return Ok(None);
  };
  if bad_branch {
    return Ok(Some(None));
  }
  let (Some(distribution), Some(branch_length_distribution)) = (subtree, &branch.distribution) else {
    return Ok(Some(None));
  };

  let leaf_weighted = if is_leaf && let Some(model) = coalescent_model {
    Some(distribution_multiply_by_fn(distribution, |time| {
      Ok(model.leaf_contribution(time))
    })?)
  } else {
    None
  };
  let outgoing = leaf_weighted.as_ref().unwrap_or(distribution);

  let negated_branch = branch_length_distribution.negate()?;
  let message = convolve_across_edge(outgoing, &negated_branch, Side::Left, EPS, GRID_POINTS)
    .wrap_err_with(|| format!("When sending the time message backward along edge {edge_key}"))?;

  Ok(Some(Some(Arc::new(message))))
}
