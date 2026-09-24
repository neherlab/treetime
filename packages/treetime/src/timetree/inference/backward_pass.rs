use crate::clock::date_constraints::DateConstraints;
use crate::coalescent::coalescent::CoalescentModel;
use crate::timetree::inference::runner::{EPS, GRID_POINTS};
use crate::timetree::timetree_state::{DateEdgeState, DateNodeState, TimetreeState};
use eyre::{Report, WrapErr};
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_distribution::NegLog;
use treetime_distribution::convolve_across_edge;
use treetime_distribution::distribution_multiplication;
use treetime_distribution::distribution_multiply_by_fn;
use treetime_distribution::distribution_product;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::pass::{GraphPassBackwardContext, GraphPassChildBackward, GraphPassNodeOutput};
use treetime_grid::Side;

pub(crate) fn propagate_distributions_backward(
  graph: &Graph,
  constraints: &DateConstraints,
  coalescent_model: Option<&CoalescentModel>,
  state: &mut TimetreeState,
) -> Result<(), Report> {
  state.map_backward(graph, |context| {
    propagate_distributions_backward_node(constraints, coalescent_model, &context)
  })
}

fn propagate_distributions_backward_node(
  constraints: &DateConstraints,
  coalescent_model: Option<&CoalescentModel>,
  context: &GraphPassBackwardContext<'_, DateNodeState, DateEdgeState, DateNodeState, DateEdgeState>,
) -> Result<GraphPassNodeOutput<DateNodeState, DateEdgeState>, Report> {
  let mut node = context.input.clone();
  let date_constraint = constraints.date_constraints.get(&context.key).cloned().flatten();
  let messages = gather_child_messages(context.children);
  let distribution = combine_child_messages(&messages)?;
  let distribution = apply_coalescent_prior(coalescent_model, context.is_root, context.children.len(), distribution)?;
  let distribution = apply_date_constraint(date_constraint.as_ref(), distribution)?;

  if !matches!(distribution, Distribution::Empty) {
    let distribution = distribution
      .normalize()
      .wrap_err_with(|| format!("When normalizing the time distribution of node {}", context.key))?;
    node.time_distribution = Some(Arc::new(distribution));
  }

  let parent_message = send_backward_message(coalescent_model, context.is_leaf, &node, context.parent_edge)?;
  Ok(GraphPassNodeOutput { node, parent_message })
}

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
fn gather_child_messages(
  children: &[GraphPassChildBackward<'_, DateNodeState, DateEdgeState>],
) -> Vec<Arc<Distribution<NegLog>>> {
  let mut messages = Vec::new();
  for child in children {
    if child.node.bad_branch {
      continue;
    }
    let edge = child.edge.expect("Non-root indexed node must own its parent edge");
    if let Some(parent_message) = &edge.msg_to_parent {
      messages.push(Arc::clone(parent_message));
    }
  }
  messages
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
  node: &DateNodeState,
  parent_edge: Option<(GraphEdgeKey, &DateEdgeState)>,
) -> Result<Option<DateEdgeState>, Report> {
  let Some((edge_key, edge)) = parent_edge else {
    return Ok(None);
  };
  let mut edge = edge.clone();
  if node.bad_branch {
    return Ok(Some(edge));
  }
  let Some(distribution) = &node.time_distribution else {
    return Ok(Some(edge));
  };
  let Some(branch_length_distribution) = &edge.branch_length_distribution else {
    return Ok(Some(edge));
  };

  let leaf_weighted = if is_leaf && let Some(model) = coalescent_model {
    Some(distribution_multiply_by_fn(distribution.as_ref(), |time| {
      Ok(model.leaf_contribution(time))
    })?)
  } else {
    None
  };
  let outgoing = leaf_weighted.as_ref().unwrap_or_else(|| distribution.as_ref());

  let negated_branch = branch_length_distribution.negate()?;
  let message = convolve_across_edge(outgoing, &negated_branch, Side::Left, EPS, GRID_POINTS)
    .wrap_err_with(|| format!("When sending the time message backward along edge {edge_key}"))?;
  edge.msg_to_parent = Some(Arc::new(message));

  Ok(Some(edge))
}
