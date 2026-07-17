use crate::coalescent::coalescent::CoalescentModel;
use crate::partition::indexed_pass::{IndexedPassSlot, with_indexed_graph_payloads};
use crate::payload::traits::{TimetreeEdge, TimetreeNode};
use eyre::Report;
use rayon::prelude::*;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_distribution::distribution_apply_neg_log_weight;
use treetime_distribution::distribution_convolution;
use treetime_distribution::distribution_multiplication;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;

/// Propagates time distributions backward from leaves to root.
///
/// If a coalescent model is provided, applies one role-specific contribution
/// after all child messages have been combined.
pub fn propagate_distributions_backward<N, E, D>(
  graph: &Graph<N, E, D>,
  coalescent_model: Option<&CoalescentModel>,
) -> Result<(), Report>
where
  N: GraphNode + TimetreeNode + Default,
  E: GraphEdge + TimetreeEdge + Default,
  D: Send + Sync,
{
  with_indexed_graph_payloads(graph, |pass| {
    pass.try_for_each_backward_frontier(|node_indices, _, _, completed, frontier| {
      frontier.par_iter_mut().try_for_each(|slot| {
        propagate_distributions_backward_slot(graph, coalescent_model, node_indices, completed, slot)
      })
    })
  })
}

/// Computes time distribution for a single internal node from its children.
fn propagate_distributions_backward_slot<N, E, D>(
  graph: &Graph<N, E, D>,
  coalescent_model: Option<&CoalescentModel>,
  node_indices: &[Option<usize>],
  completed: &[IndexedPassSlot<N, E>],
  slot: &mut IndexedPassSlot<N, E>,
) -> Result<(), Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
  D: Send + Sync,
{
  let is_leaf = graph.is_leaf(slot.key);
  let is_root = slot.parent_edge.is_none();
  let n_children = graph
    .get_node(slot.key)
    .expect("Indexed node must exist")
    .read_arc()
    .outbound()
    .len();
  let mut result: Option<Distribution> = None;

  for (child, _) in graph.children_of(&graph.get_node(slot.key).expect("Indexed node must exist").read_arc()) {
    let child_key = child.read_arc().key();
    let child_index = node_indices[child_key.as_usize()].expect("Indexed child must have a slot");
    let child = &completed[child_index];
    if child.node.bad_branch() {
      continue;
    }
    let edge = &child
      .parent_edge
      .as_ref()
      .expect("Non-root indexed node must own its parent edge")
      .1;
    if let Some(parent_message) = edge.msg_to_parent() {
      // Combine messages from all children using multiplication (intersection of constraints).
      // Normalize after each step to prevent numerical underflow in plain probability space:
      // without normalization, the peak decays as ~(peak)^N across N children, underflowing
      // to zero for large trees. Normalization (max=1.0) is safe because all downstream
      // consumers (likely_time, quantile, hpd_region) are scale-invariant.
      result = Some(if let Some(current) = result {
        distribution_multiplication(&current, parent_message)?.normalize()
      } else {
        parent_message.as_ref().clone()
      });
    }
  }

  if let (Some(model), Some(distribution)) = (coalescent_model, result.as_ref()) {
    result = Some(if is_root {
      distribution_apply_neg_log_weight(distribution, |time| model.root_contribution(time, n_children))?
    } else {
      distribution_apply_neg_log_weight(distribution, |time| model.internal_contribution(time, n_children))?
    });
  }

  if let Some(dist) = result.as_ref() {
    // Leaves retain their observed date distribution. Their coalescent factor
    // belongs only to the temporary message convolved toward the parent.
    debug_assert!(!is_leaf);
    slot.node.set_time_distribution(Some(Arc::new(dist.clone())));
  }

  let outgoing_distribution = if is_leaf {
    let distribution = slot.node.time_distribution();
    match (coalescent_model, distribution) {
      (Some(model), Some(distribution)) => Some(Arc::new(distribution_apply_neg_log_weight(
        distribution.as_ref(),
        |time| Ok(model.leaf_contribution(time)),
      )?)),
      (_, distribution) => distribution.clone(),
    }
  } else {
    slot.node.time_distribution().clone()
  };

  if !slot.node.bad_branch()
    && let Some((_, edge)) = slot.parent_edge.as_mut()
    && let (Some(branch_dist), Some(node_time_dist)) = (edge.branch_length_distribution(), outgoing_distribution)
  {
    let negated_branch_dist = branch_dist.negate()?;
    let parent_message = distribution_convolution(node_time_dist.as_ref(), &negated_branch_dist)?;
    edge.set_msg_to_parent(Some(Arc::new(parent_message)));
  }

  Ok(())
}
