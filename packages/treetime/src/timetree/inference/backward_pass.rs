use crate::partition::indexed_pass::{IndexedPassSlot, with_indexed_graph_payloads};
use crate::payload::traits::{TimetreeEdge, TimetreeNode};
use eyre::Report;
use indexmap::IndexMap;
use rayon::prelude::*;
use std::sync::Arc;
use treetime_distribution::distribution_convolution;
use treetime_distribution::distribution_multiplication;
use treetime_distribution::{Distribution, DistributionNegLog};
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey};

/// Propagates time distributions backward from leaves to root.
///
/// If coalescent_contributions is provided, multiplies each node's distribution
/// with its precalculated coalescent contribution.
pub fn propagate_distributions_backward<N, E, D>(
  graph: &Graph<N, E, D>,
  coalescent_contributions: Option<&IndexMap<GraphNodeKey, Arc<DistributionNegLog>>>,
) -> Result<(), Report>
where
  N: GraphNode + TimetreeNode + Default,
  E: GraphEdge + TimetreeEdge + Default,
  D: Send + Sync,
{
  with_indexed_graph_payloads(graph, |pass| {
    pass.try_for_each_backward_frontier(|node_indices, _, _, completed, frontier| {
      frontier.par_iter_mut().try_for_each(|slot| {
        propagate_distributions_backward_slot(graph, coalescent_contributions, node_indices, completed, slot)
      })
    })
  })
}

/// Computes time distribution for a single internal node from its children.
fn propagate_distributions_backward_slot<N, E, D>(
  graph: &Graph<N, E, D>,
  coalescent_contribs: Option<&IndexMap<GraphNodeKey, Arc<DistributionNegLog>>>,
  node_indices: &[Option<usize>],
  completed: &[IndexedPassSlot<N, E>],
  slot: &mut IndexedPassSlot<N, E>,
) -> Result<(), Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
  D: Send + Sync,
{
  let coalescent_contribution = if graph.is_leaf(slot.key) {
    None
  } else {
    coalescent_contribs.and_then(|contributions| contributions.get(&slot.key))
  };
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
      } else if let Some(contribution) = coalescent_contribution {
        let parent_message = parent_message.to_neglog();
        distribution_multiplication(contribution, &parent_message)?.to_plain_normalized()
      } else {
        parent_message.as_ref().clone()
      });
    }
  }

  if result.is_none()
    && let Some(contribution) = coalescent_contribution
  {
    result = Some(contribution.to_plain_normalized());
  }

  // Store final distribution on node
  if let Some(dist) = result {
    slot.node.set_time_distribution(Some(Arc::new(dist)));
  }

  if !slot.node.bad_branch()
    && let Some((_, edge)) = slot.parent_edge.as_mut()
    && let (Some(branch_dist), Some(node_time_dist)) =
      (edge.branch_length_distribution(), slot.node.time_distribution())
  {
    let negated_branch_dist = branch_dist.negate()?;
    let parent_message = distribution_convolution(node_time_dist.as_ref(), &negated_branch_dist)?;
    edge.set_msg_to_parent(Some(Arc::new(parent_message)));
  }

  Ok(())
}
