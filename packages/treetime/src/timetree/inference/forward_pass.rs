use crate::partition::indexed_pass::{IndexedPassDependencies, IndexedPassSlot, with_indexed_graph_payloads};
use crate::payload::traits::{TimetreeEdge, TimetreeNode};
use eyre::Report;
use std::sync::Arc;
use treetime_distribution::BoundaryBehavior;
use treetime_distribution::distribution_convolution;
use treetime_distribution::distribution_division;
use treetime_distribution::distribution_multiplication;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;

pub fn propagate_distributions_forward<N, E, D>(graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode + TimetreeNode + Default,
  E: GraphEdge + TimetreeEdge + Default,
  D: Send + Sync,
{
  with_indexed_graph_payloads(graph, |pass| {
    pass.try_for_each_forward(|dependencies, slot| propagate_distributions_forward_slot(graph, dependencies, slot))
  })
}

fn propagate_distributions_forward_slot<N, E, D>(
  graph: &Graph<N, E, D>,
  dependencies: &IndexedPassDependencies<N, E>,
  slot: &mut IndexedPassSlot<N, E>,
) -> Result<(), Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
  D: Send + Sync,
{
  refine_distribution_from_parent(graph, dependencies, slot)?;

  // Independent marginal modes can invert even though the forward pass commits each
  // parent before visiting its children. Current behavior projects non-leaf internal
  // point estimates to the committed parent time; leaves retain their observed dates.
  // This changes the point estimate without recomputing the posterior. The statistical
  // contract remains open in kb/issues/M-timetree-marginal-node-times-can-violate-topology.md.
  let parent_time = (!graph.is_leaf(slot.key))
    .then(|| parent_time(dependencies, slot))
    .flatten();
  set_likely_time(&mut slot.node, parent_time);
  Ok(())
}

/// Committed time of the slot's parent, if it has one and it is set.
fn parent_time<N, E>(dependencies: &IndexedPassDependencies<N, E>, slot: &IndexedPassSlot<N, E>) -> Option<f64>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
{
  let parent_key = slot.parent_key?;
  dependencies.node(parent_key).time()
}

fn set_likely_time(node: &mut impl TimetreeNode, parent_time: Option<f64>) {
  let time = node
    .time_distribution()
    .as_ref()
    .and_then(|time_dist| time_dist.likely_time());

  if let Some(mut time) = time {
    if let Some(parent_time) = parent_time {
      time = time.max(parent_time);
    }
    node.set_time(Some(time));
  }
}

fn refine_distribution_from_parent<N, E, D>(
  graph: &Graph<N, E, D>,
  dependencies: &IndexedPassDependencies<N, E>,
  slot: &mut IndexedPassSlot<N, E>,
) -> Result<(), Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
  D: Send + Sync,
{
  if slot.parent_key.is_none() {
    return Ok(());
  }

  // Do not overwrite leaf time_distribution (date constraint)
  if graph.is_leaf(slot.key) {
    return Ok(());
  }

  if let Some(parent_key) = slot.parent_key {
    let parent = dependencies.node(parent_key);
    let (_, edge) = slot
      .parent_edge
      .as_ref()
      .expect("Non-root indexed node must own its parent edge");

    if let (Some(parent_time_dist), Some(branch_dist), Some(subtree_dist)) = (
      parent.time_distribution(),
      edge.branch_length_distribution(),
      slot.node.time_distribution(),
    ) {
      let parent_except_subtree = if let Some(msg_to_parent) = edge.msg_to_parent() {
        distribution_division(parent_time_dist, msg_to_parent)?
      } else {
        parent_time_dist.as_ref().clone()
      };

      // Tail policy for the forward message (kb/decisions/timetree-inference-pass-boundary-tails.md).
      // The parent's time is a hard lower bound (left tail Zero); there is no upper bound from
      // the parent side on how far in the future the node could be (right tail Constant).
      let dist_from_parent = distribution_convolution(&parent_except_subtree, branch_dist)?
        .with_left_extrap(BoundaryBehavior::Zero)?
        .with_right_extrap(BoundaryBehavior::Constant)?;
      // Normalize to prevent numerical underflow: the backward pass stores normalized
      // distributions (max=1.0), and the convolution/division can produce arbitrary scales.
      // Without normalization, values accumulate downward across tree depth.
      let combined = distribution_multiplication(&dist_from_parent, subtree_dist)?.normalize();
      slot.node.set_time_distribution(Some(Arc::new(combined)));
    } else if let (Some(parent_time_dist), Some(branch_dist)) =
      (parent.time_distribution(), edge.branch_length_distribution())
    {
      // Same forward tail policy as the two-parent-message branch above. normalize() rebuilds
      // the grid, so the tail policy is applied afterwards, not before.
      let dist_from_parent = distribution_convolution(parent_time_dist, branch_dist)?
        .normalize()
        .with_left_extrap(BoundaryBehavior::Zero)?
        .with_right_extrap(BoundaryBehavior::Constant)?;
      slot.node.set_time_distribution(Some(Arc::new(dist_from_parent)));
    }
  }

  Ok(())
}
