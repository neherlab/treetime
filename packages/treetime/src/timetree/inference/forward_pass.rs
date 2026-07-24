use crate::partition::indexed_pass::{IndexedPassSlot, with_indexed_graph_payloads};
use crate::payload::traits::{TimetreeEdge, TimetreeNode};
use eyre::Report;
use rayon::prelude::*;
use std::sync::Arc;
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
    pass.try_for_each_forward_frontier(|node_indices, _, _, completed_start, frontier, completed| {
      frontier.par_iter_mut().try_for_each(|slot| {
        propagate_distributions_forward_slot(graph, node_indices, completed_start, completed, slot)
      })
    })
  })
}

fn propagate_distributions_forward_slot<N, E, D>(
  graph: &Graph<N, E, D>,
  node_indices: &[Option<usize>],
  completed_start: usize,
  completed: &[IndexedPassSlot<N, E>],
  slot: &mut IndexedPassSlot<N, E>,
) -> Result<(), Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
  D: Send + Sync,
{
  refine_distribution_from_parent(graph, node_indices, completed_start, completed, slot)?;

  // The forward pass fixes each parent's time before its children. Clamp an internal
  // child's committed time up to its parent's so the point estimates respect
  // ancestor-before-descendant ordering, which the coalescent and other consumers
  // rely on. The marginals themselves already respect the physics (strictly positive
  // branch lengths keep each parent's support below its descendant tips), but they
  // integrate out other node times, so their independent modes can invert. Leaves
  // keep their observed date and are never clamped; the same positive-branch bound
  // keeps a parent's mode below its descendant tips, so clamping cannot push a child
  // past its own tips.
  let parent_time = (!graph.is_leaf(slot.key))
    .then(|| parent_time(node_indices, completed_start, completed, slot))
    .flatten();
  set_likely_time(&mut slot.node, parent_time);
  Ok(())
}

/// Committed time of the slot's parent, if it has one and it is set.
fn parent_time<N, E>(
  node_indices: &[Option<usize>],
  completed_start: usize,
  completed: &[IndexedPassSlot<N, E>],
  slot: &IndexedPassSlot<N, E>,
) -> Option<f64>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
{
  let parent_key = slot.parent_key?;
  let parent_index = node_indices[parent_key.as_usize()].expect("Indexed parent must have a slot");
  completed[parent_index - completed_start].node.time()
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
  node_indices: &[Option<usize>],
  completed_start: usize,
  completed: &[IndexedPassSlot<N, E>],
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
    let parent_index = node_indices[parent_key.as_usize()].expect("Indexed parent must have a slot");
    let parent = &completed[parent_index - completed_start].node;
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

      let dist_from_parent = distribution_convolution(&parent_except_subtree, branch_dist)?;
      // Normalize to prevent numerical underflow: the backward pass stores normalized
      // distributions (max=1.0), and the convolution/division can produce arbitrary scales.
      // Without normalization, values accumulate downward across tree depth.
      let combined = distribution_multiplication(&dist_from_parent, subtree_dist)?.normalize();
      slot.node.set_time_distribution(Some(Arc::new(combined)));
    } else if let (Some(parent_time_dist), Some(branch_dist)) =
      (parent.time_distribution(), edge.branch_length_distribution())
    {
      let dist_from_parent = distribution_convolution(parent_time_dist, branch_dist)?.normalize();
      slot.node.set_time_distribution(Some(Arc::new(dist_from_parent)));
    }
  }

  Ok(())
}
