use crate::coalescent::coalescent::CoalescentModel;
use crate::payload::traits::{TimetreeEdge, TimetreeNode};
use crate::timetree::inference::runner::{EPS, GRID_POINTS};
use eyre::Report;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_distribution::NegLog;
use treetime_distribution::convolve_across_edge;
use treetime_distribution::distribution_multiplication;
use treetime_distribution::distribution_multiply_by_fn;
use treetime_distribution::distribution_product;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::indexed_pass::{IndexedPassDependencies, IndexedPassSlot, with_indexed_graph_payloads};
use treetime_graph::node::GraphNode;
use treetime_grid::Side;

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
    pass.try_for_each_backward(|dependencies, slot| {
      propagate_distributions_backward_slot(graph, coalescent_model, dependencies, slot)
    })
  })
}

/// Computes a node's time distribution and the backward message it sends to its parent.
///
/// The node is handled in two phases. First, the child backward messages are folded with the
/// coalescent prior and the input date constraint into the node's time distribution, which is stored
/// peak-normalized. Second, the distribution is convolved across the branch into the backward message
/// the parent folds in.
fn propagate_distributions_backward_slot<N, E, D>(
  graph: &Graph<N, E, D>,
  coalescent_model: Option<&CoalescentModel>,
  dependencies: &IndexedPassDependencies<N, E>,
  slot: &mut IndexedPassSlot<N, E>,
) -> Result<(), Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
  D: Send + Sync,
{
  // take child messages and date constraint --> determine xmin, xmax, and grid
  // evaluate coalescent on that grid (different for root, internal, child)
  // multiply messages, date constraint, and coalescent --> node distribution
  // regrid node distribution to sensible grid
  let messages = gather_child_messages(graph, dependencies, slot);
  let distribution = combine_child_messages(&messages)?;
  let distribution = apply_coalescent_prior(graph, coalescent_model, slot, distribution)?;
  let distribution = apply_date_constraint(slot, distribution)?;

  if !matches!(distribution, Distribution::Empty) {
    // Peak-normalize the combined posterior. Every downstream consumer (likely_time, quantile, and
    // the outgoing convolution via to_plain_normalized) is shift-invariant, so the peak offset
    // removed here has no effect on inferred times or likelihoods.
    let distribution = distribution.normalize();
    slot.node.set_time_distribution(Some(Arc::new(distribution)));
  }

  send_backward_message(coalescent_model, graph.is_leaf(slot.key), slot)
}

/// Gathers the backward messages from a node's good children.
///
/// Collecting the messages up front keeps the fold independent of the order the children are visited.
/// A bad-branch child carries no usable message and is skipped.
fn gather_child_messages<N, E, D>(
  graph: &Graph<N, E, D>,
  dependencies: &IndexedPassDependencies<N, E>,
  slot: &IndexedPassSlot<N, E>,
) -> Vec<Arc<Distribution<NegLog>>>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
  D: Send + Sync,
{
  let mut messages = Vec::new();
  let node = graph.get_node(slot.key).expect("Indexed node must exist");
  for (child, _) in graph.children_of(&node.read_arc()) {
    let child_key = child.read_arc().key();
    let child = dependencies.slot(child_key);
    if child.node.bad_branch() {
      continue;
    }
    let edge = &child
      .parent_edge
      .as_ref()
      .expect("Non-root indexed node must own its parent edge")
      .1;
    if let Some(parent_message) = edge.msg_to_parent() {
      messages.push(Arc::clone(parent_message));
    }
  }
  messages
}

/// Combine the backward messages of a node's children into a single time distribution.
///
/// Multiplying time distributions is a product of independent factors, which [`distribution_product`]
/// forms directly: the gridded `Function` messages are co-located on one common working grid and
/// summed once (each resampled at most once regardless of fan-out, so the result is independent of
/// child order), and the exact `Point`/`Range` messages multiply in with no grid. An empty message
/// carries no probability and is not a legitimate backward message, so it is dropped before the
/// product. Returns `Empty` when no child left a usable message.
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

/// Multiply the node's role-specific coalescent prior into its folded child messages.
///
/// The root and internal contributions differ, and both scale with the node's child count. The leaf
/// coalescent factor is deliberately not applied here: it belongs to the outgoing message only, never
/// the stored node distribution, so it is added later when the message is formed in
/// [`send_backward_message`]. A node with no coalescent model, or an `Empty` distribution (no child
/// left a message), is returned unchanged.
fn apply_coalescent_prior<N, E, D>(
  graph: &Graph<N, E, D>,
  coalescent_model: Option<&CoalescentModel>,
  slot: &IndexedPassSlot<N, E>,
  distribution: Distribution<NegLog>,
) -> Result<Distribution<NegLog>, Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
  D: Send + Sync,
{
  let Some(model) = coalescent_model else {
    return Ok(distribution);
  };
  if matches!(distribution, Distribution::Empty) {
    return Ok(distribution);
  }
  let is_root = graph.is_root(slot.key);
  let n_children = graph.degree_out(slot.key)?;
  distribution_multiply_by_fn(&distribution, |time| {
    if is_root {
      model.root_contribution(time, n_children)
    } else {
      model.internal_contribution(time, n_children)
    }
  })
}

/// Multiply the node's input date constraint into its accumulated factors.
///
/// The date constraint is an independent factor of the node's posterior, so it multiplies whatever the
/// children have to say; for a leaf, which has no children, it is the whole distribution. Applying it
/// here rather than once at load time is what keeps the input recoverable: the forward pass refines the
/// time distribution of a node whose date is uncertain in place, and sending that refined distribution
/// back to the parent on the next round would count the parent's own message toward the node a second
/// time.
fn apply_date_constraint<N, E>(
  slot: &IndexedPassSlot<N, E>,
  distribution: Distribution<NegLog>,
) -> Result<Distribution<NegLog>, Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
{
  let Some(constraint) = slot.node.date_constraint() else {
    return Ok(distribution);
  };
  if matches!(distribution, Distribution::Empty) {
    return Ok(constraint.as_ref().clone());
  }
  distribution_multiplication(&distribution, constraint)
}

/// Convolves the node's stored distribution across the branch into a backward message for the parent.
///
/// The branch is negated (parent is older than the node). The left tail is soft (parent could be
/// arbitrarily far in the past); the right tail is hard (child's sampling date bounds the parent's
/// age). A leaf weights its outgoing message by the coalescent leaf factor, which belongs to the
/// message only, not the stored distribution.
fn send_backward_message<N, E>(
  coalescent_model: Option<&CoalescentModel>,
  is_leaf: bool,
  slot: &mut IndexedPassSlot<N, E>,
) -> Result<(), Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
{
  if slot.node.bad_branch() {
    return Ok(());
  }
  let Some(distribution) = slot.node.time_distribution() else {
    return Ok(());
  };
  let Some((_, edge)) = slot.parent_edge.as_mut() else {
    return Ok(());
  };
  let Some(branch_length_distribution) = edge.branch_length_distribution() else {
    return Ok(());
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
  let message = convolve_across_edge(outgoing, &negated_branch, Side::Left, EPS, GRID_POINTS)?;
  edge.set_msg_to_parent(Some(Arc::new(message)));

  Ok(())
}
