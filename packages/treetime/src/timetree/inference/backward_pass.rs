use crate::coalescent::coalescent::CoalescentModel;
use crate::partition::indexed_pass::{IndexedPassDependencies, IndexedPassSlot, with_indexed_graph_payloads};
use crate::payload::traits::{TimetreeEdge, TimetreeNode};
use crate::timetree::inference::runner::{EPS, GRID_POINTS};
use eyre::Report;
use std::borrow::Cow;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_distribution::NegLog;
use treetime_distribution::convolve_across_edge;
use treetime_distribution::distribution_add_neg_log_weight;
use treetime_distribution::distribution_multiplication;
use treetime_distribution::distribution_product;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
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
/// The node is handled in two phases. Child multiplication folds the backward messages of the node's
/// children and multiplies in the coalescent prior and the input date constraint, giving the node's
/// time distribution ([`multiply_node_factors`]); that distribution is stored peak-normalized. The
/// node distribution is then convolved across the branch into the backward message the parent folds
/// in ([`send_message_to_parent`]).
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
  let node_distribution = multiply_node_factors(graph, coalescent_model, dependencies, slot)?;

  if let Some(distribution) = &node_distribution {
    // Peak-normalize the combined posterior and store it as-is. The child fold already produced a
    // compact grid (operand extents, finest operand pitch), so the stored node distribution needs no
    // mass sizing: that belongs to the edge crossing, where a message is regridded once as it is
    // convolved toward the parent. Every downstream consumer (likely_time, quantile, and the outgoing
    // convolution via to_plain_normalized) is shift-invariant, so the peak offset removed here has no
    // effect on inferred times or likelihoods.
    let distribution = distribution.normalize();
    slot.node.set_time_distribution(Some(Arc::new(distribution)));
  }

  send_message_to_parent(coalescent_model, graph, slot)
}

/// Multiplies the independent factors that make up a node's time distribution.
///
/// A node's posterior over its time is a product of independent factors, which under `NegLog` is a
/// sum of ordinates: the fold of the child backward messages ([`combine_child_messages`]), the
/// coalescent prior (root- or internal-specific), and the input date constraint. Returns `None` for a
/// node with neither children nor a date constraint, which therefore constrains its time not at all.
fn multiply_node_factors<N, E, D>(
  graph: &Graph<N, E, D>,
  coalescent_model: Option<&CoalescentModel>,
  dependencies: &IndexedPassDependencies<N, E>,
  slot: &IndexedPassSlot<N, E>,
) -> Result<Option<Distribution<NegLog>>, Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
  D: Send + Sync,
{
  let messages = gather_child_messages(graph, dependencies, slot);
  let distribution = combine_child_messages(&messages)?;
  let distribution = apply_coalescent_prior(graph, coalescent_model, slot, distribution)?;
  apply_date_constraint(slot, distribution)
}

/// Multiply the node's role-specific coalescent prior into its folded child messages.
///
/// The root and internal contributions differ, and both scale with the node's child count. The leaf
/// coalescent factor is deliberately not applied here: it belongs to the outgoing message only, never
/// the stored node distribution, so it is added later when the message is formed in
/// [`send_message_to_parent`]. A node with no coalescent model, or none of whose children left a
/// message, is returned unchanged.
fn apply_coalescent_prior<N, E, D>(
  graph: &Graph<N, E, D>,
  coalescent_model: Option<&CoalescentModel>,
  slot: &IndexedPassSlot<N, E>,
  distribution: Option<Distribution<NegLog>>,
) -> Result<Option<Distribution<NegLog>>, Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
  D: Send + Sync,
{
  match (coalescent_model, distribution) {
    (Some(model), Some(current)) => {
      let n_children = graph
        .get_node(slot.key)
        .expect("Indexed node must exist")
        .read_arc()
        .outbound()
        .len();
      let weighted = if slot.parent_edge.is_none() {
        distribution_add_neg_log_weight(&current, |time| model.root_contribution(time, n_children))?
      } else {
        distribution_add_neg_log_weight(&current, |time| model.internal_contribution(time, n_children))?
      };
      Ok(Some(weighted))
    },
    (_, distribution) => Ok(distribution),
  }
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
  distribution: Option<Distribution<NegLog>>,
) -> Result<Option<Distribution<NegLog>>, Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
{
  let Some(constraint) = slot.node.date_constraint() else {
    return Ok(distribution);
  };
  Ok(Some(match distribution {
    Some(current) => distribution_multiplication(&current, constraint)?,
    None => constraint.as_ref().clone(),
  }))
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
/// product; a node whose children all fell away constrains its time not at all and returns `None`.
fn combine_child_messages(messages: &[Arc<Distribution<NegLog>>]) -> Result<Option<Distribution<NegLog>>, Report> {
  let factors: Vec<&Distribution<NegLog>> = messages
    .iter()
    .map(Arc::as_ref)
    .filter(|message| !matches!(message, Distribution::Empty))
    .collect();
  if factors.is_empty() {
    return Ok(None);
  }
  Ok(Some(distribution_product(&factors)?))
}

/// Convolves a node's distribution across the branch to its parent and sets the backward message.
///
/// This is the edge crossing of the backward pass. The node's outgoing distribution is convolved with
/// the negated branch-length distribution. A node with no parent edge, a bad branch, or a missing
/// branch-length or node distribution sends no message.
///
/// The branch is negated because the parent is older than the node, so its age is the node's time minus
/// the branch length. The convolution then picks its output grid by probability mass from the operands
/// and lands the message on it in a single regrid (see [`convolve_across_edge`]).
///
/// Tail policy (kb/decisions/distribution-tails-and-arithmetic.md): the parent could be arbitrarily far
/// in the past, so the left side is soft (`Side::Left`) -- a fitted log-linear tail that decays with
/// finite mass, keeping the quantile and HPD integrals well-defined. The child's sampling date is a hard
/// upper bound on the parent's age, so the right tail stays `Hard`.
fn send_message_to_parent<N, E, D>(
  coalescent_model: Option<&CoalescentModel>,
  graph: &Graph<N, E, D>,
  slot: &mut IndexedPassSlot<N, E>,
) -> Result<(), Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
  D: Send + Sync,
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

  let outgoing = outgoing_distribution(coalescent_model, graph.is_leaf(slot.key), distribution)?;

  let negated_branch = branch_length_distribution.negate()?;
  let message = convolve_across_edge(&outgoing, &negated_branch, Side::Left, EPS, GRID_POINTS)?;
  edge.set_msg_to_parent(Some(Arc::new(message)));

  Ok(())
}

/// The node's distribution as it leaves toward its parent.
///
/// A leaf weights its outgoing message by the leaf coalescent factor. This factor captures how the leaf
/// informs its parent, not the leaf's own time, so it stays out of the stored distribution the forward
/// pass reuses (`multiply_node_factors` sets that one without it). Every other node -- internal, or any
/// node without a coalescent model -- sends its stored distribution unchanged, borrowed rather than cloned.
fn outgoing_distribution<'a>(
  coalescent_model: Option<&CoalescentModel>,
  is_leaf: bool,
  distribution: &'a Arc<Distribution<NegLog>>,
) -> Result<Cow<'a, Distribution<NegLog>>, Report> {
  if is_leaf && let Some(model) = coalescent_model {
    let weighted = distribution_add_neg_log_weight(distribution.as_ref(), |time| Ok(model.leaf_contribution(time)))?;
    Ok(Cow::Owned(weighted))
  } else {
    Ok(Cow::Borrowed(distribution.as_ref()))
  }
}
