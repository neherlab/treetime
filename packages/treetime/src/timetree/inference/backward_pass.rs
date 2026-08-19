use crate::coalescent::coalescent::CoalescentModel;
use crate::partition::indexed_pass::{IndexedPassDependencies, IndexedPassSlot, with_indexed_graph_payloads};
use crate::payload::traits::{TimetreeEdge, TimetreeNode};
use crate::timetree::inference::runner::{EPS, GRID_POINTS};
use eyre::{Report, WrapErr};
use std::cmp::Ordering;
use std::sync::Arc;
use treetime_distribution::BoundaryBehavior;
use treetime_distribution::Distribution;
use treetime_distribution::DistributionFunction;
use treetime_distribution::NegLog;
use treetime_distribution::convolve_across_edge;
use treetime_distribution::distribution_add_neg_log_weight;
use treetime_distribution::distribution_multiplication;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_grid::{DEFAULT_TAIL_FIT_POINTS, GridFn, Side, SoftTailLaw};
use treetime_utils::{make_internal_error, make_internal_report};

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
    slot
      .node
      .set_time_distribution(Some(Arc::new(distribution.normalize())));
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
  let mut distribution = combine_child_messages(&messages)?;

  // Multiply in the coalescent prior once the child messages are folded. The root and internal
  // contributions differ, and both scale with the node's child count. The leaf coalescent factor is
  // deliberately not applied here: it belongs to the outgoing message only, never the stored node
  // distribution, so it is added later in `outgoing_node_distribution`.
  if let (Some(model), Some(current)) = (coalescent_model, distribution.as_ref()) {
    let n_children = graph
      .get_node(slot.key)
      .expect("Indexed node must exist")
      .read_arc()
      .outbound()
      .len();
    distribution = Some(if slot.parent_edge.is_none() {
      distribution_add_neg_log_weight(current, |time| model.root_contribution(time, n_children))?
    } else {
      distribution_add_neg_log_weight(current, |time| model.internal_contribution(time, n_children))?
    });
  }

  // Multiply in the input date constraint. It is an independent factor of the node's posterior, so it
  // multiplies whatever the children have to say; for a leaf, which has no children, it is the whole
  // distribution. Doing this here rather than once at load time is what keeps the input recoverable:
  // the forward pass refines the time distribution of a node whose date is uncertain in place, and
  // sending that refined distribution back to the parent on the next round would count the parent's
  // own message toward the node a second time.
  if let Some(constraint) = slot.node.date_constraint() {
    distribution = Some(match distribution {
      Some(current) => distribution_multiplication(&current, constraint)?,
      None => constraint.as_ref().clone(),
    });
  }

  Ok(distribution)
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
/// Multiplying time distributions is pointwise, and under `NegLog` that is addition of ordinates:
/// exact, associative, and independent of the order the operands are combined. Gridded `Function`
/// messages must share a grid before they can be added, so they are summed on one common working
/// grid, each resampled exactly once (see [`sum_function_messages`]). The remaining `Point` and
/// `Range` messages are exact factors that need no grid; their product is formed directly and then
/// multiplied into the function sum in a single step.
///
/// This replaces the old left-to-right fold, which multiplied and normalized once per child. That
/// fold resampled the accumulator on every step (interpolation drift that grew with the fan-out)
/// and imposed a resampling order. Here every message is resampled at most once regardless of
/// fan-out, and the result no longer depends on child order.
fn combine_child_messages(messages: &[Arc<Distribution<NegLog>>]) -> Result<Option<Distribution<NegLog>>, Report> {
  let mut functions = Vec::new();
  let mut factors = Vec::new();
  for message in messages {
    match message.as_ref() {
      Distribution::Function(function) => functions.push(function),
      // An empty message carries no probability and is not a legitimate backward message; skip it.
      Distribution::Empty => {},
      point_or_range => factors.push(point_or_range),
    }
  }

  let function_sum = if functions.is_empty() {
    None
  } else {
    Some(sum_function_messages(&functions)?)
  };

  // Exact product of all point/range factors, formed with no grid so it introduces no resampling.
  let mut factor_product: Option<Distribution<NegLog>> = None;
  for factor in factors {
    factor_product = Some(match factor_product {
      Some(current) => distribution_multiplication(&current, factor)?,
      None => factor.clone(),
    });
  }

  Ok(match (function_sum, factor_product) {
    (Some(summed), Some(product)) => Some(distribution_multiplication(&summed, &product)?),
    (Some(distribution), None) | (None, Some(distribution)) => Some(distribution),
    (None, None) => None,
  })
}

/// Sum neg-log `Function` messages on one common working grid.
///
/// The grid is chosen once, following the hard/soft boundary rule of the log-space design
/// (`kb/proposals/distribution-log-space-and-hard-soft-boundaries.md`, Parts B and D):
///
/// - hard (or undeclared) side: the tightest bound wins, because a hard boundary is a fact that
///   probability is zero beyond it, so the product can be non-zero only where every operand is;
/// - soft side: the loosest bound wins, because a soft tail law continues past the grid edge, so
///   every operand stays evaluable out to the farthest edge;
/// - spacing: the finest operand's `dx`, so no operand loses resolution.
///
/// Each message is resampled onto this grid exactly once (tail-preserving, via
/// [`DistributionFunction::resample_range_n_points`], which carries each operand's own tail law to
/// the points beyond its support). Because every message then lands on the same grid, the fold is a
/// plain elementwise sum of ordinates: exact, and independent of child order.
fn sum_function_messages(functions: &[&DistributionFunction<f64, NegLog>]) -> Result<Distribution<NegLog>, Report> {
  let left = combine_lower_bound(functions.iter().map(|f| (f.x_min(), f.left_extrap())));
  let right = combine_upper_bound(functions.iter().map(|f| (f.x_max(), f.right_extrap())));
  let dx = functions
    .iter()
    .map(|f| f.dx())
    .reduce(f64::min)
    .expect("at least one function message");

  if left.partial_cmp(&right) != Some(Ordering::Less) {
    // A soft boundary never separates domains, so with the backward messages' soft-left tails this
    // is unreachable; only genuinely disjoint hard domains could reach it. Guarding it keeps an
    // empty product from silently poisoning every ancestor to the root (the motivating defect).
    return make_internal_error!(
      "Backward child fold produced an empty working grid [{left}, {right}]: the messages' hard domains are disjoint"
    );
  }

  // Discretize `[left, right]` at the finest operand spacing `dx`, matching how
  // `multiply_function_function` builds its intersection grid (`round(range / dx) + 1` points over an
  // inclusive grid). This intermediate fold grid is a transient fidelity grid, not a stored
  // distribution: it stays at least as fine as its sharpest operand so no operand loses resolution
  // before the sum (design D3). The combined posterior is only peak-normalized after the fold; mass
  // sizing happens later, once per message, when the outgoing message is regridded across an edge.
  let resampled = functions
    .iter()
    .map(|f| f.resample_range_dx((left, right), dx))
    .collect::<Result<Vec<_>, Report>>()?;

  // Every resampled message shares the same [left, x_max] / dx grid, so the fold is elementwise.
  let mut summed = resampled[0].y().clone();
  for function in &resampled[1..] {
    summed += function.y();
  }

  // Rebuild on the resampled grid's own extent (`resample_range_dx` ends at `left + (n - 1) * dx`,
  // which can fall a rounding step short of `right`), not on the nominal `[left, right]`, so the fold
  // result keeps the exact grid the operands were resampled onto.
  let summed = DistributionFunction::from_range_values((resampled[0].x_min(), resampled[0].x_max()), summed)?;
  let left_tail = derive_summed_tail(functions.iter().map(|f| f.left_extrap()), summed.grid_fn(), Side::Left)?;
  let right_tail = derive_summed_tail(
    functions.iter().map(|f| f.right_extrap()),
    summed.grid_fn(),
    Side::Right,
  )?;
  let summed = summed.with_left_extrap(left_tail)?.with_right_extrap(right_tail)?;
  Ok(Distribution::Function(summed))
}

/// Lower (left) working-grid bound: the tightest hard bound when any operand terminates the domain
/// on the left, otherwise the loosest soft bound.
fn combine_lower_bound(operands: impl Iterator<Item = (f64, BoundaryBehavior)>) -> f64 {
  let mut hard: Option<f64> = None;
  let mut soft: Option<f64> = None;
  for (x_min, tail) in operands {
    if is_restricting(tail) {
      hard = Some(hard.map_or(x_min, |h| h.max(x_min)));
    } else {
      soft = Some(soft.map_or(x_min, |s| s.min(x_min)));
    }
  }
  hard.or(soft).expect("at least one function message")
}

/// Upper (right) working-grid bound: the tightest hard bound when any operand terminates the domain
/// on the right, otherwise the loosest soft bound.
fn combine_upper_bound(operands: impl Iterator<Item = (f64, BoundaryBehavior)>) -> f64 {
  let mut hard: Option<f64> = None;
  let mut soft: Option<f64> = None;
  for (x_max, tail) in operands {
    if is_restricting(tail) {
      hard = Some(hard.map_or(x_max, |h| h.min(x_max)));
    } else {
      soft = Some(soft.map_or(x_max, |s| s.max(x_max)));
    }
  }
  hard.or(soft).expect("at least one function message")
}

/// Whether a tail terminates the evaluable domain on its side. `Hard`/`HardApproach` is zero
/// probability beyond the edge and `Error` is undefined beyond it, so all restrict the working grid;
/// the soft `Linear` law continues past the edge and does not.
fn is_restricting(tail: BoundaryBehavior) -> bool {
  tail.is_hard()
}

/// Derive the concrete tail policy of a summed message on one side from the operands' tail classes.
///
/// Mirrors the multiplication tail rule (`treetime_distribution::distribution_ops::multiply`) as a
/// class lattice `Error > Hard > Linear`: the product is soft only when every operand is soft, a
/// hard bound dominates any soft one, and an undeclared (`Error`) bound dominates all. The backward
/// messages carry no fitted law, so the hard and undeclared classes need none here (`Hard`,
/// `Error`). Only the soft log-linear class needs a law, which is fit from the summed data on
/// `side` -- a grid too degenerate to fit is an error rather than a silent flat fallback.
///
/// Backward messages carry a fitted `Linear` tail on the left (soft) side and `Hard` on the right,
/// so a fan-in of such messages reaches the soft-fit branch here and re-fits the summed left tail
/// from the combined grid; the right side stays `Hard`.
fn derive_summed_tail(
  tails: impl Iterator<Item = BoundaryBehavior>,
  summed: &GridFn<f64>,
  side: Side,
) -> Result<BoundaryBehavior, Report> {
  let mut any_error = false;
  let mut any_hard = false;
  let mut any_linear = false;
  for tail in tails {
    match tail {
      BoundaryBehavior::Error => any_error = true,
      BoundaryBehavior::Hard | BoundaryBehavior::HardApproach(_) => any_hard = true,
      BoundaryBehavior::Linear(_) => any_linear = true,
    }
  }

  if any_error {
    return Ok(BoundaryBehavior::Error);
  }
  if any_hard {
    return Ok(BoundaryBehavior::Hard);
  }
  if any_linear {
    let law = SoftTailLaw::fit(summed, side, DEFAULT_TAIL_FIT_POINTS).wrap_err_with(|| {
      make_internal_report!("Backward child fold cannot fit a soft tail on the {side:?} side from the summed grid")
    })?;
    return Ok(BoundaryBehavior::Linear(law));
  }
  // Every operand tail falls into one of the three classes above, so reaching here means the fold
  // received no operand messages -- the caller guarantees at least one, so this is an internal bug.
  make_internal_error!("Backward child fold on the {side:?} side received no operand tails")
}

/// Convolves a node's distribution across the branch to its parent and sets the backward message.
///
/// This is the edge crossing of the backward pass. The node's outgoing distribution
/// ([`outgoing_node_distribution`]) is convolved with the negated branch-length distribution
/// ([`convolve_across_branch`]). A node with no parent edge, a bad branch, or a missing branch-length
/// or node distribution sends no message.
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
  let outgoing = outgoing_node_distribution(coalescent_model, graph.is_leaf(slot.key), &slot.node)?;

  if !slot.node.bad_branch()
    && let Some((_, edge)) = slot.parent_edge.as_mut()
    && let (Some(branch_length_distribution), Some(outgoing)) = (edge.branch_length_distribution(), outgoing)
  {
    let message = convolve_across_branch(&outgoing, branch_length_distribution)?;
    edge.set_msg_to_parent(Some(Arc::new(message)));
  }

  Ok(())
}

/// The node's distribution as it leaves toward its parent.
///
/// A leaf's coalescent factor belongs only to this outgoing message, never to the distribution stored
/// on the node, so it is added here rather than in [`multiply_node_factors`]. An internal node sends
/// its stored distribution unchanged.
fn outgoing_node_distribution<N: TimetreeNode>(
  coalescent_model: Option<&CoalescentModel>,
  is_leaf: bool,
  node: &N,
) -> Result<Option<Arc<Distribution<NegLog>>>, Report> {
  if is_leaf && let (Some(model), Some(distribution)) = (coalescent_model, node.time_distribution()) {
    // Leaf coalescent factor: added to this outgoing message only, so it never reaches the stored node
    // distribution (set by `multiply_node_factors` without it). It weights how the leaf informs its
    // parent, not the leaf's own time, so it stays out of the stored distribution the forward pass reuses.
    let with_leaf = distribution_add_neg_log_weight(distribution.as_ref(), |time| Ok(model.leaf_contribution(time)))?;
    return Ok(Some(Arc::new(with_leaf)));
  }
  Ok(node.time_distribution().clone())
}

/// Convolves an outgoing node distribution across the branch to form the backward message.
///
/// The branch is negated because the parent is older than the node, so its age is the node's time
/// minus the branch length. The convolution then picks its output grid by probability mass from the
/// operands and lands the message on it in a single regrid (see [`convolve_across_edge`]).
///
/// Tail policy (kb/decisions/distribution-tails-and-arithmetic.md): the parent could be arbitrarily
/// far in the past, so the left side is soft (`Side::Left`) -- a fitted log-linear tail that decays
/// with finite mass, keeping the quantile and HPD integrals well-defined. The child's sampling date is
/// a hard upper bound on the parent's age, so the right tail stays `Hard`.
fn convolve_across_branch(
  outgoing: &Distribution<NegLog>,
  branch_length_distribution: &Distribution<NegLog>,
) -> Result<Distribution<NegLog>, Report> {
  let negated_branch = branch_length_distribution.negate()?;
  convolve_across_edge(outgoing, &negated_branch, Side::Left, EPS, GRID_POINTS)
}
