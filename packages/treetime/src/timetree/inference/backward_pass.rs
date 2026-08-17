use crate::coalescent::coalescent::CoalescentModel;
use crate::partition::indexed_pass::{IndexedPassDependencies, IndexedPassSlot, with_indexed_graph_payloads};
use crate::payload::traits::{TimetreeEdge, TimetreeNode};
use crate::timetree::inference::runner::{EPS, GRID_POINTS};
use crate::timetree::inference::tail_fit::fit_message_soft_tail;
use eyre::{Report, WrapErr};
use std::cmp::Ordering;
use std::sync::Arc;
use treetime_distribution::BoundaryBehavior;
use treetime_distribution::Distribution;
use treetime_distribution::DistributionFunction;
use treetime_distribution::NegLog;
use treetime_distribution::distribution_add_neg_log_weight;
use treetime_distribution::distribution_convolution;
use treetime_distribution::distribution_multiplication;
use treetime_distribution::rewindow_to_mass;
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

/// Computes time distribution for a single internal node from its children.
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
  let is_leaf = graph.is_leaf(slot.key);
  let is_root = slot.parent_edge.is_none();
  let n_children = graph
    .get_node(slot.key)
    .expect("Indexed node must exist")
    .read_arc()
    .outbound()
    .len();

  // Gather the backward messages from all good children, then combine them on one common grid (see
  // fn combine_child_messages). Collecting the messages first keeps the fold independent of the
  // order the children are visited.
  let mut messages: Vec<Arc<Distribution<NegLog>>> = Vec::new();
  for (child, _) in graph.children_of(&graph.get_node(slot.key).expect("Indexed node must exist").read_arc()) {
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

  let mut result = combine_child_messages(&messages)?;

  if let (Some(model), Some(distribution)) = (coalescent_model, result.as_ref()) {
    result = Some(if is_root {
      distribution_add_neg_log_weight(distribution, |time| model.root_contribution(time, n_children))?
    } else {
      distribution_add_neg_log_weight(distribution, |time| model.internal_contribution(time, n_children))?
    });
  }

  // Lift the input date constraint into the time distribution. It is an independent factor of the
  // node's posterior, so it multiplies whatever the children have to say; for a leaf, which has no
  // children, it is the whole message. Doing this here rather than once at load time is what keeps
  // the input recoverable: the forward pass refines the time distribution of a node whose date is
  // uncertain in place, and sending that refined distribution back to the parent on the next round
  // would count the parent's own message toward the node a second time.
  if let Some(constraint) = slot.node.date_constraint().clone() {
    result = Some(match result {
      Some(dist) => distribution_multiplication(&dist, &constraint)?,
      None => constraint.as_ref().clone(),
    });
  }

  if let Some(dist) = result.as_ref() {
    // Re-window once at the end: size the combined posterior's grid by probability mass (design D3),
    // which peak-normalizes as its first step (subsuming the shift-only normalize this replaces).
    // Every downstream consumer (likely_time, quantile, hpd_region, and the outgoing convolution via
    // to_plain_normalized) is shift-invariant, so the peak offset removed here has no effect on
    // inferred times or likelihoods.
    slot
      .node
      .set_time_distribution(Some(Arc::new(rewindow_to_mass(dist, EPS, GRID_POINTS)?)));
  }

  // A leaf's coalescent factor belongs only to the temporary message convolved toward the parent,
  // never to the distribution stored on the node.
  let outgoing_distribution = if is_leaf {
    let distribution = slot.node.time_distribution();
    match (coalescent_model, distribution) {
      (Some(model), Some(distribution)) => Some(Arc::new(distribution_add_neg_log_weight(
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
    // Tail policy for the backward message (kb/decisions/distribution-tails-and-arithmetic.md).
    // The parent could be arbitrarily far in the past, so the left side is soft: a fitted log-linear
    // Linear tail that decays with finite mass. A flat Constant tail is non-integrable and corrupts
    // the quantile and HPD integrals, so it is retired here. The child's sampling date is a hard
    // upper bound on the parent's age, so the right tail stays Hard.
    let message = distribution_convolution(node_time_dist.as_ref(), &negated_branch_dist)?;
    let left_tail = fit_message_soft_tail(&message, Side::Left)?;
    let parent_message = message
      .with_left_extrap(left_tail)?
      .with_right_extrap(BoundaryBehavior::Hard)?;
    edge.set_msg_to_parent(Some(Arc::new(parent_message)));
  }

  Ok(())
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
  // before the sum (design D3). The mass re-window to `GRID_POINTS` happens once on the combined
  // posterior after the fold.
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
/// the soft laws (`Constant`, `Linear`) continue past the edge and do not.
fn is_restricting(tail: BoundaryBehavior) -> bool {
  tail.is_hard()
}

/// Derive the concrete tail policy of a summed message on one side from the operands' tail classes.
///
/// Mirrors the multiplication tail rule (`treetime_distribution::distribution_ops::multiply`) as a
/// class lattice `Error > Hard > Linear > Constant`: the product is soft only when every operand is
/// soft, a hard bound dominates any soft one, and an undeclared (`Error`) bound dominates all. The
/// backward messages carry no fitted law, so the hard and undeclared classes need none here
/// (`Hard`, `Error`); the flat class stays `Constant`. Only the soft log-linear class needs a law,
/// which is fit from the summed data on `side` -- a grid too degenerate to fit is an error rather
/// than a silent flat fallback.
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
      BoundaryBehavior::Constant => {},
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
  Ok(BoundaryBehavior::Constant)
}
