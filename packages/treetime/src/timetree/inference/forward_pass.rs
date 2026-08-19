use crate::partition::indexed_pass::{IndexedPassDependencies, IndexedPassSlot, with_indexed_graph_payloads};
use crate::payload::traits::{TimetreeEdge, TimetreeNode};
use crate::timetree::inference::tail_fit::fit_message_soft_tail;
use eyre::Report;
use log::{Level, debug, log_enabled, warn};
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use treetime_distribution::BoundaryBehavior;
use treetime_distribution::Distribution;
use treetime_distribution::NegLog;
use treetime_distribution::distribution_convolution;
use treetime_distribution::distribution_division;
use treetime_distribution::distribution_multiplication;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_grid::Side;

pub fn propagate_distributions_forward<N, E, D>(graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode + Named + TimetreeNode + Default,
  E: GraphEdge + TimetreeEdge + Default,
  D: Send + Sync,
{
  let contradicted = AtomicUsize::new(0);
  with_indexed_graph_payloads(graph, |pass| {
    pass.try_for_each_forward(|dependencies, slot| {
      propagate_distributions_forward_slot(graph, dependencies, slot, &contradicted)
    })
  })?;

  // Reported once for the pass rather than per node: on a tree whose clock or topology is off,
  // every dated node in a subtree disagrees at once, and a few thousand identical warnings say no
  // more than their count does. The per-node detail is in the debug trace.
  let contradicted = contradicted.load(Ordering::Relaxed);
  if contradicted > 0 {
    warn!(
      "Timetree forward pass: {contradicted} node(s) carry a date that the rest of the tree gives \
       no probability at all, so their posterior came out empty and each kept the date it was \
       given, unrefined. The usual cause is a sequence whose divergence implies a date far from \
       the one it is stamped with, which the clock filter reports separately. Run with \
       `--verbosity=debug` to see which nodes and where the tree puts each of them."
    );
  }

  Ok(())
}

fn propagate_distributions_forward_slot<N, E, D>(
  graph: &Graph<N, E, D>,
  dependencies: &IndexedPassDependencies<N, E>,
  slot: &mut IndexedPassSlot<N, E>,
  contradicted: &AtomicUsize,
) -> Result<(), Report>
where
  N: GraphNode + Named + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
  D: Send + Sync,
{
  if refine_distribution_from_parent(dependencies, slot)? == Refinement::ContradictedGivenDate {
    contradicted.fetch_add(1, Ordering::Relaxed);
  }

  let has_exact_date = has_exact_date(&slot.node);

  // Independent marginal modes can invert even though the forward pass commits each parent before
  // visiting its children. Current behavior projects the point estimate of a node whose time was
  // inferred to the committed parent time; a node given an exact date keeps it, so that a date
  // conflicting with the fitted clock stays visible to `commit_clock_branch_lengths` rather than
  // being clamped away. This changes the point estimate without recomputing the posterior. The
  // statistical contract remains open in
  // kb/issues/M-timetree-marginal-node-times-can-violate-topology.md.
  let parent_time = (!has_exact_date).then(|| parent_time(dependencies, slot)).flatten();

  // A node whose time distribution is empty (or degenerate) yields no likely time, so no date is
  // assigned. Under exact arithmetic that cannot happen: the messages being combined all went into
  // the parent's own posterior, so some time always carries probability. Under NegLog the ordinate
  // is -ln(probability), so even astronomically small posteriors are represented exactly and never
  // collapse to zero; an empty result therefore means the messages meeting here have genuinely
  // disjoint support (their hard domains do not overlap). Surface it instead of dropping the date
  // silently: it means the subtree below this node says something the rest of the tree gives no
  // weight to at all. An undated leaf under an undated parent has nothing to infer from and is not
  // worth reporting: it was already accounted for when its date was found missing.
  let is_dateable = !graph.is_leaf(slot.key) || slot.node.date_constraint().is_some();
  if set_likely_time(&mut slot.node, parent_time).is_none() && is_dateable {
    let name = slot.node.name();
    let name = name.as_ref().map_or("<unnamed>", |name| name.as_ref());
    warn!(
      "Timetree forward pass: node '{name}' has an empty time distribution; no date was assigned. \
       The messages meeting at this node leave no time with any probability: the dates below it \
       and the times the rest of the tree implies have disjoint support."
    );
  }
  Ok(())
}

/// Whether the node was given an exact date, which pins it to a single time.
///
/// Anything else -- an uncertain or ranged date, or no date at all -- leaves the node's time to be
/// inferred, and so to be refined by the message coming down from its parent.
fn has_exact_date(node: &impl TimetreeNode) -> bool {
  node.date_constraint().as_ref().is_some_and(|dist| dist.is_point())
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

/// Assign the node's committed time from the peak of its time distribution, clamped to be no
/// earlier than the parent's committed time. Returns the assigned time, or `None` when the time
/// distribution is empty or degenerate and no time could be determined (the node is left undated).
fn set_likely_time(node: &mut impl TimetreeNode, parent_time: Option<f64>) -> Option<f64> {
  let time = node
    .time_distribution()
    .as_ref()
    .and_then(|time_dist| time_dist.likely_time())?;

  let time = parent_time.map_or(time, |parent_time| time.max(parent_time));
  node.set_time(Some(time));
  Some(time)
}

/// Trace what the refinement of one node did to its grid, for `RUST_LOG=debug`.
///
/// The grid a node ends up on is set by the operands it was built from, so a grid that grows or a
/// support that drifts is read off the pair: `parent` is the distribution the message came from and
/// `refined` is what the node now carries.
fn log_refinement<N: TimetreeNode + Named>(node: &N, parent: &Distribution<NegLog>, refined: &Distribution<NegLog>) {
  if !log_enabled!(Level::Debug) {
    return;
  }
  let name = node.name();
  let name = name.as_ref().map_or("<unnamed>", |name| name.as_ref());
  debug!(
    "Timetree forward pass: node '{name}': parent {} -> refined {}",
    describe_grid(parent),
    describe_grid(refined)
  );
}

/// Name the node whose given date the rest of the tree contradicts, and what the tree implied.
fn log_kept_given_date<N: TimetreeNode + Named>(node: &N, dist_from_parent: &Distribution<NegLog>) {
  if !log_enabled!(Level::Debug) {
    return;
  }
  let name = node.name();
  let name = name.as_ref().map_or("<unnamed>", |name| name.as_ref());
  let given = node
    .date_constraint()
    .as_ref()
    .map_or_else(|| "none".to_owned(), |constraint| describe_grid(constraint.as_ref()));
  debug!(
    "Timetree forward pass: node '{name}' keeps the date it was given, {given}: the rest of the \
     tree puts it at {}, which leaves no probability on that date",
    describe_grid(dist_from_parent)
  );
}

/// Grid size and support of a distribution: `n=<points>, [<t_min>, <t_max>]`.
///
/// The point count is what the distribution would be evaluated on: one for a point date, two for a
/// range or a formula's end points, and the full grid for a function.
fn describe_grid(dist: &Distribution<NegLog>) -> String {
  let n_points = dist.t().len();
  match dist.time_bounds() {
    Some((t_min, t_max)) => format!("n={n_points}, [{t_min:.6}, {t_max:.6}]"),
    None => "empty".to_owned(),
  }
}

/// What [`refine_distribution_from_parent`] did to a node.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Refinement {
  /// The node's time distribution now carries the message from its parent, or there was nothing to
  /// refine it with: the node is the root, its date is exact, or one of the distributions the
  /// refinement needs is missing.
  Done,
  /// The message from the parent put the node outside the date it was given. The date was kept as
  /// it is, so the node stays dated where the input says.
  ContradictedGivenDate,
}

fn refine_distribution_from_parent<N, E>(
  dependencies: &IndexedPassDependencies<N, E>,
  slot: &mut IndexedPassSlot<N, E>,
) -> Result<Refinement, Report>
where
  N: GraphNode + Named + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
{
  if slot.parent_key.is_none() {
    return Ok(Refinement::Done);
  }

  // An exactly dated node has nothing to refine: its time is observed, not inferred. Leaves whose
  // date is uncertain or missing are refined like any other node -- the message from the parent is
  // the only thing that can narrow an ambiguous date down.
  if has_exact_date(&slot.node) {
    return Ok(Refinement::Done);
  }

  if let Some(parent_key) = slot.parent_key {
    let parent = dependencies.node(parent_key);
    let (_, edge) = slot
      .parent_edge
      .as_ref()
      .expect("Non-root indexed node must own its parent edge");

    // case where parent, branch, and constraint/distribution on the node all exist.
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
      let parent_except_subtree = restrict_to_reachable(parent_except_subtree, subtree_dist, branch_dist)?;

      // Tail policy for the forward message (kb/decisions/distribution-tails-and-arithmetic.md).
      // The parent's time is a hard lower bound (left tail Hard). There is no upper bound from the
      // parent side on how far in the future the node could be, so the right side is soft: a fitted
      // log-linear Linear tail with finite mass, so the quantile and HPD integrals stay well-defined.
      let forward_message = distribution_convolution(&parent_except_subtree, branch_dist)?;
      let right_tail = fit_message_soft_tail(&forward_message, Side::Right)?;
      let dist_from_parent = forward_message
        .with_left_extrap(BoundaryBehavior::Hard)?
        .with_right_extrap(right_tail)?;
      // Peak-normalize the refined posterior and store it as-is. Multiplying the parent message by
      // the subtree evidence is a pointwise operation that does not resize the grid, so the stored
      // posterior needs no mass sizing: that belongs to the edge crossing (the convolution that
      // produced the parent message), not to this product.
      let combined = distribution_multiplication(&dist_from_parent, subtree_dist)?.normalize();
      log_refinement(&slot.node, parent_time_dist, &combined);

      // The product is empty when the message from the parent leaves no probability on the date
      // the node was given -- in exact arithmetic never, and under NegLog only when the parent
      // message and the subtree distribution have disjoint support. Refining onto that would leave
      // the node undated, which is strictly worse than the date the input carries: keep the given
      // date and let the caller report the disagreement. A node with no given date has nothing to
      // fall back on and is still left undated, which is the pre-existing contract.
      if combined.likely_time().is_none() && slot.node.date_constraint().is_some() {
        log_kept_given_date(&slot.node, &dist_from_parent);
        return Ok(Refinement::ContradictedGivenDate);
      }
      slot.node.set_time_distribution(Some(Arc::new(combined)));
    } else if let (Some(parent_time_dist), Some(branch_dist)) =
      (parent.time_distribution(), edge.branch_length_distribution())
    {
      // Same forward tail policy as the two-parent-message branch above: Hard left, fitted Linear
      // right. Fit after normalize(), a pure ordinate shift that preserves the grid the fit reads
      // and leaves the slope unchanged (the soft-tail slope is shift-invariant).
      let forward_message = distribution_convolution(parent_time_dist, branch_dist)?.normalize();
      let right_tail = fit_message_soft_tail(&forward_message, Side::Right)?;
      let dist_from_parent = forward_message
        .with_left_extrap(BoundaryBehavior::Hard)?
        .with_right_extrap(right_tail)?;
      log_refinement(&slot.node, parent_time_dist, &dist_from_parent);
      slot.node.set_time_distribution(Some(Arc::new(dist_from_parent)));
    }
  }

  Ok(Refinement::Done)
}

/// Cut the parent's distribution down to the times from which the branch can reach `target`.
///
/// The refinement multiplies the message from the parent by `target`, which is zero outside its own
/// support, so a parent time the branch cannot carry into that support contributes nothing to the
/// result. Dropping those times before the convolution leaves the result unchanged where it is not
/// zero.
///
/// This is what makes refining a node with a narrow date range affordable. The convolution costs
/// the product of the two grid sizes, and a posterior deep in a large tree can carry tens of
/// thousands of points spanning centuries, where a date range spans weeks. The parent is returned
/// untouched when it is not gridded, when either support is unbounded, or when the window would cut
/// little enough that the resampling would cost more than it saves.
fn restrict_to_reachable(
  parent: Distribution<NegLog>,
  target: &Distribution<NegLog>,
  branch: &Distribution<NegLog>,
) -> Result<Distribution<NegLog>, Report> {
  /// Fraction of the parent's support the window has to get under to be worth resampling.
  const WORTHWHILE: f64 = 0.9;
  /// Grid steps the window is widened to when the node's support is narrower than the parent's
  /// resolution. Below two steps there is no grid to resample onto at all.
  const MIN_STEPS: f64 = 4.0;

  let (Distribution::Function(parent_fn), Some((t_min, t_max)), Some((d_min, d_max))) =
    (&parent, target.time_bounds(), branch.time_bounds())
  else {
    return Ok(parent);
  };

  let dx = parent_fn.dx();
  let span = parent_fn.x_max() - parent_fn.x_min();
  // Nothing to save on a grid that is already a handful of points wide, and nothing to compute a
  // window from if either quantity is not a number.
  if !(dx.is_finite() && span.is_finite()) || dx <= 0.0 || span <= MIN_STEPS * dx {
    return Ok(parent);
  }

  // Times outside the parent's own grid are not ours to invent, so the window is clamped to it.
  // A date range can be narrower than the parent's resolution, which leaves too few grid points to
  // resample onto, so the window is held to a few steps around the times that can reach the node.
  let reach_min = (t_min - d_max).max(parent_fn.x_min());
  let reach_max = (t_max - d_min).min(parent_fn.x_max());
  let centre = f64::midpoint(reach_min, reach_max);
  let half_width = f64::max((reach_max - reach_min) / 2.0, MIN_STEPS * dx / 2.0);
  let window_min = (centre - half_width).max(parent_fn.x_min());
  let window_max = (centre + half_width).min(parent_fn.x_max());

  // Snapped outwards onto the parent's own grid points, so that the window holds the values it
  // already has rather than interpolations of them, and the convolution below lands on the grid it
  // would have landed on had nothing been cut. Without this the cut moves node times by a fraction
  // of a grid step -- no worse an estimate, but a gratuitously different one.
  let steps_from_start = |t: f64| (t - parent_fn.x_min()) / dx;
  let window_min = parent_fn.x_min() + steps_from_start(window_min).floor().max(0.0) * dx;
  let window_max = (parent_fn.x_min() + steps_from_start(window_max).ceil() * dx).min(parent_fn.x_max());

  // A window that does not meet the parent's grid, or that saves too little to pay for the
  // resampling, leaves the parent as it is. Where the two cannot be reconciled at all the caller
  // sees the same disagreement in the empty product, and reports it there.
  let width = window_max - window_min;
  if !width.is_finite() || width < 2.0 * dx || width > WORTHWHILE * span {
    return Ok(parent);
  }

  // Resampling evaluates the grid's own end points, which rounding can push a hair outside it, so
  // the clamped resample holds the boundary value at that overshoot and the parent's own tails are
  // put back afterwards. Same treatment as `DistributionFunction::resample_dx`.
  let restricted = parent_fn
    .resample_range_dx_clamped((window_min, window_max), parent_fn.dx())?
    .with_left_extrap(parent_fn.left_extrap())?
    .with_right_extrap(parent_fn.right_extrap())?;
  Ok(Distribution::Function(restricted))
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::payload::timetree::NodeTimetree;
  use crate::pretty_assert_ulps_eq;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_distribution::Distribution;

  #[test]
  fn test_forward_pass_set_likely_time_empty_distribution_returns_none() {
    let mut node = node_with_distribution(Some(Distribution::empty()));
    assert_eq!(None, set_likely_time(&mut node, None));
    assert_eq!(None, node.time());
  }

  #[test]
  fn test_forward_pass_set_likely_time_missing_distribution_returns_none() {
    let mut node = node_with_distribution(None);
    assert_eq!(None, set_likely_time(&mut node, None));
    assert_eq!(None, node.time());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::no_parent(      None,      5.0)]
  #[case::parent_earlier( Some(3.0), 5.0)]
  #[case::parent_clamps(  Some(8.0), 8.0)]
  #[trace]
  fn test_forward_pass_set_likely_time_uses_distribution_peak(
    #[case] parent_time: Option<f64>,
    #[case] expected: f64,
  ) {
    let mut node = node_with_distribution(Some(Distribution::point(5.0, 1.0)));
    let assigned = set_likely_time(&mut node, parent_time).expect("a time should be assigned");
    pretty_assert_ulps_eq!(assigned, expected, max_ulps = 4);
    let committed = node.time().expect("node time should be committed");
    pretty_assert_ulps_eq!(committed, expected, max_ulps = 4);
  }

  mod helpers {
    use super::*;

    pub(super) fn node_with_distribution(distribution: Option<Distribution<NegLog>>) -> NodeTimetree {
      NodeTimetree {
        time_distribution: distribution.map(Arc::new),
        ..NodeTimetree::default()
      }
    }
  }

  use helpers::*;
}
