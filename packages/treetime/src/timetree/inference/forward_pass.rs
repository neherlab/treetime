use crate::partition::indexed_pass::{IndexedPassDependencies, IndexedPassSlot, with_indexed_graph_payloads};
use crate::payload::traits::{TimetreeEdge, TimetreeNode};
use crate::timetree::inference::runner::{EPS, GRID_POINTS};
use eyre::Report;
use log::{Level, debug, log_enabled, warn};
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use treetime_distribution::Distribution;
use treetime_distribution::NegLog;
use treetime_distribution::convolve_across_edge;
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

  // Once per pass, not per node: a broken clock or topology makes a whole subtree disagree at once.
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

/// Refine the node's time distribution from its parent, then commit its point-estimate time.
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
  commit_node_time(graph, dependencies, slot);
  Ok(())
}

/// What [`refine_distribution_from_parent`] did to a node.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Refinement {
  /// Refined from the parent, or nothing to refine (root, exact date, or a missing operand).
  Done,
  /// The parent message ruled out the given date; the date was kept unrefined.
  ContradictedGivenDate,
}

/// Refine the node's time distribution with the message coming down from its parent.
///
/// The message -- the rest of the tree's opinion, carried across the branch -- is multiplied into the
/// node's own subtree evidence. When the two have disjoint support the product is empty; the given
/// date is then kept rather than refined away, and the caller told so it can report the disagreement.
fn refine_distribution_from_parent<N, E>(
  dependencies: &IndexedPassDependencies<N, E>,
  slot: &mut IndexedPassSlot<N, E>,
) -> Result<Refinement, Report>
where
  N: GraphNode + Named + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
{
  // Nothing to refine:
  // - root -- no parent message
  // - exact date -- observed, not inferred
  // A leaf with an uncertain date is refined like any inferred node.
  let Some(parent_key) = slot.parent_key else {
    return Ok(Refinement::Done);
  };
  if has_exact_date(&slot.node) {
    return Ok(Refinement::Done);
  }

  let parent = dependencies.node(parent_key);
  let (_, edge) = slot
    .parent_edge
    .as_ref()
    .expect("Non-root indexed node must own its parent edge");

  let (Some(parent_time_dist), Some(branch_dist)) = (parent.time_distribution(), edge.branch_length_distribution())
  else {
    return Ok(Refinement::Done);
  };

  // No distribution of its own: take the parent message across the branch as is (nothing to divide
  // out, multiply back, or contradict).
  let Some(subtree_dist) = slot.node.time_distribution() else {
    let dist_from_parent = convolve_across_edge(parent_time_dist, branch_dist, Side::Right, EPS, GRID_POINTS)?;
    log_refinement(&slot.node, parent_time_dist, &dist_from_parent);
    slot.node.set_time_distribution(Some(Arc::new(dist_from_parent)));
    return Ok(Refinement::Done);
  };

  let dist_from_parent = message_from_parent(
    parent_time_dist,
    edge.msg_to_parent().as_deref(),
    subtree_dist,
    branch_dist,
  )?;

  // Peak-normalize and store. The multiply is pointwise and does not resize the grid, so no further
  // mass sizing is needed.
  let combined = distribution_multiplication(&dist_from_parent, subtree_dist)?.normalize();
  log_refinement(&slot.node, parent_time_dist, &combined);

  // Empty product: the parent message and the given date have disjoint support. Refining onto it
  // would leave the node undated -- worse than the input date -- so keep the given date and report
  // it. A node with no given date has nothing to fall back on and stays undated.
  if combined.likely_time().is_none() && slot.node.date_constraint().is_some() {
    log_kept_given_date(&slot.node, &dist_from_parent);
    return Ok(Refinement::ContradictedGivenDate);
  }
  slot.node.set_time_distribution(Some(Arc::new(combined)));
  Ok(Refinement::Done)
}

/// The message the parent sends down to this node across the branch.
///
/// This node's own contribution was folded into the parent's posterior during the backward pass, so
/// it is divided back out first -- the cavity distribution -- then convolved across the branch. With
/// no message on record there is nothing to divide out and the parent posterior is used as is.
fn message_from_parent(
  parent_time_dist: &Distribution<NegLog>,
  msg_to_parent: Option<&Distribution<NegLog>>,
  subtree_dist: &Distribution<NegLog>,
  branch_dist: &Distribution<NegLog>,
) -> Result<Distribution<NegLog>, Report> {
  let parent_except_subtree = match msg_to_parent {
    Some(msg_to_parent) => distribution_division(parent_time_dist, msg_to_parent)?,
    None => parent_time_dist.clone(),
  };

  // TODO: remove or review (richard) -- analogous to the backward-pass regridding after message
  // combination.
  let parent_except_subtree = restrict_to_reachable(parent_except_subtree, subtree_dist, branch_dist)?;

  // Tail policy (kb/decisions/distribution-tails-and-arithmetic.md): parent time is a hard lower
  // bound (left Hard); the node's time is unbounded forward, so the right tail is soft and fitted.
  convolve_across_edge(&parent_except_subtree, branch_dist, Side::Right, EPS, GRID_POINTS)
}

/// Cut the parent's distribution down to the times from which the branch can reach `target`.
///
/// The refinement multiplies the parent message by `target`, which is zero outside its own support,
/// so parent times the branch cannot carry into that support contribute nothing: dropping them
/// leaves the result unchanged. This keeps a narrow date range affordable against a parent grid that
/// can span centuries. The parent is returned untouched when:
/// - it is not gridded, or either support is unbounded
/// - the window would save too little to pay for the resampling
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
  //TODO: we should be sampling into the  branch's tail if necessary.

  let dx = parent_fn.dx();
  let span = parent_fn.x_max() - parent_fn.x_min();
  // Nothing to save on a grid already a few points wide, or if the geometry is not finite.
  if !(dx.is_finite() && span.is_finite()) || dx <= 0.0 || span <= MIN_STEPS * dx {
    return Ok(parent);
  }

  // Window around the reachable times, clamped to the parent's grid and widened to a few steps when
  // the date range is narrower than the parent's resolution.
  let reach_min = (t_min - d_max).max(parent_fn.x_min());
  let reach_max = (t_max - d_min).min(parent_fn.x_max());
  let centre = f64::midpoint(reach_min, reach_max);
  let half_width = f64::max((reach_max - reach_min) / 2.0, MIN_STEPS * dx / 2.0);
  let window_min = (centre - half_width).max(parent_fn.x_min());
  let window_max = (centre + half_width).min(parent_fn.x_max());

  // Snapped outwards onto the parent's own grid points, so the window holds values it already has
  // and the convolution lands where it would have with nothing cut. Otherwise the cut shifts node
  // times by a fraction of a grid step -- no worse, but gratuitously different.
  let steps_from_start = |t: f64| (t - parent_fn.x_min()) / dx;
  let window_min = parent_fn.x_min() + steps_from_start(window_min).floor().max(0.0) * dx;
  let window_max = (parent_fn.x_min() + steps_from_start(window_max).ceil() * dx).min(parent_fn.x_max());

  // Too narrow to grid, or saving too little to pay for the resampling: leave the parent as is.
  let width = window_max - window_min;
  if !width.is_finite() || width < 2.0 * dx || width > WORTHWHILE * span {
    return Ok(parent);
  }

  // Clamped resample holds the boundary value where rounding overshoots the grid's end points; the
  // parent's own tails are restored afterwards. Same treatment as `DistributionFunction::resample_dx`.
  let restricted = parent_fn
    .resample_range_dx_clamped((window_min, window_max), parent_fn.dx())?
    .with_left_extrap(parent_fn.left_extrap())?
    .with_right_extrap(parent_fn.right_extrap())?;
  Ok(Distribution::Function(restricted))
}

/// Commit the node's point-estimate time from the peak of its refined time distribution.
///
/// Projected to be no earlier than the parent's committed time -- except for an exact date, which is
/// kept so a clock conflict stays visible to `commit_clock_branch_lengths`. An empty distribution
/// gets no time, and warns when the node was dateable.
fn commit_node_time<N, E, D>(
  graph: &Graph<N, E, D>,
  dependencies: &IndexedPassDependencies<N, E>,
  slot: &mut IndexedPassSlot<N, E>,
) where
  N: GraphNode + Named + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
  D: Send + Sync,
{
  // Project the inferred point estimate onto the committed parent time (an exact date keeps its
  // own). This adjusts the estimate without recomputing the posterior; the statistical contract is
  // open in kb/issues/M-timetree-marginal-node-times-can-violate-topology.md.
  let parent_time = (!has_exact_date(&slot.node))
    .then(|| parent_time(dependencies, slot))
    .flatten();

  // An empty distribution yields no time. Under NegLog (ordinate -ln p) even tiny posteriors survive
  // exactly, so empty means genuinely disjoint hard domains -- the subtree disagrees with the rest of
  // the tree. Surface that rather than drop the date silently. An undated leaf under an undated parent
  // has nothing to infer from and was already reported when its date was found missing.
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
pub(super) fn set_likely_time(node: &mut impl TimetreeNode, parent_time: Option<f64>) -> Option<f64> {
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
