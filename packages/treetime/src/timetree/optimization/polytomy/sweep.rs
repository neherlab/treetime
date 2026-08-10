//! Stochastic coalescent-with-mutations sweep over the children of one polytomy.
//!
//! Pure simulation: no graph, no partitions, no I/O. Takes a summary of the children and
//! returns a plan of mergers for [`super::apply`] to realise. Keeping it separate makes the
//! part that is easy to get wrong -- the competing-clocks bookkeeping -- testable with a
//! seeded RNG and no tree.
//!
//! # Model
//!
//! The children of a polytomy are lineages whose branches carry a known number of
//! substitutions. Sweeping backwards in time from the most recent child toward the parent,
//! two kinds of event compete:
//!
//! - **Mutation**: one of the substitutions on a branch is placed at the current time.
//! - **Coalescence**: two lineages merge under a new node.
//!
//! A lineage may only coalesce once every substitution on its branch has been placed, and only
//! once the sweep has reached it (a lineage is not available above its own node time). Writing
//! $\mu$ for the whole-alignment mutation rate, $\kappa(t)$ for the per-branch coalescent
//! merger rate, $M$ for the substitutions still unplaced across live lineages, and $R$ for the
//! set of live mutation-free lineages:
//!
//! $$R_{\text{mut}} = \mu M, \qquad R_{\text{coal}} = \max(0, \lvert R\rvert - 1)\,\kappa(t)$$
//!
//! The waiting time is $\operatorname{Exp}(R_{\text{mut}} + R_{\text{coal}})$; the event is a
//! mutation with probability $R_{\text{mut}}/(R_{\text{mut}} + R_{\text{coal}})$; the mutating
//! branch is drawn $\propto m_b$ and the coalescing pair uniformly from $R$.
//!
//! # Divergence from v0
//!
//! v0's `generate_subtree` adds $\mu$ into both rate channels, which leaves its
//! branch-selection weights inconsistent with the rate that gated them; see
//! [kb/v0-errata/timetree-stochastic-resolve-rate-selection-mismatch.md]. Two further v0
//! defects are avoided here: events drawn past the parent bound are discarded rather than
//! committed (which in v0 yields negative branch lengths), and a draw that crosses a lineage
//! arrival resumes at the arrival time rather than at the drawn time (which in v0 skips an
//! interval entirely). See the sibling errata entries.
//!
//! # Time convention
//!
//! Calendar time throughout: the parent is older and has the *smaller* value, so the sweep
//! runs from `max(child.time)` **downwards** to `t_stop`. v0 works in `time_before_present`
//! and runs upwards; every comparison here is inverted relative to it.

use eyre::Report;
use rand::Rng;
use rand_distr::{Distribution as _, Exp};
use std::cmp::Ordering;
use std::collections::VecDeque;
use treetime_utils::make_internal_error;

/// One child branch of the polytomy, as seen by the sweep.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Lineage {
  /// Calendar time of the child node. Larger is more recent.
  pub time: f64,
  /// Substitutions mapped to the branch above this node.
  pub mutations: u32,
}

/// One coalescence: `left` and `right` become the children of a new node at `time`.
///
/// Lineage ids below the child count index the `children` slice passed to
/// [`simulate_subtree`]. Id `children.len() + j` refers to the node created by the `j`-th
/// entry of [`SubtreePlan::mergers`], so a merger may only reference earlier mergers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Merger {
  pub time: f64,
  pub left: usize,
  pub right: usize,
}

/// What the sweep produced.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SubtreePlan {
  /// Mergers in creation order, oldest last.
  pub mergers: Vec<Merger>,
  /// Lineages still unmerged when the sweep stopped. These stay children of the polytomy
  /// parent, so a plan with fewer than two mergers short of full resolution leaves a residual
  /// multifurcation -- the expected outcome when the time window runs out.
  pub roots: Vec<usize>,
}

impl SubtreePlan {
  /// A plan that changes nothing: every child stays a direct child of the parent.
  fn unresolved(n_children: usize) -> Self {
    Self {
      mergers: Vec::new(),
      roots: (0..n_children).collect(),
    }
  }
}

/// A lineage the sweep is tracking, live or pending.
#[derive(Clone, Copy, Debug)]
struct Tracked {
  id: usize,
  time: f64,
  mutations: u32,
}

/// Run the sweep over one polytomy's children.
///
/// `t_stop` is the parent's calendar time: the sweep may not place events at or before it,
/// because doing so would put a node older than its own parent. Returns early without
/// mergers when there is nothing to resolve (fewer than three children) or no time to
/// resolve it in (`max(child.time) <= t_stop`).
///
/// `merger_rate` is the per-branch coalescent merger rate $\kappa(t)$ at a calendar time. The
/// sweep scales it by its own $\lvert R\rvert - 1$ to obtain the total; where $\kappa$ comes
/// from is the caller's choice, so no coalescent model reaches this module.
pub fn simulate_subtree(
  children: &[Lineage],
  t_stop: f64,
  mutation_rate: f64,
  merger_rate: &dyn Fn(f64) -> Result<f64, Report>,
  rng: &mut dyn rand::RngCore,
) -> Result<SubtreePlan, Report> {
  let n_children = children.len();
  if n_children < 3 {
    return Ok(SubtreePlan::unresolved(n_children));
  }

  // Most recent first. Ties break by id so the sweep is a deterministic function of the RNG
  // stream rather than of the input ordering.
  let mut ordered: Vec<Tracked> = children
    .iter()
    .enumerate()
    .map(|(id, child)| Tracked {
      id,
      time: child.time,
      mutations: child.mutations,
    })
    .collect();
  ordered.sort_by(|a, b| {
    b.time
      .partial_cmp(&a.time)
      .unwrap_or(Ordering::Equal)
      .then(a.id.cmp(&b.id))
  });

  let t_start = ordered[0].time;
  if !(t_start > t_stop) {
    return Ok(SubtreePlan::unresolved(n_children));
  }

  let mut to_come: VecDeque<Tracked> = ordered.into_iter().collect();
  let mut alive: Vec<Tracked> = Vec::with_capacity(n_children);
  let mut t = t_start;
  admit_arrivals(&mut alive, &mut to_come, t);

  let mut mergers: Vec<Merger> = Vec::new();

  while alive.len() + to_come.len() > 2 && t > t_stop {
    let n_ready = alive.iter().filter(|lineage| lineage.mutations == 0).count();
    let total_mutations: u64 = alive.iter().map(|lineage| u64::from(lineage.mutations)).sum();

    let kappa = merger_rate(t)?;
    let rate_mut = mutation_rate * total_mutations as f64;
    let rate_coal = n_ready.saturating_sub(1) as f64 * kappa;
    let rate_total = rate_mut + rate_coal;

    // No event can occur in the current configuration: a single ready lineage cannot merge
    // and mutation-free lineages cannot mutate. Advance to the next arrival, which changes
    // the configuration, or stop if there is none. v0 instead relies on `Exp(0)` returning
    // infinity and the loop condition catching it on the next pass.
    if !(rate_total > 0.0) {
      let Some(next) = to_come.front().copied() else {
        break;
      };
      t = next.time;
      admit_arrivals(&mut alive, &mut to_come, t);
      continue;
    }

    let dt = Exp::new(rate_total)?.sample(rng);
    let t_event = t - dt;

    // The rate is only valid until the next lineage becomes available. If the draw crosses
    // that boundary no event occurred before it, so resume there and redraw under the new
    // rate -- the exponential is memoryless, so discarding the draw is unbiased. Resuming at
    // `t_event` instead (as v0 does) would skip the interval between the arrival and the
    // drawn time.
    if let Some(next) = to_come.front().copied() {
      if next.time >= t_event {
        t = next.time;
        admit_arrivals(&mut alive, &mut to_come, t);
        continue;
      }
    }

    // The draw ran past the parent: no event occurred within the available window. Stop
    // without committing it. v0 commits the event first and lets the loop condition notice
    // afterwards, which places a node older than its parent and yields a negative branch
    // length.
    if t_event <= t_stop {
      break;
    }
    t = t_event;

    if rng.gen_bool((rate_mut / rate_total).clamp(0.0, 1.0)) {
      place_mutation(&mut alive, total_mutations, rng)?;
    } else {
      let merger = coalesce_pair(&mut alive, t, n_children + mergers.len(), rng)?;
      mergers.push(merger);
    }
  }

  let roots = alive.iter().chain(to_come.iter()).map(|lineage| lineage.id).collect();

  Ok(SubtreePlan { mergers, roots })
}

/// Move every lineage whose node is at or more recent than `t` into the live set.
///
/// `to_come` is ordered most-recent-first, so this is a prefix pop.
fn admit_arrivals(alive: &mut Vec<Tracked>, to_come: &mut VecDeque<Tracked>, t: f64) {
  while to_come.front().is_some_and(|next| next.time >= t) {
    let arrived = to_come.pop_front().expect("front was just inspected");
    alive.push(arrived);
  }
}

/// Consume one unplaced substitution, choosing the branch in proportion to how many it has
/// left.
///
/// The weights sum to `total_mutations`, which is exactly the quantity that produced
/// `rate_mut`, so the scan cannot run off the end. v0's equivalent can, because its mutation
/// rate carries a term with no branch behind it; it then silently mutates an arbitrary
/// branch. An unreachable state here is reported instead.
fn place_mutation(alive: &mut [Tracked], total_mutations: u64, rng: &mut dyn rand::RngCore) -> Result<(), Report> {
  let mut remaining = rng.gen_range(0.0..total_mutations as f64);
  for lineage in alive.iter_mut() {
    remaining -= f64::from(lineage.mutations);
    if remaining < 0.0 {
      lineage.mutations -= 1;
      return Ok(());
    }
  }
  make_internal_error!(
    "Polytomy sweep drew a mutation event but no branch carried it: {total_mutations} mutations were expected across {} live branches",
    alive.len()
  )
}

/// Merge two uniformly chosen mutation-free lineages into a new lineage at `time`.
fn coalesce_pair(
  alive: &mut Vec<Tracked>,
  time: f64,
  new_id: usize,
  rng: &mut dyn rand::RngCore,
) -> Result<Merger, Report> {
  let ready: Vec<usize> = alive
    .iter()
    .enumerate()
    .filter(|(_, lineage)| lineage.mutations == 0)
    .map(|(index, _)| index)
    .collect();

  if ready.len() < 2 {
    // Unreachable: `rate_coal` is zero below two ready lineages, so a coalescence cannot be
    // selected. Reported rather than asserted so a broken rate surfaces as an error.
    return make_internal_error!(
      "Polytomy sweep drew a coalescence with {} lineages ready to merge",
      ready.len()
    );
  }

  let first = rng.gen_range(0..ready.len());
  let mut second = rng.gen_range(0..ready.len() - 1);
  if second >= first {
    second += 1;
  }

  // Remove the higher index first so the lower one stays valid.
  let (lo, hi) = if ready[first] < ready[second] {
    (ready[first], ready[second])
  } else {
    (ready[second], ready[first])
  };
  let right = alive.remove(hi);
  let left = alive.remove(lo);

  alive.push(Tracked {
    id: new_id,
    time,
    mutations: 0,
  });

  // Order the pair by id so a plan is comparable across runs that drew the same pair in the
  // opposite order.
  let (left_id, right_id) = if left.id <= right.id {
    (left.id, right.id)
  } else {
    (right.id, left.id)
  };

  Ok(Merger {
    time,
    left: left_id,
    right: right_id,
  })
}
