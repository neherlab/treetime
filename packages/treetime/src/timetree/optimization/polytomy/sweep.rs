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
//! A unit-exponential hazard threshold is integrated across every interval where these rates
//! are constant. The event is a mutation with probability
//! $R_{\text{mut}}/(R_{\text{mut}} + R_{\text{coal}})$; the mutating branch is drawn
//! $\propto m_b$ and the coalescing pair uniformly from $R$.
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
use rand_distr::{Distribution as _, Exp1};
use std::collections::VecDeque;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
use treetime_utils::{make_error, make_internal_error};

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
  mutations: u32,
}

/// A lineage that becomes live when the sweep reaches `elapsed`.
#[derive(Clone, Copy, Debug)]
struct Pending {
  lineage: Tracked,
  elapsed: f64,
}

/// Run the sweep over one polytomy's children.
///
/// `t_stop` is the parent's calendar time: the sweep may not place events at or before it,
/// because doing so would put a node older than its own parent. Returns early without
/// mergers when there is nothing to resolve (fewer than three children) or no time to
/// resolve it in (`max(child.time) <= t_stop`).
///
/// `merger_rate` is the piecewise-constant per-branch coalescent merger rate $\kappa(t)$ in
/// calendar time. The sweep scales it by its own $\lvert R\rvert - 1$ to obtain the total and
/// integrates across every schedule breakpoint.
pub fn simulate_subtree(
  children: &[Lineage],
  t_stop: f64,
  mutation_rate: f64,
  merger_rate: &PiecewiseConstantFn,
  rng: &mut dyn rand::RngCore,
) -> Result<SubtreePlan, Report> {
  // Finite inputs make every time comparison a total-order comparison.
  validate_inputs(children, t_stop, mutation_rate)?;

  let n_children = children.len();
  if n_children < 3 {
    return Ok(SubtreePlan::unresolved(n_children));
  }

  // Most recent first. Ties break by id so the sweep is a deterministic function of the RNG
  // stream rather than of the input ordering.
  let mut ordered: Vec<(usize, &Lineage)> = children.iter().enumerate().collect();
  ordered.sort_by(|(left_id, left), (right_id, right)| right.time.total_cmp(&left.time).then(left_id.cmp(right_id)));

  let t_start = ordered[0].1.time;
  if t_start <= t_stop {
    return Ok(SubtreePlan::unresolved(n_children));
  }

  let mut to_come: VecDeque<Pending> = ordered
    .into_iter()
    .map(|(id, child)| Pending {
      lineage: Tracked {
        id,
        mutations: child.mutations,
      },
      elapsed: t_start - child.time,
    })
    .collect();
  let mut alive: Vec<Tracked> = Vec::with_capacity(n_children);
  let mut elapsed = 0.0;
  let stop_elapsed = t_start - t_stop;
  admit_arrivals(&mut alive, &mut to_come, elapsed);

  // Descending calendar breakpoints become ascending elapsed-time boundaries.
  let mut rate_boundaries: VecDeque<f64> = merger_rate
    .breakpoints()
    .iter()
    .rev()
    .copied()
    .filter(|&time| t_stop < time && time < t_start)
    .map(|time| t_start - time)
    .collect();

  let mut mergers: Vec<Merger> = Vec::new();

  'sweep: while alive.len() + to_come.len() > 2 && elapsed < stop_elapsed {
    // One threshold must survive every deterministic boundary before the next event.
    let mut hazard_left: f64 = Exp1.sample(rng);

    loop {
      let boundary = next_boundary(elapsed, stop_elapsed, &to_come, &rate_boundaries);
      let rate_time = t_start - f64::midpoint(elapsed, boundary);
      let (rate_mut, rate_total, total_mutations) =
        event_rates(&alive, mutation_rate, merger_rate.eval(rate_time), rate_time)?;
      let interval_hazard = rate_total * (boundary - elapsed);

      if hazard_left >= interval_hazard {
        if interval_hazard.is_finite() {
          hazard_left -= interval_hazard;
        }
        elapsed = boundary;
        if elapsed >= stop_elapsed {
          break 'sweep;
        }
        while rate_boundaries.front().is_some_and(|&next| next <= elapsed) {
          rate_boundaries.pop_front();
        }
        admit_arrivals(&mut alive, &mut to_come, elapsed);
        continue;
      }

      elapsed += hazard_left / rate_total;
      let event_time = t_start - elapsed;

      // Valid component rates make this ratio a probability without clamping.
      if rng.gen_bool(rate_mut / rate_total) {
        place_mutation(&mut alive, total_mutations, rng)?;
      } else {
        let merger = coalesce_pair(&mut alive, event_time, n_children + mergers.len(), rng)?;
        mergers.push(merger);
      }
      continue 'sweep;
    }
  }

  let roots = alive
    .iter()
    .map(|lineage| lineage.id)
    .chain(to_come.iter().map(|pending| pending.lineage.id))
    .collect();

  Ok(SubtreePlan { mergers, roots })
}

fn event_rates(alive: &[Tracked], mutation_rate: f64, kappa: f64, time: f64) -> Result<(f64, f64, u64), Report> {
  // A model rate must be valid before it enters event-rate arithmetic.
  if !kappa.is_finite() || kappa < 0.0 {
    return make_error!(
      "Polytomy merger rate must be finite and non-negative at calendar time {time:.6e}, got {kappa:.6e}"
    );
  }

  let n_ready = alive.iter().filter(|lineage| lineage.mutations == 0).count();
  let total_mutations = alive.iter().map(|lineage| u64::from(lineage.mutations)).sum();
  let rate_mut = mutation_rate * total_mutations as f64;
  let rate_coal = n_ready.saturating_sub(1) as f64 * kappa;
  // Finite component rates prevent overflow from becoming zero-hazard control flow.
  if !rate_mut.is_finite() || !rate_coal.is_finite() {
    return make_error!(
      "Polytomy event rates must be finite at calendar time {time:.6e}, got mutation rate {rate_mut:.6e} and merger rate {rate_coal:.6e}"
    );
  }
  let rate_total = rate_mut + rate_coal;
  if !rate_total.is_finite() {
    return make_error!("Polytomy total event rate must be finite at calendar time {time:.6e}, got {rate_total:.6e}");
  }

  Ok((rate_mut, rate_total, total_mutations))
}

fn next_boundary(elapsed: f64, stop_elapsed: f64, to_come: &VecDeque<Pending>, rate_boundaries: &VecDeque<f64>) -> f64 {
  let arrival = to_come.front().map_or(stop_elapsed, |pending| pending.elapsed);
  let rate_change = rate_boundaries.front().copied().unwrap_or(stop_elapsed);
  f64::min(stop_elapsed, f64::min(arrival, rate_change)).max(elapsed)
}

fn validate_inputs(children: &[Lineage], t_stop: f64, mutation_rate: f64) -> Result<(), Report> {
  if !t_stop.is_finite() {
    return make_error!("Polytomy parent time must be finite, got {t_stop}");
  }
  if !mutation_rate.is_finite() || mutation_rate < 0.0 {
    return make_error!("Polytomy mutation rate must be finite and non-negative, got {mutation_rate}");
  }
  if let Some((index, child)) = children.iter().enumerate().find(|(_, child)| !child.time.is_finite()) {
    return make_error!("Polytomy child {index} time must be finite, got {}", child.time);
  }
  Ok(())
}

/// Move every lineage whose arrival is at or before `elapsed` into the live set.
///
/// `to_come` is ordered most-recent-first, so this is a prefix pop.
fn admit_arrivals(alive: &mut Vec<Tracked>, to_come: &mut VecDeque<Pending>, elapsed: f64) {
  while to_come.front().is_some_and(|next| next.elapsed <= elapsed) {
    let arrived = to_come.pop_front().expect("front was just inspected");
    alive.push(arrived.lineage);
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
  // Integer draws preserve exact branch weights above f64's integer precision limit.
  let mut remaining = rng.gen_range(0..total_mutations);
  for lineage in alive.iter_mut() {
    let mutations = u64::from(lineage.mutations);
    if remaining < mutations {
      lineage.mutations -= 1;
      return Ok(());
    }
    remaining -= mutations;
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
