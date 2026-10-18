//! Clock rate step toward the maximum of the marginal likelihood of the node times.
//!
//! The branch-time distributions are a rate-independent likelihood on a branch-length grid, divided by the clock rate
//! (`fn compute_branch_length_distribution()`). A different rate therefore rescales their time axis, and each
//! candidate rate costs one backward pass, not a new branch-likelihood computation.

use crate::clock::date_constraints::DateConstraints;
use crate::coalescent::mean_field::CoalescentPrior;
use crate::make_report;
use crate::timetree::inference::backward_pass::propagate_distributions_backward;
use crate::timetree::inference::result::BranchLikelihood;
use eyre::Report;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_grid::MaxGridPoints;

/// The marginal likelihood is probed at the reference rate divided and multiplied by this factor
const PROBE_FACTOR: f64 = 1.25;
/// The parabola's maximum is taken at most this factor away from the reference rate
const MAX_TARGET_FACTOR: f64 = 2.0;
/// Fraction of the way, in log rate, that the rate moves toward the target each round
pub(crate) const STEP_FRACTION: f64 = 2.0 / 3.0;

pub(crate) struct MarginalRateInputs<'a> {
  pub graph: &'a Graph,
  pub date_constraints: &'a DateConstraints,
  pub coalescent: Option<&'a CoalescentPrior>,
  pub bad_branches: &'a BTreeMap<GraphNodeKey, bool>,
  /// Branch-time distributions computed at `reference_rate`
  pub branches: &'a BTreeMap<GraphEdgeKey, BranchLikelihood>,
  pub reference_rate: f64,
  pub max_points: MaxGridPoints,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct MarginalRateStep {
  /// Rate `STEP_FRACTION` of the way, in log rate, from the reference rate to `target_rate`
  pub clock_rate: f64,
  /// Maximum of the parabola through the probes, within `MAX_TARGET_FACTOR` of the reference rate
  pub target_rate: f64,
  /// Log marginal likelihood at the reference rate
  pub log_lh: f64,
}

/// Branch-time distributions at `clock_rate`, from those at `reference_rate`: times scale by their ratio
pub(crate) fn branches_at_rate(
  branches: &BTreeMap<GraphEdgeKey, BranchLikelihood>,
  reference_rate: f64,
  clock_rate: f64,
) -> Result<BTreeMap<GraphEdgeKey, BranchLikelihood>, Report> {
  let factor = reference_rate / clock_rate;
  branches
    .iter()
    .map(|(&key, branch)| {
      let distribution = branch
        .distribution
        .as_deref()
        .map(|distribution| distribution.scale_time(factor).map(Arc::new))
        .transpose()?;
      let time_length = branch.time_length.map(|time| time * factor);
      Ok((
        key,
        BranchLikelihood {
          distribution,
          time_length,
        },
      ))
    })
    .collect()
}

/// Log marginal likelihood of the node times at `clock_rate`; `None` when the root has no distribution
pub(crate) fn marginal_log_lh_at_rate(inputs: &MarginalRateInputs<'_>, clock_rate: f64) -> Result<Option<f64>, Report> {
  let branches = branches_at_rate(inputs.branches, inputs.reference_rate, clock_rate)?;
  let backward = propagate_distributions_backward(
    inputs.graph,
    inputs.date_constraints,
    inputs.coalescent,
    inputs.bad_branches,
    &branches,
    inputs.max_points,
  )?;
  backward.log_marginal_likelihood(inputs.graph.root_key()?)
}

/// One damped step of the clock rate toward the maximum of the log marginal likelihood. A parabola in log rate through
/// the reference rate and the rates `PROBE_FACTOR` below and above gives a target, at most `MAX_TARGET_FACTOR` away;
/// without a maximum the target is that bound on the side of the better probe. The log rate x moves the fraction a =
/// `STEP_FRACTION` of the way to the target T: x' = x + a (T - x). The refinement loop repeats the step every round,
/// with the times, Tc and the reconstruction re-estimated in between, so the rate needs no precise optimum. Near the
/// fixed point x* the target follows the rate, T = x* + s (x - x*), and the distance to x* shrinks by 1 - a (1 - s) per
/// round: a full step (a = 1) oscillates where the target overshoots (s near -1), and a < 1 damps that at the cost of
/// slower convergence where the target lags behind the rate (0 < s < 1). A probe whose backward pass fails or leaves
/// the root without a distribution counts as impossible.
pub(crate) fn step_marginal_clock_rate(inputs: &MarginalRateInputs<'_>) -> Result<MarginalRateStep, Report> {
  let reference = inputs.reference_rate;
  if !(reference.is_finite() && reference > 0.0) {
    return Err(make_report!(
      "Marginal clock rate step needs a positive reference rate, got {reference}"
    ));
  }
  let probe = |factor: f64| {
    marginal_log_lh_at_rate(inputs, reference * factor)
      .ok()
      .flatten()
      .unwrap_or(f64::NEG_INFINITY)
  };
  let log_lh = probe(1.0);
  if !log_lh.is_finite() {
    return Err(make_report!(
      "Marginal clock rate step found no finite marginal likelihood at the reference rate {reference}"
    ));
  }
  let target = log_rate_target(probe(1.0 / PROBE_FACTOR), log_lh, probe(PROBE_FACTOR));
  Ok(MarginalRateStep {
    clock_rate: reference * (STEP_FRACTION * target).exp(),
    target_rate: reference * target.exp(),
    log_lh,
  })
}

/// Log of the target rate relative to the reference: the maximum of the parabola through `(-h, below)`, `(0, at)` and
/// `(h, above)` with h = ln(PROBE_FACTOR), bounded by ln(MAX_TARGET_FACTOR)
fn log_rate_target(below: f64, at: f64, above: f64) -> f64 {
  let h = PROBE_FACTOR.ln();
  let bound = MAX_TARGET_FACTOR.ln();
  let curvature = above + below - 2.0 * at;
  let target = if curvature.is_finite() && curvature < 0.0 {
    h * (below - above) / (2.0 * curvature)
  } else if above.max(below) <= at {
    0.0
  } else if above > below {
    bound
  } else {
    -bound
  };
  target.clamp(-bound, bound)
}
