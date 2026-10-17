//! Clock rate that maximizes the marginal likelihood of the node times.
//!
//! The branch-time distributions are a rate-independent likelihood on a branch-length grid, divided by the clock rate
//! (`fn compute_branch_length_distribution()`). A different rate therefore rescales their time axis, and each
//! candidate rate costs one backward pass, not a new branch-likelihood computation.

use crate::brent::minimize_brent;
use crate::clock::date_constraints::DateConstraints;
use crate::coalescent::coalescent::CoalescentModel;
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

/// The search covers the reference rate divided and multiplied by this factor; the refinement loop repeats the search
/// every round, so a rate beyond it is reached in a later round
const RATE_SEARCH_FACTOR: f64 = 3.0;
/// Relative tolerance on the square root of the rate
const RATE_SEARCH_TOLERANCE: f64 = 1e-3;
const RATE_SEARCH_MAX_ITERS: u64 = 20;

pub(crate) struct MarginalRateInputs<'a> {
  pub graph: &'a Graph,
  pub date_constraints: &'a DateConstraints,
  pub coalescent: Option<&'a CoalescentModel>,
  pub bad_branches: &'a BTreeMap<GraphNodeKey, bool>,
  /// Branch-time distributions computed at `reference_rate`
  pub branches: &'a BTreeMap<GraphEdgeKey, BranchLikelihood>,
  pub reference_rate: f64,
  pub max_points: MaxGridPoints,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct MarginalRateFit {
  pub clock_rate: f64,
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

/// Maximizes the log marginal likelihood over the clock rate by Brent's method on its square root, which keeps the rate
/// positive. A rate whose backward pass fails or leaves the root without a distribution scores as impossible.
pub(crate) fn optimize_marginal_clock_rate(inputs: &MarginalRateInputs<'_>) -> Result<MarginalRateFit, Report> {
  let reference = inputs.reference_rate;
  if !(reference.is_finite() && reference > 0.0) {
    return Err(make_report!(
      "Marginal clock rate search needs a positive reference rate, got {reference}"
    ));
  }
  let (lower, upper) = (
    (reference / RATE_SEARCH_FACTOR).sqrt(),
    (reference * RATE_SEARCH_FACTOR).sqrt(),
  );
  let minimum = minimize_brent(
    lower,
    upper,
    Some((RATE_SEARCH_TOLERANCE, f64::EPSILON)),
    RATE_SEARCH_MAX_ITERS,
    |s| {
      let log_lh = marginal_log_lh_at_rate(inputs, s * s).ok().flatten();
      Ok(log_lh.map_or(f64::INFINITY, |log_lh| -log_lh))
    },
  )?;
  if !minimum.value.is_finite() {
    return Err(make_report!(
      "Marginal clock rate search found no rate between {} and {} with a finite marginal likelihood",
      lower * lower,
      upper * upper
    ));
  }
  Ok(MarginalRateFit {
    clock_rate: minimum.x * minimum.x,
    log_lh: -minimum.value,
  })
}
