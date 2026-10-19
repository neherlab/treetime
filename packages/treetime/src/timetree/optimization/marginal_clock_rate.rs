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

const PROBE_FACTOR: f64 = 1.25;
const MAX_TARGET_FACTOR: f64 = 2.0;
pub(crate) const STEP_FRACTION: f64 = 2.0 / 3.0;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum RateWindow {
  FollowSlope,
  InteriorMaximum,
}

pub(crate) struct MarginalRateInputs<'a> {
  pub graph: &'a Graph,
  pub date_constraints: &'a DateConstraints,
  pub coalescent: Option<&'a CoalescentPrior>,
  pub bad_branches: &'a BTreeMap<GraphNodeKey, bool>,
  pub branches: &'a BTreeMap<GraphEdgeKey, BranchLikelihood>,
  pub reference_rate: f64,
  pub window: RateWindow,
  pub max_points: MaxGridPoints,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct MarginalRateStep {
  pub clock_rate: f64,
  pub target_rate: f64,
  pub log_lh: f64,
}

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
  let target = log_rate_target(probe(1.0 / PROBE_FACTOR), log_lh, probe(PROBE_FACTOR), inputs.window);
  Ok(MarginalRateStep {
    clock_rate: reference * (STEP_FRACTION * target).exp(),
    target_rate: reference * target.exp(),
    log_lh,
  })
}

fn log_rate_target(below: f64, at: f64, above: f64, window: RateWindow) -> f64 {
  let h = PROBE_FACTOR.ln();
  let bound = MAX_TARGET_FACTOR.ln();
  let curvature = above + below - 2.0 * at;
  let maximum = (curvature.is_finite() && curvature < 0.0).then(|| h * (below - above) / (2.0 * curvature));
  match (window, maximum) {
    (RateWindow::InteriorMaximum, Some(target)) if target.abs() <= h => target,
    (RateWindow::InteriorMaximum, _) => 0.0,
    (RateWindow::FollowSlope, Some(target)) => target.clamp(-bound, bound),
    (RateWindow::FollowSlope, None) if above.max(below) <= at => 0.0,
    (RateWindow::FollowSlope, None) if above > below => bound,
    (RateWindow::FollowSlope, None) => -bound,
  }
}
