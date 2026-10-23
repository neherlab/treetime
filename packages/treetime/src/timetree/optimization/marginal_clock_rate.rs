use crate::clock::date_constraints::DateConstraints;
use crate::coalescent::mean_field::CoalescentPrior;
use crate::make_report;
use crate::timetree::inference::backward_pass::propagate_distributions_backward;
use crate::timetree::inference::result::{BranchLikelihood, log_marginal_likelihood};
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_grid::MaxGridPoints;

const PROBE_FACTOR: f64 = 1.25;
const MAX_TARGET_FACTOR: f64 = 2.0;
pub(crate) const STEP_FRACTION: f64 = 2.0 / 3.0;

pub(crate) fn step_marginal_clock_rate(inputs: &MarginalRateInputs<'_>) -> Result<MarginalRateStep, Report> {
  let reference = inputs.reference_rate;
  if !(reference.is_finite() && reference > 0.0) {
    return Err(make_report!(
      "Marginal clock rate step needs a positive reference rate, got {reference}"
    ));
  }
  let probe = |factor: f64| -> Result<f64, Report> {
    let rate = reference * factor;
    let log_lh = marginal_log_lh_at_rate(inputs, rate)
      .wrap_err_with(|| format!("When computing the marginal likelihood of the node times at clock rate {rate}"))?;
    Ok(log_lh.unwrap_or(f64::NEG_INFINITY))
  };
  let log_lh = probe(1.0)?;
  if !log_lh.is_finite() {
    return Err(make_report!(
      "Marginal clock rate step found no finite marginal likelihood at the reference rate {reference}"
    ));
  }
  let research = Research::from_env();
  let window = if research.guard {
    inputs.window
  } else {
    RateWindow::FollowSlope
  };
  let target = log_rate_target(
    probe(1.0 / PROBE_FACTOR)?,
    log_lh,
    probe(PROBE_FACTOR)?,
    window,
    research.power,
    research.bound,
  );
  Ok(MarginalRateStep {
    clock_rate: reference * (research.step * target).exp(),
    target_rate: reference * target.exp(),
    log_lh,
  })
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct MarginalRateStep {
  pub clock_rate: f64,
  pub target_rate: f64,
  pub log_lh: f64,
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
  log_marginal_likelihood(&backward, inputs.graph.root_key()?, inputs.coalescent)
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

// RESEARCH (not for commit): coordinate, guards and step fraction of the rate step from the environment
struct Research {
  power: f64,
  guard: bool,
  step: f64,
  bound: f64,
}

impl Research {
  fn from_env() -> Self {
    let power = match std::env::var("TREETIME_RATE_COORD").as_deref() {
      Ok("sqrt") => 0.5,
      Ok("cbrt") => 1.0 / 3.0,
      Ok("rate") => 1.0,
      _ => 0.0,
    };
    let guard = std::env::var("TREETIME_RATE_GUARD").as_deref() != Ok("none");
    let step = std::env::var("TREETIME_RATE_STEP")
      .ok()
      .and_then(|s| s.parse().ok())
      .unwrap_or(STEP_FRACTION);
    let bound = std::env::var("TREETIME_RATE_BOUND")
      .ok()
      .and_then(|s| s.parse().ok())
      .unwrap_or(MAX_TARGET_FACTOR);
    Self {
      power,
      guard,
      step,
      bound,
    }
  }
}

/// Target in log rate from the parabola through the probes in u = exp(power x) (u = x for power 0)
#[expect(
  clippy::suspicious_operation_groupings,
  reason = "the coefficients of the parabola through three points pair each squared abscissa with the other two ordinates"
)]
fn log_rate_target(below: f64, at: f64, above: f64, window: RateWindow, power: f64, max_factor: f64) -> f64 {
  let h = PROBE_FACTOR.ln();
  let bound = max_factor.ln();
  let maximum = if power == 0.0 {
    let curvature = above + below - 2.0 * at;
    (curvature.is_finite() && curvature < 0.0).then(|| h * (below - above) / (2.0 * curvature))
  } else {
    let [u1, u2, u3] = [-h, 0.0, h].map(|x| (power * x).exp());
    let denom = (u1 - u2) * (u1 - u3) * (u2 - u3);
    let a = (u3 * (at - below) + u2 * (below - above) + u1 * (above - at)) / denom;
    let b = (u3 * u3 * (below - at) + u2 * u2 * (above - below) + u1 * u1 * (at - above)) / denom;
    let vertex = -b / (2.0 * a);
    match (a.is_finite() && a < 0.0, vertex > 0.0) {
      (true, true) => Some(vertex.ln() / power),
      (true, false) => Some(-bound), // the maximum lies at zero rate: as far down as the bound allows
      (false, _) => None,
    }
  };
  match (window, maximum) {
    (RateWindow::InteriorMaximum, Some(target)) if target.abs() <= h => target,
    (RateWindow::InteriorMaximum, _) => 0.0,
    (RateWindow::FollowSlope, Some(target)) => target.clamp(-bound, bound),
    (RateWindow::FollowSlope, None) if above.max(below) <= at => 0.0,
    (RateWindow::FollowSlope, None) if above > below => bound,
    (RateWindow::FollowSlope, None) => -bound,
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum RateWindow {
  FollowSlope,
  InteriorMaximum,
}
