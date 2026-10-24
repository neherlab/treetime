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
const MIN_PROBE_SHRINK: f64 = 8.0;
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
  let h = inputs.window.half_width;
  let (target, window) = log_rate_target(probe((-h).exp())?, log_lh, probe(h.exp())?, inputs.window);
  Ok(MarginalRateStep {
    clock_rate: reference * (STEP_FRACTION * target).exp(),
    target_rate: reference * target.exp(),
    log_lh,
    window,
  })
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct MarginalRateStep {
  pub clock_rate: f64,
  pub target_rate: f64,
  pub log_lh: f64,
  pub window: RateWindow,
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

fn log_rate_target(below: f64, at: f64, above: f64, window: RateWindow) -> (f64, RateWindow) {
  let h = window.half_width;
  let bound = MAX_TARGET_FACTOR.ln();
  let curvature = above + below - 2.0 * at;
  let maximum = (curvature.is_finite() && curvature < 0.0).then(|| h * (below - above) / (2.0 * curvature));
  let held = RateWindow {
    follow_slope: false,
    ..window
  };
  match maximum {
    Some(target) if target.abs() <= h => (target, held.shrunk()),
    _ if !window.follow_slope => (0.0, held),
    Some(target) => (target.clamp(-bound, bound), held),
    None if above.max(below) <= at => (0.0, held),
    None if above > below => (bound, held),
    None => (-bound, held),
  }
}

/// Probe window of the rate step in log rate. The first step follows the slope of log Z out of the window, up to the
/// 2-fold bound, to correct a poor start rate; later steps move only to a maximum between the probes, and otherwise
/// keep the rate. The half-width halves after each step to a maximum between the probes, down to 1/8 of the start:
/// a parabola in log rate through probes at distance h stops h^2/6 below the maximum of a Poisson-shaped log Z.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct RateWindow {
  pub half_width: f64,
  pub follow_slope: bool,
}

impl RateWindow {
  pub(crate) fn initial() -> Self {
    Self {
      half_width: PROBE_FACTOR.ln(),
      follow_slope: true,
    }
  }

  fn shrunk(self) -> Self {
    Self {
      half_width: (self.half_width / 2.0).max(PROBE_FACTOR.ln() / MIN_PROBE_SHRINK),
      ..self
    }
  }
}
