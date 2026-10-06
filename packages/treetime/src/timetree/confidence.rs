use crate::clock::clock_model::{ClockModel, ClockModelStats};
use crate::coalescent::coalescent::CoalescentModel;
use crate::make_error;
use crate::progress::LogSink;
use crate::timetree::inference::result::{NodePosterior, TimeInference};
use crate::timetree::inference::runner::{TimeInferenceInputs, run_timetree};
use crate::{progress_info, progress_warn};
use deser::Serialize;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use ordered_float::OrderedFloat;
use statrs::function::erf::erf_inv;
use std::collections::BTreeMap;
use std::f64::consts::SQRT_2;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

pub const CI_FRACTION: f64 = 0.9;

const CI_LOWER_QUANTILE: f64 = (1.0 - CI_FRACTION) * 0.5;
const CI_UPPER_QUANTILE: f64 = 1.0 - (1.0 - CI_FRACTION) * 0.5;

pub(crate) fn compute_rate_susceptibility(
  time_inputs: &TimeInferenceInputs<'_>,
  coalescent: Option<&CoalescentModel>,
  rate_std: f64,
  log: &dyn LogSink,
) -> Result<RateSusceptibility, Report> {
  let graph = time_inputs.graph;
  let current_rate = time_inputs.clock_model.clock_rate();

  let upper_rate = current_rate + rate_std;
  let lower_rate = (0.1 * current_rate).max(current_rate - rate_std);

  let run_scaled = |scale: f64| {
    let scaled_gammas = time_inputs
      .gammas
      .iter()
      .map(|(key, gamma)| (*key, gamma * scale))
      .collect::<BTreeMap<_, _>>();
    let scaled_inputs = TimeInferenceInputs {
      gammas: &scaled_gammas,
      ..*time_inputs
    };
    run_timetree(&scaled_inputs, coalescent, log)
  };

  progress_info!(log, "Rate susceptibility: running with upper rate {upper_rate:.6e}");
  let upper = run_scaled(upper_rate / current_rate).wrap_err("Rate susceptibility: timetree at upper rate failed")?;

  progress_info!(log, "Rate susceptibility: running with lower rate {lower_rate:.6e}");
  let lower = run_scaled(lower_rate / current_rate).wrap_err("Rate susceptibility: timetree at lower rate failed")?;

  progress_info!(log, "Rate susceptibility: running with central rate {current_rate:.6e}");
  let central = run_scaled(1.0).wrap_err("Rate susceptibility: timetree at central rate failed")?;

  let dates = graph
    .get_nodes()
    .filter_map(|node_ref| {
      let key = node_ref.key();
      let central_date = central.posterior[&key].time?;
      let upper_date = upper.posterior[&key].time?;
      let lower_date = lower.posterior[&key].time?;
      let mut dates = [lower_date, central_date, upper_date];
      dates.sort_by_key(|date| OrderedFloat(*date));
      Some((key, dates))
    })
    .collect();

  progress_info!(log, "Rate susceptibility analysis completed");
  Ok(RateSusceptibility { dates, central })
}

pub(crate) struct RateSusceptibility {
  pub dates: BTreeMap<GraphNodeKey, [f64; 3]>,
  pub central: TimeInference,
}

pub fn extract_confidence_intervals(
  graph: &Graph,
  posterior: &BTreeMap<GraphNodeKey, NodePosterior>,
  rate_susceptibility_dates: &BTreeMap<GraphNodeKey, [f64; 3]>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Vec<NodeConfidenceInterval> {
  graph
    .get_nodes()
    .filter_map(|node| {
      let key = node.key();
      let node_posterior = &posterior[&key];
      let name = names[&key].clone().unwrap_or_default();
      let date = node_posterior.time?;

      let mutation_contribution: Option<(f64, f64)> = None;

      let rate_contribution = rate_susceptibility_dates
        .get(&key)
        .map(|dates| date_uncertainty_due_to_rate(*dates, (CI_LOWER_QUANTILE, CI_UPPER_QUANTILE)));

      let (lower, upper) = if rate_contribution.is_none() && mutation_contribution.is_none() {
        (date, date)
      } else {
        let limits = node_posterior
          .distribution
          .as_ref()
          .and_then(|dist| dist.time_bounds())
          .unwrap_or((f64::NEG_INFINITY, f64::INFINITY));
        combine_confidence(date, limits, rate_contribution, mutation_contribution)
      };

      let lower = lower.min(date);
      let upper = upper.max(date);

      Some(NodeConfidenceInterval {
        key,
        name,
        date,
        lower,
        upper,
      })
    })
    .sorted_by_key(|ci| ci.key)
    .collect_vec()
}

pub(crate) fn date_uncertainty_due_to_rate(dates: [f64; 3], interval: (f64, f64)) -> (f64, f64) {
  let [lower, central, upper] = dates;
  let z_lower = quantile_to_zscore(interval.0);
  let z_upper = quantile_to_zscore(interval.1);
  let ci_lower = central + z_lower * (lower - central).abs();
  let ci_upper = central + z_upper * (upper - central).abs();
  (ci_lower, ci_upper)
}

#[derive(Debug, Clone, Serialize)]
pub struct NodeConfidenceInterval {
  #[deser(skip)]
  pub key: GraphNodeKey,
  pub name: String,
  pub date: f64,
  pub lower: f64,
  pub upper: f64,
}

pub(crate) fn combine_confidence(
  center: f64,
  limits: (f64, f64),
  c1: Option<(f64, f64)>,
  c2: Option<(f64, f64)>,
) -> (f64, f64) {
  let (min_val, max_val) = match (c1, c2) {
    (None, None) => return limits,
    (Some(c), None) | (None, Some(c)) => c,
    (Some(c1), Some(c2)) => {
      let min_val = center - (c1.0 - center).hypot(c2.0 - center);
      let max_val = center + (c1.1 - center).hypot(c2.1 - center);
      (min_val, max_val)
    },
  };

  (limits.0.max(min_val), limits.1.min(max_val))
}

pub(crate) fn determine_rate_std(
  clock_std_dev: Option<f64>,
  covariation: bool,
  clock_model: &ClockModel,
  log: &dyn LogSink,
) -> Result<Option<f64>, Report> {
  if let Some(std_dev) = clock_std_dev {
    if std_dev <= 0.0 {
      return make_error!("--clock-std-dev must be positive, got {std_dev}");
    }
    return Ok(Some(std_dev));
  }

  if !covariation {
    return Ok(None);
  }

  match clock_model.stats() {
    ClockModelStats::Estimated(stats) => {
      let rate_variance = stats.cov[[0, 0]];
      if rate_variance <= 0.0 {
        progress_warn!(
          log,
          "Rate variance from regression covariance is non-positive ({rate_variance:.4e}), skipping rate susceptibility"
        );
        return Ok(None);
      }
      Ok(Some(rate_variance.sqrt()))
    },
    ClockModelStats::Fixed => Ok(None),
  }
}

pub(crate) fn quantile_to_zscore(p: f64) -> f64 {
  if p * (1.0 - p) == 0.0 {
    return 0.0;
  }
  SQRT_2 * erf_inv(2.0 * p - 1.0)
}
