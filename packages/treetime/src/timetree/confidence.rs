use crate::clock::clock_model::{ClockModel, ClockModelStats};
use crate::clock::clock_state::ClockState;
use crate::clock::date_constraints::DateConstraints;
use crate::coalescent::coalescent::CoalescentModel;
use crate::make_error;
use crate::partition::timetree::partition::PartitionTimetree;
use crate::timetree::inference::runner::run_timetree;
use crate::timetree::timetree_state::TimetreeState;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use log::{info, warn};
use ordered_float::OrderedFloat;
use serde::Serialize;
use statrs::function::erf::erf_inv;
use std::collections::BTreeMap;
use std::f64::consts::SQRT_2;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

const CI_FRACTION: f64 = 0.9;

const CI_LOWER_QUANTILE: f64 = (1.0 - CI_FRACTION) * 0.5;
const CI_UPPER_QUANTILE: f64 = 1.0 - (1.0 - CI_FRACTION) * 0.5;

#[allow(clippy::too_many_arguments)]
pub fn compute_rate_susceptibility(
  graph: &mut Graph,
  constraints: &DateConstraints,
  partitions: &[PartitionTimetree],
  clock_model: &ClockModel,
  coalescent: Option<&CoalescentModel>,
  rate_std: f64,
  no_indels: bool,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  state: &mut TimetreeState,
  clock_state: &mut ClockState,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<BTreeMap<GraphNodeKey, [f64; 3]>, Report> {
  let current_rate = clock_model.clock_rate();

  let upper_rate = current_rate + rate_std;
  let lower_rate = (0.1 * current_rate).max(current_rate - rate_std);

  let original_gammas = save_gammas(state);

  let run_branch_lengths = branch_lengths;
  let run_names = names;

  scale_gammas(state, &original_gammas, upper_rate / current_rate);
  info!("Rate susceptibility: running with upper rate {upper_rate:.6e}");
  *state = run_timetree(
    graph,
    constraints,
    partitions,
    run_branch_lengths,
    run_names,
    clock_model,
    coalescent,
    no_indels,
    std::mem::take(state),
    clock_state,
  )
  .wrap_err("Rate susceptibility: timetree at upper rate failed")?;
  let upper_dates = collect_node_times(state);

  scale_gammas(state, &original_gammas, lower_rate / current_rate);
  info!("Rate susceptibility: running with lower rate {lower_rate:.6e}");
  *state = run_timetree(
    graph,
    constraints,
    partitions,
    run_branch_lengths,
    run_names,
    clock_model,
    coalescent,
    no_indels,
    std::mem::take(state),
    clock_state,
  )
  .wrap_err("Rate susceptibility: timetree at lower rate failed")?;
  let lower_dates = collect_node_times(state);

  scale_gammas(state, &original_gammas, 1.0);
  info!("Rate susceptibility: running with central rate {current_rate:.6e}");
  *state = run_timetree(
    graph,
    constraints,
    partitions,
    run_branch_lengths,
    run_names,
    clock_model,
    coalescent,
    no_indels,
    std::mem::take(state),
    clock_state,
  )
  .wrap_err("Rate susceptibility: timetree at central rate failed")?;

  let mut rate_susceptibility_dates = BTreeMap::new();
  for node_ref in graph.get_nodes() {
    let key = node_ref.key();

    let central_date = state.node(key).time;
    let upper_date = upper_dates.get(&key).copied();
    let lower_date = lower_dates.get(&key).copied();

    if let (Some(c), Some(u), Some(l)) = (central_date, upper_date, lower_date) {
      let mut dates = [l, c, u];
      dates.sort_by_key(|x| OrderedFloat(*x));
      rate_susceptibility_dates.insert(key, dates);
    }
  }

  info!("Rate susceptibility analysis completed");
  Ok(rate_susceptibility_dates)
}

pub(crate) fn date_uncertainty_due_to_rate(dates: [f64; 3], interval: (f64, f64)) -> (f64, f64) {
  let [lower, central, upper] = dates;
  let z_lower = quantile_to_zscore(interval.0);
  let z_upper = quantile_to_zscore(interval.1);
  let ci_lower = central + z_lower * (lower - central).abs();
  let ci_upper = central + z_upper * (upper - central).abs();
  (ci_lower, ci_upper)
}

pub fn extract_confidence_intervals(
  graph: &Graph,
  state: &TimetreeState,
  rate_susceptibility_dates: &BTreeMap<GraphNodeKey, [f64; 3]>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Vec<NodeConfidenceInterval> {
  graph
    .get_nodes()
    .filter_map(|node_ref| {
      let node = node_ref;
      let key = node.key();
      let node_state = state.node(key);
      let name = names[&key].clone().unwrap_or_default();
      let date = node_state.time?;

      let mutation_contribution: Option<(f64, f64)> = None;

      let rate_contribution = rate_susceptibility_dates
        .get(&key)
        .map(|dates| date_uncertainty_due_to_rate(*dates, (CI_LOWER_QUANTILE, CI_UPPER_QUANTILE)));

      let (lower, upper) = if rate_contribution.is_none() && mutation_contribution.is_none() {
        (date, date)
      } else {
        let limits = node_state
          .time_distribution
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

#[derive(Debug, Clone, Serialize)]
pub struct NodeConfidenceInterval {
  #[serde(skip)]
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
        warn!(
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

fn save_gammas(state: &TimetreeState) -> Vec<(GraphEdgeKey, f64)> {
  state.edges.iter().map(|(key, edge)| (*key, edge.gamma)).collect_vec()
}

fn scale_gammas(state: &mut TimetreeState, original_gammas: &[(GraphEdgeKey, f64)], scale_factor: f64) {
  for &(key, orig_gamma) in original_gammas {
    if let Some(edge) = state.edges.get_mut(&key) {
      edge.gamma = orig_gamma * scale_factor;
    }
  }
}

fn collect_node_times(state: &TimetreeState) -> BTreeMap<GraphNodeKey, f64> {
  state
    .nodes
    .iter()
    .filter_map(|(key, node)| node.time.map(|time| (*key, time)))
    .collect()
}
