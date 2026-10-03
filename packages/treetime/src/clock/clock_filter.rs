use crate::branch_lengths::branch_length_or_zero;
use crate::clock::clock_model::ClockLine;
use crate::clock::clock_state::ClockInputs;
use crate::clock::divergence::root_to_node_divergences;
use crate::error::input_error;
use crate::progress::LogSink;
use crate::progress_info;
use eyre::Report;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

#[expect(
  clippy::integer_division,
  clippy::integer_division_remainder_used,
  reason = "the quartile indices are floor divisions of the leaf count"
)]
pub(crate) fn clock_filter(
  graph: &Graph,
  inputs: &ClockInputs,
  clock_line: &(impl ClockLine + Sync),
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  threshold: f64,
  log: &dyn LogSink,
) -> Result<ClockFilterResult, Report> {
  progress_info!(log, "### Filtering outliers (threshold={threshold})");
  log::debug!(
    "Clock model for filtering: rate={:.6e}, intercept={:.4}",
    clock_line.clock_rate(),
    clock_line.intercept()
  );

  let divergences = root_to_node_divergences(graph, |edge_key| branch_length_or_zero(branch_lengths, edge_key))?;

  let deviations_by_leaf: Vec<(GraphNodeKey, f64)> = graph
    .get_leaves()
    .collect::<Vec<_>>()
    .into_par_iter()
    .filter_map(|leaf| {
      let key = leaf.key();
      let time = inputs.likely_time(key)?;
      Some((key, clock_line.clock_deviation(time, divergences[&key])))
    })
    .collect();
  let leaf_clock_deviations: Vec<f64> = deviations_by_leaf
    .iter()
    .map(|(_, deviation)| OrderedFloat(*deviation))
    .sorted()
    .map(OrderedFloat::into_inner)
    .collect();

  let n = leaf_clock_deviations.len();
  if n == 0 {
    return Err(input_error("Clock filtering requires at least one dated leaf"));
  }
  let iq75 = (3 * n) / 4;
  let iq25 = n / 4;
  let iqd = leaf_clock_deviations[iq75] - leaf_clock_deviations[iq25];

  let outliers: BTreeSet<GraphNodeKey> = deviations_by_leaf
    .iter()
    .filter(|(_, deviation)| deviation.abs() > iqd * threshold)
    .map(|(key, _)| *key)
    .collect();

  progress_info!(
    log,
    "Outlier filtering: {} leaves flagged as outliers, IQD={iqd:.6e}",
    outliers.len()
  );
  log::debug!(
    "Leaf clock deviations (n={}): min={:.6e}, Q1={:.6e}, Q3={:.6e}, max={:.6e}",
    n,
    leaf_clock_deviations.first().copied().unwrap_or(0.0),
    leaf_clock_deviations[iq25],
    leaf_clock_deviations[iq75],
    leaf_clock_deviations.last().copied().unwrap_or(0.0)
  );

  Ok(ClockFilterResult {
    outliers,
    divergences,
    iqd,
  })
}

#[derive(Debug)]
pub(crate) struct ClockFilterResult {
  pub(crate) outliers: BTreeSet<GraphNodeKey>,
  pub(crate) divergences: BTreeMap<GraphNodeKey, f64>,
  pub(crate) iqd: f64,
}
