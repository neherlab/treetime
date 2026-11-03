use crate::branch_lengths::branch_length_or_zero;
use crate::clock::clock_model::ClockLine;
use crate::clock::clock_state::ClockInputs;
use crate::clock::divergence::root_to_node_divergences;
use crate::error::input_error;
use crate::progress::LogSink;
use crate::progress_info;
use eyre::{Report, WrapErr};
use ndarray::Array1;
use ndarray_stats::Quantile1dExt;
use ndarray_stats::interpolate::Linear;
use noisy_float::types::n64;
use ordered_float::OrderedFloat;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::BTreeSet;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::node::GraphNodeKey;

pub(crate) fn clock_filter(
  graph: &Graph,
  inputs: &ClockInputs,
  clock_line: &(impl ClockLine + Sync),
  branch_lengths: &EdgeMap<Option<f64>>,
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
  let mut leaf_clock_deviations: Array1<OrderedFloat<f64>> = deviations_by_leaf
    .iter()
    .map(|(_, deviation)| OrderedFloat(*deviation))
    .collect();

  let n = leaf_clock_deviations.len();
  if n == 0 {
    return Err(input_error("Clock filtering requires at least one dated leaf"));
  }
  let q25 = quartile(&mut leaf_clock_deviations, 0.25)?;
  let q75 = quartile(&mut leaf_clock_deviations, 0.75)?;
  let iqd = q75 - q25;

  let root_children = root_children(graph);
  let outliers: BTreeSet<GraphNodeKey> = deviations_by_leaf
    .iter()
    .filter(|(key, deviation)| deviation.abs() > iqd * threshold && !root_children.contains(key))
    .map(|(key, _)| *key)
    .collect();

  progress_info!(
    log,
    "Outlier filtering: {} leaves flagged as outliers, IQD={iqd:.6e}",
    outliers.len()
  );
  log::debug!(
    "Leaf clock deviations (n={n}): min={:.6e}, Q1={q25:.6e}, Q3={q75:.6e}, max={:.6e}",
    leaf_clock_deviations.iter().min().map_or(0.0, |min| min.into_inner()),
    leaf_clock_deviations.iter().max().map_or(0.0, |max| max.into_inner()),
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
  pub(crate) divergences: NodeMap<f64>,
  pub(crate) iqd: f64,
}

fn quartile(deviations: &mut Array1<OrderedFloat<f64>>, q: f64) -> Result<f64, Report> {
  let quartile = deviations
    .quantile_mut(n64(q), &Linear)
    .wrap_err_with(|| format!("When computing the {q} quantile of the leaf clock deviations"))?;
  Ok(quartile.into_inner())
}

fn root_children(graph: &Graph) -> BTreeSet<GraphNodeKey> {
  graph
    .get_roots()
    .flat_map(|root| graph.children_keys_of(root).map(|(child, _)| child))
    .collect()
}
