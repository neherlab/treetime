use crate::clock::clock_model::ClockLine;
use crate::clock::clock_state::ClockState;
use crate::make_error;
use eyre::Report;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use rayon::prelude::*;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pass::GraphPassNodeOutput;

#[derive(Debug, Clone, Copy)]
pub struct ClockFilterResult {
  pub new_outliers: i32,
  pub iqd: f64,
}

/// Filter outliers based on clock model residuals.
///
/// Marks leaves as outliers if their clock deviation exceeds `threshold * IQD`
/// where IQD is the interquartile distance of clock deviations.
///
/// Accepts any `ClockLine` implementor: both validated `ClockModel` (positive
/// rate) and raw `ClockRegression` (any rate sign, used in pre-filter path).
#[allow(clippy::integer_division_remainder_used)]
pub fn clock_filter_inplace<D>(
  graph: &Graph<D>,
  state: &mut ClockState,
  clock_line: &(impl ClockLine + Sync),
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  threshold: f64,
) -> Result<ClockFilterResult, Report>
where
  D: Send + Sync,
{
  log::info!("### Filtering outliers (threshold={threshold})");
  log::debug!(
    "Clock model for filtering: rate={:.6e}, intercept={:.4}",
    clock_line.clock_rate(),
    clock_line.intercept()
  );

  // Assign divergence to each node: div = parent.div + branch_length, parents before children.
  state.map_forward(graph, |context| {
    let mut node = context.input;
    let parent_message = if let Some((edge_key, edge)) = context.parent_edge {
      let parent = context.parent.expect("Non-root node must have a parent");
      node.div = parent.div + edge_branch_length(edge_key, branch_lengths);
      Some(edge)
    } else {
      node.div = 0.0;
      None
    };
    Ok(GraphPassNodeOutput { node, parent_message })
  })?;

  // collect clock_deviation of leaf nodes into a vector
  let leaf_clock_deviations: Vec<f64> = graph
    .get_leaves()
    .par_iter()
    .filter_map(|leaf| {
      let node = state.node(leaf.read_arc().key());
      let div = node.div;
      let time = node.likely_time();
      time.map(|time| clock_line.clock_deviation(time, div))
    })
    .collect::<Vec<_>>()
    .into_iter()
    .map(OrderedFloat)
    .sorted()
    .map(OrderedFloat::into_inner)
    .collect();

  // calculate the interquartile range by taking the difference between the 3/4 and 1/4 quantile
  let n = leaf_clock_deviations.len();
  if n == 0 {
    return make_error!("Clock filtering requires at least one dated leaf");
  }
  let iq75 = (3 * n) / 4;
  let iq25 = n / 4;
  let iqd = leaf_clock_deviations[iq75] - leaf_clock_deviations[iq25];

  // Compute each leaf's outlier decision in parallel from the clock state (read-only), then apply
  // the results serially. A `BTreeMap` cannot be mutated concurrently, so the write phase is split
  // out; the decisions are per-leaf independent and the count is order-free, so this is bit-identical
  // to the payload-based per-node locked write.
  let outlier_updates: Vec<(GraphNodeKey, bool, i32)> = graph
    .get_leaves()
    .par_iter()
    .filter_map(|leaf| {
      let key = leaf.read_arc().key();
      let node = state.node(key);
      let div = node.div;
      let was_outlier = node.is_outlier;
      node.likely_time().map(|time| {
        let clock_deviation = clock_line.clock_deviation(time, div);
        let is_outlier = clock_deviation.abs() > iqd * threshold;
        (key, is_outlier, i32::from(was_outlier != is_outlier))
      })
    })
    .collect();

  let new_outliers = outlier_updates.iter().map(|(_, _, changed)| *changed).sum();
  for (key, is_outlier, _) in outlier_updates {
    state.node_mut(key).is_outlier = is_outlier;
  }

  log::info!("Outlier filtering: {new_outliers} leaves changed status, IQD={iqd:.6e}");
  log::debug!(
    "Leaf clock deviations (n={}): min={:.6e}, Q1={:.6e}, Q3={:.6e}, max={:.6e}",
    n,
    leaf_clock_deviations.first().copied().unwrap_or(0.0),
    leaf_clock_deviations[iq25],
    leaf_clock_deviations[iq75],
    leaf_clock_deviations.last().copied().unwrap_or(0.0)
  );

  Ok(ClockFilterResult { new_outliers, iqd })
}

/// Branch length of an edge (an input), read from the value map, defaulting to `0.0` when unset.
fn edge_branch_length(edge_key: GraphEdgeKey, branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> f64 {
  branch_lengths[&edge_key].unwrap_or_default()
}
