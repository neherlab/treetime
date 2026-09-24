use crate::clock::clock_model::ClockLine;
use crate::clock::clock_state::{ClockInputs, ClockState};
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

#[allow(
  clippy::expect_used,
  clippy::integer_division,
  reason = "expect on a value an upstream invariant guarantees is present; integer division is the intended floor division"
)]
#[expect(
  clippy::integer_division_remainder_used,
  reason = "the remainder distributes work evenly across chunks"
)]
pub(crate) fn clock_filter_inplace(
  graph: &Graph,
  inputs: &ClockInputs,
  state: &mut ClockState,
  clock_line: &(impl ClockLine + Sync),
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  threshold: f64,
) -> Result<ClockFilterResult, Report> {
  log::info!("### Filtering outliers (threshold={threshold})");
  log::debug!(
    "Clock model for filtering: rate={:.6e}, intercept={:.4}",
    clock_line.clock_rate(),
    clock_line.intercept()
  );

  state.map_forward(graph, |context| {
    let mut node = context.input.clone();
    let parent_message = if let Some((edge_key, edge)) = context.parent_edge {
      let parent = context.parent.expect("Non-root node must have a parent");
      node.div = parent.div + edge_branch_length(edge_key, branch_lengths);
      Some(edge.clone())
    } else {
      node.div = 0.0;
      None
    };
    Ok(GraphPassNodeOutput { node, parent_message })
  })?;

  let leaf_clock_deviations: Vec<f64> = graph
    .get_leaves()
    .collect::<Vec<_>>()
    .into_par_iter()
    .filter_map(|leaf| {
      let key = leaf.key();
      let div = state.node(key).div;
      let time = inputs.likely_time(key);
      time.map(|time| clock_line.clock_deviation(time, div))
    })
    .collect::<Vec<_>>()
    .into_iter()
    .map(OrderedFloat)
    .sorted()
    .map(OrderedFloat::into_inner)
    .collect();

  let n = leaf_clock_deviations.len();
  if n == 0 {
    return make_error!("Clock filtering requires at least one dated leaf");
  }
  let iq75 = (3 * n) / 4;
  let iq25 = n / 4;
  let iqd = leaf_clock_deviations[iq75] - leaf_clock_deviations[iq25];

  let outlier_updates: Vec<(GraphNodeKey, bool, i32)> = graph
    .get_leaves()
    .collect::<Vec<_>>()
    .into_par_iter()
    .filter_map(|leaf| {
      let key = leaf.key();
      let node = state.node(key);
      let div = node.div;
      let was_outlier = node.is_outlier;
      inputs.likely_time(key).map(|time| {
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

#[derive(Debug, Clone, Copy)]
pub struct ClockFilterResult {
  pub(crate) new_outliers: i32,
  pub(crate) iqd: f64,
}

fn edge_branch_length(edge_key: GraphEdgeKey, branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> f64 {
  branch_lengths[&edge_key].unwrap_or_default()
}
