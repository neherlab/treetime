use crate::clock::clock_model::ClockLine;
use crate::make_error;
use crate::payload::traits::{ClockEdge, ClockNode};
use eyre::Report;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use rayon::prelude::*;
use treetime_graph::breadth_first::GraphTraversalContinuation;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;

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
pub fn clock_filter_inplace<N, E, D>(
  graph: &Graph<N, E, D>,
  clock_line: &(impl ClockLine + Sync),
  threshold: f64,
) -> Result<ClockFilterResult, Report>
where
  N: GraphNode + ClockNode,
  E: GraphEdge + ClockEdge,
  D: Send + Sync,
{
  log::info!("### Filtering outliers (threshold={threshold})");
  log::debug!(
    "Clock model for filtering: rate={:.6e}, intercept={:.4}",
    clock_line.clock_rate(),
    clock_line.intercept()
  );

  // assign divergence to each node
  graph.par_iter_breadth_first_forward(|mut node| {
    let div = node.get_exactly_one_parent().map_or(0.0, |(parent, edge)| {
      let branch_length = edge.read_arc().branch_length().unwrap_or_default();
      let parent_div = parent.read_arc().div();
      parent_div + branch_length
    });
    node.payload.set_div(div);
    Ok(GraphTraversalContinuation::Continue)
  })?;

  // collect clock_deviation of leaf nodes into a vector
  let leaf_clock_deviations: Vec<f64> = graph
    .get_leaves()
    .par_iter()
    .filter_map(|leaf| {
      let node = leaf.read_arc();
      let payload_arc = node.payload();
      let payload = payload_arc.read();
      let div = payload.div();
      let time = payload.likely_time();
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

  // loop over the leaf nodes and mark the outliers if the absolute value of the deviation is greater than the threshold
  let new_outliers = graph
    .get_leaves()
    .par_iter()
    .map(|leaf| {
      let node = leaf.read_arc();
      let payload_arc = node.payload();
      let (div, time, was_outlier) = {
        let payload = payload_arc.read();
        (payload.div(), payload.likely_time(), payload.is_outlier())
      };
      if let Some(time) = time {
        let clock_deviation = clock_line.clock_deviation(time, div);
        let is_outlier = clock_deviation.abs() > iqd * threshold;
        payload_arc.write().set_is_outlier(is_outlier);
        i32::from(was_outlier != is_outlier)
      } else {
        0
      }
    })
    .sum();

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
