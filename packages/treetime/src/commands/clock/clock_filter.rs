use crate::commands::clock::clock_graph::ClockGraph;
use crate::commands::clock::clock_model::ClockModel;
use crate::graph::breadth_first::GraphTraversalContinuation;
use itertools::Itertools;
use ordered_float::OrderedFloat;

/// Get results of the root-to-tip clock inference.
#[allow(clippy::integer_division_remainder_used)]
pub fn clock_filter_inplace(graph: &ClockGraph, clock_model: &ClockModel, clock_filter_threshold: f64) -> i32 {
  // assign divergence to each node
  graph.par_iter_breadth_first_forward(|mut node| {
    node.payload.div = node
      .get_exactly_one_parent()
      .map(|(parent, edge)| {
        let branch_length = edge.read_arc().branch_length.unwrap_or_default();
        let parent_div = parent.read_arc().div;
        parent_div + branch_length
      })
      .unwrap_or(0.0);
    GraphTraversalContinuation::Continue
  });

  // collect clock_deviation of leaf nodes into a vector
  let leaf_clock_deviations: Vec<f64> = graph
    .get_leaves()
    .iter()
    .map(|leaf| {
      let div = leaf.read_arc().payload().read().div;
      let date = leaf.read_arc().payload().read().date.unwrap();
      clock_model.clock_deviation(date, div)
    })
    .map(OrderedFloat)
    .sorted()
    .map(OrderedFloat::into_inner)
    .collect();

  // calculate the interquartile range by taking the difference between the 3/4 and 1/4 quantile
  let n = leaf_clock_deviations.len();
  let iq75 = (3 * n) / 4;
  let iq25 = n / 4;
  let iqd = leaf_clock_deviations[iq75] - leaf_clock_deviations[iq25];

  let mut new_outliers: i32 = 0;
  // loop over the leaf nodes and mark the outliers if the absolute value of the deviation is greater than the threshold
  graph.get_leaves().iter().for_each(|leaf| {
    let div = leaf.read_arc().payload().read().div;
    let date = leaf.read_arc().payload().read().date.unwrap();
    let was_outlier = leaf.read_arc().payload().write().is_outlier;
    let clock_deviation = clock_model.clock_deviation(date, div);
    let is_outlier = clock_deviation.abs() > iqd * clock_filter_threshold;
    if was_outlier != is_outlier {
      new_outliers += 1;
    }
    leaf.read_arc().payload().write().is_outlier = is_outlier;
  });
  new_outliers
}
