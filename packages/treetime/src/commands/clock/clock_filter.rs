use crate::commands::clock::clock_graph::GraphClock;
use crate::commands::clock::clock_model::ClockModel;
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::edge::HasBranchLength;
use crate::graph::node::Outlier;
use itertools::Itertools;
use ordered_float::OrderedFloat;

/// Get results of the root-to-tip clock inference.
#[allow(clippy::integer_division_remainder_used)]
pub fn clock_filter_inplace(graph: &GraphClock, clock_model: &ClockModel, clock_filter_threshold: f64) -> i32 {
  log::info!("### Filtering outliers (threshold={clock_filter_threshold})");
  log::debug!(
    "Clock model for filtering: rate={:.6e}, intercept={:.4}",
    clock_model.clock_rate(),
    clock_model.intercept()
  );

  // assign divergence to each node
  graph.par_iter_breadth_first_forward(|mut node| {
    node.payload.div = node
      .get_exactly_one_parent()
      .map(|(parent, edge)| {
        let branch_length = edge.read_arc().branch_length().unwrap_or_default();
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
    .filter_map(|leaf| {
      let div = leaf.read_arc().payload().read().div;
      let date = leaf.read_arc().payload().read().date;
      date.map(|date| clock_model.clock_deviation(date, div))
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
    let date = leaf.read_arc().payload().read().date;
    let was_outlier = leaf.read_arc().payload().read().is_outlier();
    if let Some(date) = date {
      let clock_deviation = clock_model.clock_deviation(date, div);
      let is_outlier = clock_deviation.abs() > iqd * clock_filter_threshold;
      if was_outlier != is_outlier {
        new_outliers += 1;
      }
      leaf.read_arc().payload().write().set_is_outlier(is_outlier);
    }
  });

  log::info!("Outlier filtering: {new_outliers} leaves changed status, IQD={iqd:.6e}");
  log::debug!(
    "Leaf clock deviations (n={}): min={:.6e}, Q1={:.6e}, Q3={:.6e}, max={:.6e}",
    n,
    leaf_clock_deviations.first().copied().unwrap_or(0.0),
    leaf_clock_deviations[iq25],
    leaf_clock_deviations[iq75],
    leaf_clock_deviations.last().copied().unwrap_or(0.0)
  );

  new_outliers
}
