use std::collections::BTreeMap;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pair_by_name::pair_by_name;
use treetime_primitives::date::DateConstraint;
use treetime_utils::least_squares::LineFit;

pub(crate) fn half_residual_sum_of_squares(dates: &[f64], divs: &[f64]) -> f64 {
  let fit = LineFit::least_squares(dates, divs);
  0.5
    * dates
      .iter()
      .zip(divs)
      .map(|(date, div)| (div - fit.slope * date - fit.intercept).powi(2))
      .sum::<f64>()
}

pub(crate) fn dates_by_node(
  dates: impl IntoIterator<Item = (String, Option<DateConstraint>)>,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> BTreeMap<GraphNodeKey, DateConstraint> {
  pair_by_name(graph.get_nodes().map(|node| node.key()), names, dates)
    .by_node
    .into_iter()
    .filter_map(|(key, date)| Some((key, date?)))
    .collect()
}
