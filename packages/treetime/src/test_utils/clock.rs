use crate::clock::clock_regression::{ClockFitParams, ClockFitResult, ClockTree, fit_clock_model};
use crate::progress::NoopProgress;
use crate::reroot::driver::reroot_at_best_root;
use crate::reroot::params::RerootParams;
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_graph::pair_by_name::pair_by_name;
use treetime_graph::reroot::RerootResult;
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
  names: &NodeNames,
) -> BTreeMap<GraphNodeKey, DateConstraint> {
  pair_by_name(graph.get_nodes().map(|node| node.key()), names, dates)
    .by_node
    .into_iter()
    .filter_map(|(key, date)| Some((key, date?)))
    .collect()
}

pub(crate) fn reroot_and_fit_clock(
  tree: &mut ClockTree,
  outliers: &BTreeSet<GraphNodeKey>,
  clock_rate: Option<f64>,
  reroot_params: &RerootParams,
  names: &mut NodeNames,
) -> Result<(RerootResult, ClockFitResult), Report> {
  let ClockTree {
    graph,
    branch_lengths,
    inputs,
  } = tree;
  let reroot = reroot_at_best_root(
    graph,
    branch_lengths,
    inputs,
    outliers,
    &reroot_params.clone().with_clock_rate(clock_rate),
    names,
    &NoopProgress,
  )?;
  let fit = fit_clock_model(
    tree,
    &ClockFitParams {
      outliers,
      clock_rate,
      prev_clock_rate: None,
    },
    &NoopProgress,
  )?;
  Ok((reroot, fit))
}
