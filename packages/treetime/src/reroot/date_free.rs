use crate::clock::clock_state::ClockInputs;
use crate::progress::LogSink;
use crate::reroot::driver::reroot_at_best_root;
use crate::reroot::params::{RerootParams, RerootSpec};
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node_names::NodeNames;
use treetime_graph::reroot::RerootResult;

pub(crate) fn reroot_undated(
  graph: &mut Graph,
  names: &mut NodeNames,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  spec: &RerootSpec,
  log: &dyn LogSink,
) -> Result<RerootResult, Report> {
  let tip_dates = graph.get_leaves().map(|leaf| (leaf.key(), Some(0.0))).collect();
  let mut inputs = ClockInputs::from_times(graph, &tip_dates, &BTreeMap::new());
  let reroot_params = RerootParams::new(spec.clone(), false);
  reroot_at_best_root(
    graph,
    branch_lengths,
    &mut inputs,
    &BTreeSet::new(),
    &reroot_params,
    names,
    log,
  )
}
