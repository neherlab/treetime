use crate::branch_lengths::invalid_branch_length_descriptions;
use crate::clock::clock_state::ClockInputs;
use crate::error::input_error;
use crate::progress::LogSink;
use crate::reroot::driver::reroot_at_best_root;
use crate::reroot::params::{RerootParams, RerootSpec};
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_graph::reroot::RerootResult;

pub(crate) fn reroot_date_free(
  graph: &mut Graph,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  spec: &RerootSpec,
  is_dated: &dyn Fn(GraphNodeKey) -> bool,
  advice: LengthsAdvice,
  names: &mut NodeNames,
  log: &dyn LogSink,
) -> Result<RerootResult, Report> {
  require_valid_lengths(graph, branch_lengths, spec, advice, names)?;
  let times = graph
    .get_nodes()
    .map(|node| {
      let key = node.key();
      (key, (graph.is_leaf(key) || is_dated(key)).then_some(0.0))
    })
    .collect();
  let mut inputs = ClockInputs::from_times(graph, &times, &BTreeMap::new());
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

#[derive(Clone, Copy, Debug)]
pub(crate) enum LengthsAdvice {
  OptimizeTheTree,
  InitialBranchLengths,
}

impl LengthsAdvice {
  fn text(self) -> &'static str {
    match self {
      Self::OptimizeTheTree => "Give the tree branch lengths, for example with `treetime optimize`",
      Self::InitialBranchLengths => {
        "Use --initial-branch-lengths=auto to fill in missing values or 'always' to recompute all branch lengths"
      },
    }
  }
}

fn require_valid_lengths(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  spec: &RerootSpec,
  advice: LengthsAdvice,
  names: &NodeNames,
) -> Result<(), Report> {
  let RerootSpec::Method(method) = spec else {
    return Ok(());
  };
  let invalid = invalid_branch_length_descriptions(graph, branch_lengths, names)?;
  if invalid.is_empty() {
    return Ok(());
  }
  Err(input_error(format!(
    "--reroot={method} places the root from branch lengths and requires every edge to have a finite, non-negative \
     branch length. Invalid edges:\n  {}\n{}",
    invalid.join("\n  "),
    advice.text()
  )))
}
