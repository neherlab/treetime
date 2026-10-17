use crate::clock::clock_regression::RegressedClockTree;
use crate::clock::clock_set::ClockSet;
use crate::clock::find_best_root::cost_function::BranchPointCostFunction;
use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RootObjective};
use crate::clock::find_best_root::{method_brent, method_golden_section, method_grid_search};
use crate::node_label::node_label;
use crate::progress::LogSink;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;

pub(crate) fn find_best_split(
  regressed: RegressedClockTree<'_>,
  edge: GraphEdgeKey,
  params: &BranchPointOptimizationParams,
  objective: RootObjective,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> Result<FindRootResult, Report> {
  let cost_fn = BranchPointCostFunction::new(regressed, edge, objective)?;
  let branch = node_label(names, regressed.tree.graph.get_target_node_key(edge)?);

  match params {
    BranchPointOptimizationParams::Grid(params) => {
      method_grid_search::optimize_grid_search(edge, &cost_fn, params, &branch, log)
    },
    BranchPointOptimizationParams::Brent(params) => method_brent::optimize_brent(edge, &cost_fn, params, &branch, log),
    BranchPointOptimizationParams::GoldenSection(params) => {
      method_golden_section::optimize_golden_section(edge, &cost_fn, params, &branch, log)
    },
  }
}

#[derive(Debug)]
pub struct FindRootResult {
  pub(crate) edge: Option<GraphEdgeKey>,

  pub(crate) split: f64,

  pub(crate) clock_set: ClockSet,

  pub(crate) chisq: f64,
}
