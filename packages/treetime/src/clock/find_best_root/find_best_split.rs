use crate::clock::clock_regression::ClockVarianceParams;
use crate::clock::clock_set::ClockSet;
use crate::clock::clock_state::{ClockInputs, ClockState};
use crate::clock::find_best_root::cost_function::BranchPointCostFunction;
use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RootObjective};
use crate::clock::find_best_root::{method_brent, method_golden_section, method_grid_search};
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;

#[derive(Debug, Serialize, Deserialize)]
pub struct FindRootResult {
  pub edge: Option<GraphEdgeKey>,

  pub split: f64,

  pub clock_set: ClockSet,

  pub chisq: f64,
}

pub(crate) fn find_best_split(
  graph: &Graph,
  inputs: &ClockInputs,
  state: &ClockState,
  edge: GraphEdgeKey,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &ClockVarianceParams,
  params: &BranchPointOptimizationParams,
  objective: RootObjective,
) -> Result<FindRootResult, Report> {
  let cost_fn = BranchPointCostFunction::new(graph, inputs, state, edge, branch_lengths, options, objective)?;

  match params {
    BranchPointOptimizationParams::Grid(params) => method_grid_search::optimize_grid_search(edge, &cost_fn, params),
    BranchPointOptimizationParams::Brent(params) => method_brent::optimize_brent(edge, &cost_fn, params),
    BranchPointOptimizationParams::GoldenSection(params) => {
      method_golden_section::optimize_golden_section(edge, &cost_fn, params)
    },
  }
}
