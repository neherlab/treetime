use crate::make_report;
use crate::reroot::cost_function::EdgeCostFn;
use crate::reroot::method_brent::optimize_brent;
use crate::reroot::params::BrentParams;
use crate::reroot::traits::RootStats;
use crate::reroot::variance::VarianceModel;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;

pub(crate) fn find_best_split<S>(
  graph: &Graph,
  edge: GraphEdgeKey,
  edge_stats: &BTreeMap<GraphEdgeKey, (S, S)>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  variance: &VarianceModel,
  params: &BrentParams,
) -> Result<FindRootResult, Report>
where
  S: RootStats,
{
  let edge_obj = graph
    .get_edge(edge)
    .ok_or_else(|| make_report!("Edge not found: {edge}"))?;
  let target_key = edge_obj.target();
  let branch_length = branch_lengths[&edge].ok_or_else(|| make_report!("Edge {edge} has no branch length"))?;

  let is_leaf = graph
    .get_node(target_key)
    .ok_or_else(|| make_report!("Target node not found for edge {edge}"))?
    .is_leaf();

  let (to_parent, to_child) = edge_stats
    .get(&edge)
    .ok_or_else(|| make_report!("No root statistics for edge {edge}"))?
    .clone();

  let cost_fn = EdgeCostFn {
    to_parent,
    to_child,
    branch_length,
    branch_variance: variance.branch(branch_length),
    is_leaf,
    leaf_time: None,
    variance_offset_leaf: variance.variance_offset_leaf,
  };

  optimize_brent(edge, &cost_fn, params)
}

#[derive(Debug, Clone)]
pub struct FindRootResult {
  pub(crate) edge: Option<GraphEdgeKey>,

  pub(crate) split: f64,

  pub(crate) score: f64,
}
