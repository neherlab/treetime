use crate::reroot::params::BrentParams;
use crate::reroot::split::{FindRootResult, find_best_split};
use crate::reroot::traits::RootStats;
use crate::reroot::variance::VarianceModel;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;

pub(crate) fn find_best_root<S>(
  graph: &Graph,
  edge_stats: &BTreeMap<GraphEdgeKey, (S, S)>,
  root_stats: &S,
  variance: &VarianceModel,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  params: &BrentParams,
) -> Result<FindRootResult<S>, Report>
where
  S: RootStats,
{
  let mut best = FindRootResult {
    edge: None,
    split: 0.0,
    stats: root_stats.clone(),
    score: root_stats.score(),
  };

  for edge_obj in graph.get_edges() {
    let edge_key = edge_obj.key();
    let res = find_best_split(graph, edge_key, edge_stats, branch_lengths, variance, params)?;
    if res.score < best.score {
      best = res;
    }
  }

  Ok(best)
}
