use crate::make_report;
use crate::reroot::params::BrentParams;
use crate::reroot::search::find_best_root;
use crate::reroot::traits::RootStats;
use crate::reroot::variance::VarianceModel;
use approx::ulps_eq;
use eyre::Report;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::reroot::{
  RerootResult, apply_reroot_topology, record_merge, record_split, remove_node_if_trivial, split_edge,
  trivial_node_branch_lengths,
};

pub(crate) fn reroot_in_place<S, F>(
  graph: &mut Graph,
  edge_stats: &BTreeMap<GraphEdgeKey, (S, S)>,
  root_stats: &S,
  variance: &VarianceModel,
  opt_params: &BrentParams,
  topo: RerootTopologyParams,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  fixup: F,
) -> Result<RerootResult, Report>
where
  S: RootStats,
  F: FnMut(&mut Graph, &[GraphEdgeKey]) -> Result<(), Report>,
{
  let best = find_best_root(graph, edge_stats, root_stats, variance, branch_lengths, opt_params)?;
  apply_root_at_edge(graph, best.edge, best.split, topo, branch_lengths, fixup)
}

pub(crate) fn reroot_at_node<F>(
  graph: &mut Graph,
  node_key: GraphNodeKey,
  topo: RerootTopologyParams,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  fixup: F,
) -> Result<RerootResult, Report>
where
  F: FnMut(&mut Graph, &[GraphEdgeKey]) -> Result<(), Report>,
{
  let edge = graph.parent_inbound_edge(node_key)?;
  apply_root_at_edge(graph, edge, 0.5, topo, branch_lengths, fixup)
}

fn apply_root_at_edge<F>(
  graph: &mut Graph,
  edge: Option<GraphEdgeKey>,
  split: f64,
  topo: RerootTopologyParams,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  mut fixup: F,
) -> Result<RerootResult, Report>
where
  F: FnMut(&mut Graph, &[GraphEdgeKey]) -> Result<(), Report>,
{
  let old_root_key = graph.get_exactly_one_root()?.key();

  let Some(edge_key) = edge else {
    return Ok(RerootResult {
      new_root_key: old_root_key,
      edge_split: None,
      edge_merge: None,
      inverted_edge_keys: vec![],
    });
  };

  let (source_key, target_key) = {
    let edge = graph
      .get_edge(edge_key)
      .ok_or_else(|| make_report!("Edge not found: {edge_key}"))?;
    (edge.source(), edge.target())
  };

  let (new_root_key, edge_split) = if ulps_eq!(split, 0.0, max_ulps = 5) {
    (source_key, None)
  } else if ulps_eq!(split, 1.0, max_ulps = 5) {
    (target_key, None)
  } else if topo.split_edge {
    let length = branch_lengths.get(&edge_key).copied().flatten();
    let info = split_edge(graph, edge_key, split, length)?;
    record_split(branch_lengths, &info);
    (info.new_node_key, Some(info))
  } else {
    (if split < 0.5 { source_key } else { target_key }, None)
  };

  let (inverted_edge_keys, edge_merge) = if new_root_key != old_root_key {
    let mut inverted = apply_reroot_topology(graph, old_root_key, new_root_key)?;
    fixup(graph, &inverted)?;

    let merge = if topo.remove_trivial_root {
      let (parent_branch, child_branch) = trivial_node_branch_lengths(graph, old_root_key, branch_lengths);
      remove_node_if_trivial(graph, old_root_key, parent_branch, child_branch)?
    } else {
      None
    };

    if let Some(merge) = &merge {
      inverted.retain(|k| *k != merge.parent_edge_key && *k != merge.child_edge_key);
      record_merge(branch_lengths, merge);
    }

    (inverted, merge)
  } else {
    (vec![], None)
  };

  Ok(RerootResult {
    new_root_key,
    edge_split,
    edge_merge,
    inverted_edge_keys,
  })
}

#[derive(Debug, Clone, Copy, SmartDefault, Serialize, Deserialize)]
pub struct RerootTopologyParams {
  #[default = true]
  pub(crate) split_edge: bool,

  #[default = true]
  pub(crate) remove_trivial_root: bool,
}
