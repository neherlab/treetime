use crate::branch_lengths::is_valid_branch_length_value;
use crate::clock::clock_state::ClockInputs;
use crate::error::input_error;
use crate::reroot::find_root::EdgePoint;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::common_ancestor::common_ancestor;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::split_fraction::SplitFraction;

pub(crate) fn named_node_point(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  node_key: GraphNodeKey,
) -> Result<Option<EdgePoint>, Report> {
  let Some((parent_key, edge)) = graph.node_parent(node_key)? else {
    return Ok(None);
  };
  let middle = EdgePoint {
    edge,
    split: SplitFraction::MIDDLE,
  };
  if graph.node_parent(parent_key)?.is_some() {
    return Ok(Some(middle));
  }
  let [first, second] = graph[parent_key].outbound() else {
    return Ok(Some(middle));
  };
  let sibling = if *first == edge { *second } else { *first };
  let valid_length =
    |edge_key: GraphEdgeKey| branch_lengths[&edge_key].filter(|&length| is_valid_branch_length_value(length));
  let (Some(a), Some(b)) = (valid_length(edge), valid_length(sibling)) else {
    return Ok(None);
  };
  Ok(if a > b {
    Some(EdgePoint {
      edge,
      split: SplitFraction::new((a - b) / (2.0 * a))?,
    })
  } else if b > a {
    Some(EdgePoint {
      edge: sibling,
      split: SplitFraction::new((b - a) / (2.0 * b))?,
    })
  } else {
    None
  })
}

pub(crate) fn oldest_tip(graph: &Graph, inputs: &ClockInputs) -> Result<GraphNodeKey, Report> {
  graph
    .get_leaves()
    .filter_map(|node| {
      let key = node.key();
      let time = inputs.likely_time(key)?;
      Some((time, key))
    })
    .min_by(|(lhs, _), (rhs, _)| lhs.total_cmp(rhs))
    .map(|(_, key)| key)
    .ok_or_else(|| input_error("Cannot reroot to oldest tip because no dated leaves were found"))
}

pub(crate) fn tip_group_mrca(graph: &Graph, tips: &[GraphNodeKey]) -> Result<GraphNodeKey, Report> {
  if tips.is_empty() {
    return Err(input_error("--outgroup requires at least one tip name"));
  }
  common_ancestor(graph, tips)
}
