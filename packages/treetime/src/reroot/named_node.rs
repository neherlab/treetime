use crate::clock::clock_state::ClockInputs;
use crate::error::input_error;
use crate::reroot::find_root::EdgePoint;
use eyre::Report;
use treetime_graph::common_ancestor::common_ancestor;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::split_fraction::SplitFraction;

pub(crate) fn named_node_point(graph: &Graph, node_key: GraphNodeKey) -> Result<Option<EdgePoint>, Report> {
  Ok(graph.node_parent(node_key)?.map(|(_, edge)| EdgePoint {
    edge,
    split: SplitFraction::MIDDLE,
  }))
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
