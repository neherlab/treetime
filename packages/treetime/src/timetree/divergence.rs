use crate::clock::divergence::root_to_node_divergences;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

pub(crate) fn final_divergences(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  filter_divergences: Option<&BTreeMap<GraphNodeKey, f64>>,
) -> Result<BTreeMap<GraphNodeKey, f64>, Report> {
  let divergences = root_to_node_divergences(graph, |edge_key| branch_lengths[&edge_key].unwrap_or_default())?;
  Ok(
    graph
      .get_nodes()
      .map(|node| {
        let key = node.key();
        let div = if names[&key].is_some() {
          divergences[&key]
        } else {
          filter_divergences
            .and_then(|filter_divergences| filter_divergences.get(&key))
            .copied()
            .unwrap_or(0.0)
        };
        (key, div)
      })
      .collect(),
  )
}
