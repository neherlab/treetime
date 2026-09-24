use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

pub(crate) fn live_node_keys(graph: &Graph) -> BTreeSet<GraphNodeKey> {
  graph.get_nodes().map(|node| node.key()).collect()
}

pub(crate) fn reconcile_node_states<State>(
  node_states: BTreeMap<GraphNodeKey, State>,
  live_nodes: &BTreeSet<GraphNodeKey>,
  empty: impl Fn() -> State,
) -> BTreeMap<GraphNodeKey, State> {
  let mut node_states = node_states;
  for &key in live_nodes {
    node_states.entry(key).or_insert_with(&empty);
  }
  node_states.retain(|key, _| live_nodes.contains(key));
  node_states
}
