use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

/// The keys of every node the graph currently holds.
pub fn live_node_keys(graph: &Graph) -> BTreeSet<GraphNodeKey> {
  graph.get_nodes().iter().map(|node| node.read_arc().key()).collect()
}

/// Reconcile a node-state map to a node set after a structural change: seed a placeholder state for
/// every live node absent from the map and drop the states of nodes that are gone. Leaf seeds and
/// other surviving states are carried across unchanged; the next marginal update recomputes the
/// evolving values.
pub fn reconcile_node_states<State>(
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
