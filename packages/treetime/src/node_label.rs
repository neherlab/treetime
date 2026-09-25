use std::collections::BTreeMap;
use treetime_graph::node::GraphNodeKey;

pub(crate) fn node_label(names: &BTreeMap<GraphNodeKey, Option<String>>, key: GraphNodeKey) -> String {
  names
    .get(&key)
    .cloned()
    .flatten()
    .unwrap_or_else(|| format!("node {key}"))
}
