use crate::graph::Graph;
use crate::node::GraphNodeKey;
use std::collections::BTreeSet;

/// Whether a directed path runs from `start` to `finish` along edge directions (roots to leaves).
///
/// A node reaches itself, so `start == finish` is `true`. Serial depth-first walk over the outbound
/// edges; the visited set makes it terminate on cycles.
pub fn exists_forward_path_between(graph: &Graph, start: GraphNodeKey, finish: GraphNodeKey) -> bool {
  let mut visited = BTreeSet::new();
  let mut stack = vec![start];

  while let Some(node_key) = stack.pop() {
    if node_key == finish {
      return true;
    }
    if !visited.insert(node_key) {
      continue;
    }
    if let Some(node) = graph.get_node(node_key) {
      for (child_key, _) in graph.children_keys_of(node) {
        stack.push(child_key);
      }
    }
  }

  false
}
