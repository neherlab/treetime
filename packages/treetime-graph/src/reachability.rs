use crate::graph::Graph;
use crate::node::Node;
use parking_lot::RwLock;
use std::collections::BTreeSet;
use std::sync::Arc;

/// Whether a directed path runs from `start` to `finish` along edge directions (roots to leaves).
///
/// A node reaches itself, so `start == finish` is `true`. Serial depth-first walk over the outbound
/// edges; the visited set makes it terminate on cycles.
pub fn exists_forward_path_between<D>(graph: &Graph<D>, start: &Arc<RwLock<Node>>, finish: &Arc<RwLock<Node>>) -> bool
where
  D: Sync + Send,
{
  let finish_key = finish.read().key();
  let mut visited = BTreeSet::new();
  let mut stack = vec![Arc::clone(start)];

  while let Some(node) = stack.pop() {
    let key = node.read().key();
    if key == finish_key {
      return true;
    }
    if !visited.insert(key) {
      continue;
    }
    for (child, _) in graph.children_of(&node.read()) {
      stack.push(child);
    }
  }

  false
}
