use maplit::btreemap;
use std::collections::BTreeMap;
use treetime_graph::edge::EdgeOptimizeOps;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNodeKey, NodeOptimizeOps};
use treetime_utils::collections::container::get_exactly_one;

#[derive(Debug, Default, Copy, Clone)]
pub struct OnlyLeaves(pub bool);

/// Calculate mapping of node name to node divergence (accumulated by summing branch lengths).
/// Only nodes with names are included in the result.
pub fn compute_divs<N: NodeOptimizeOps, E: EdgeOptimizeOps, D: Send + Sync>(
  graph: &Graph<N, E, D>,
  only_leaves: OnlyLeaves,
) -> BTreeMap<String, f64> {
  // Track divergence by node key (always available) for internal computation
  let mut divs_by_key: BTreeMap<GraphNodeKey, f64> = btreemap! {};
  let mut result: BTreeMap<String, f64> = btreemap! {};

  graph.iter_depth_first_preorder_forward(|node| {
    let div = if node.is_root {
      0.0
    } else {
      let (parent_key, edge_key) = get_exactly_one(&node.parent_keys).unwrap();
      let parent_div = divs_by_key.get(parent_key).copied().unwrap_or_default();
      let edge = graph.get_edge(*edge_key).unwrap();
      let branch_length = edge.read_arc().payload().read_arc().branch_length().unwrap_or_default();
      parent_div + branch_length
    };

    divs_by_key.insert(node.key, div);

    // Add to result only if node has a name and matches filter criteria
    if node.is_leaf || !only_leaves.0 {
      if let Some(name) = node.payload.name() {
        result.insert(name.as_ref().to_owned(), div);
      }
    }
  });

  result
}
