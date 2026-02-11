use crate::graph::edge::EdgeOptimizeOps;
use crate::graph::graph::Graph;
use crate::graph::node::NodeOptimizeOps;
use maplit::btreemap;
use std::collections::BTreeMap;

#[derive(Debug, Default, Copy, Clone)]
pub struct OnlyLeaves(pub bool);

/// Calculate mapping of node name to node divergence (accumulated by summing branch lengths)
pub fn compute_divs<N: NodeOptimizeOps, E: EdgeOptimizeOps, D: Send + Sync>(
  graph: &Graph<N, E, D>,
  only_leaves: OnlyLeaves,
) -> BTreeMap<String, f64> {
  let mut divs: BTreeMap<String, f64> = btreemap! {};
  let mut result: BTreeMap<String, f64> = btreemap! {};
  graph.iter_depth_first_preorder_forward(|node| {
    let name = node.payload.name().unwrap().as_ref().to_owned();
    if node.is_root {
      divs.insert(name.clone(), 0.0);
      if !only_leaves.0 {
        result.insert(name, 0.0);
      }
    } else {
      let (parent, edge) = node.get_exactly_one_parent().unwrap();
      let parent_div: f64 = {
        let parent_name = parent.read_arc().name().unwrap().as_ref().to_owned();
        divs.get(&parent_name).copied().unwrap_or_default()
      };
      let branch_length = edge.read_arc().branch_length().unwrap_or_default();
      let div = parent_div + branch_length;
      divs.insert(name.clone(), div);
      if node.is_leaf || !only_leaves.0 {
        result.insert(name, div);
      }
    }
  });
  result
}
