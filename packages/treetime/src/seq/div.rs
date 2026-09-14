use crate::seq::mutation::Sub;
use eyre::Report;
use maplit::btreemap;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::collections::container::get_exactly_one;

#[derive(Debug, Default, Copy, Clone)]
pub struct OnlyLeaves(pub bool);

/// Calculate mapping of node name to node divergence (accumulated by summing branch lengths).
/// Only nodes with names are included in the result.
pub fn compute_divs(
  graph: &Graph,
  only_leaves: OnlyLeaves,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<BTreeMap<String, f64>, Report> {
  // Track divergence by node key (always available) for internal computation
  let mut divs_by_key: BTreeMap<GraphNodeKey, f64> = btreemap! {};
  let mut result: BTreeMap<String, f64> = btreemap! {};

  graph.iter_depth_first_preorder_forward(|node| {
    let div = if node.is_root {
      0.0
    } else {
      let (parent_key, edge_key) = get_exactly_one(&node.parent_keys).unwrap();
      let parent_div = divs_by_key.get(parent_key).copied().unwrap_or_default();
      let branch_length = branch_lengths[edge_key].unwrap_or_default();
      parent_div + branch_length
    };

    divs_by_key.insert(node.key, div);

    if node.is_leaf || !only_leaves.0 {
      if let Some(name) = &names[&node.key] {
        result.insert(name.clone(), div);
      }
    }
    Ok(())
  })?;

  Ok(result)
}

/// Count reconstructed substitutions per edge.
///
/// Returns a map from edge key to the number of canonical (non-gap, non-ambiguous) substitutions on
/// that edge, read from a pre-gathered per-edge substitution map.
pub fn compute_edge_mutation_counts(
  graph: &Graph,
  edge_subs: &BTreeMap<GraphEdgeKey, Vec<Sub>>,
) -> BTreeMap<GraphEdgeKey, usize> {
  graph
    .get_edges()
    .map(|edge| {
      let edge_key = edge.key();
      (edge_key, edge_subs[&edge_key].len())
    })
    .collect()
}
