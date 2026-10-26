use crate::graph::Graph;
use crate::key_map::EdgeMap;
use crate::node::GraphNodeKey;
use eyre::Report;
use std::collections::BTreeSet;
use treetime_utils::{make_error, make_internal_error};

pub fn remove_leaves(
  graph: &mut Graph,
  branch_lengths: &mut EdgeMap<Option<f64>>,
  leaves: &BTreeSet<GraphNodeKey>,
) -> Result<(), Report> {
  leaves
    .iter()
    .try_for_each(|&leaf| remove_leaf(graph, branch_lengths, leaf))?;
  graph.build()
}

fn remove_leaf(graph: &mut Graph, branch_lengths: &mut EdgeMap<Option<f64>>, leaf: GraphNodeKey) -> Result<(), Report> {
  let Some(node) = graph.get_node(leaf) else {
    return make_internal_error!("Node {leaf} to remove is not in the graph");
  };
  if !node.is_leaf() {
    return make_internal_error!("Node {leaf} to remove is not a leaf");
  }
  let parent = match node.inbound() {
    [] => return make_error!("Removing the tips would leave an empty tree"),
    [edge] => graph[*edge].source(),
    _ => return make_internal_error!("Leaf {leaf} to remove has more than one parent"),
  };
  remove_node_and_edges(graph, branch_lengths, leaf)?;
  remove_trivial_node(graph, branch_lengths, parent)
}

fn remove_trivial_node(
  graph: &mut Graph,
  branch_lengths: &mut EdgeMap<Option<f64>>,
  node_key: GraphNodeKey,
) -> Result<(), Report> {
  let node = &graph[node_key];
  match (node.inbound(), node.outbound()) {
    ([_], []) | ([], []) => remove_leaf(graph, branch_lengths, node_key),
    ([parent_edge], [child_edge]) => {
      let (parent_edge, child_edge) = (*parent_edge, *child_edge);
      let grandparent = graph[parent_edge].source();
      branch_lengths.insert(
        child_edge,
        merged_length(branch_lengths[&parent_edge], branch_lengths[&child_edge]),
      );
      graph.reparent_edge(child_edge, grandparent)?;
      remove_node_and_edges(graph, branch_lengths, node_key)
    },
    ([], [_]) => remove_node_and_edges(graph, branch_lengths, node_key),
    _ => Ok(()),
  }
}

fn remove_node_and_edges(
  graph: &mut Graph,
  branch_lengths: &mut EdgeMap<Option<f64>>,
  node_key: GraphNodeKey,
) -> Result<(), Report> {
  let (_, edges) = graph.remove_node(node_key)?;
  for edge in &edges {
    branch_lengths.remove(edge.key());
  }
  Ok(())
}

fn merged_length(parent: Option<f64>, child: Option<f64>) -> Option<f64> {
  match (parent, child) {
    (Some(parent), Some(child)) => Some(parent + child),
    (parent, child) => parent.or(child),
  }
}
