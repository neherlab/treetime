use crate::commands::clock::clock_regression::ClockOptions;
use crate::commands::clock::clock_set::ClockSet;
use crate::commands::clock::clock_traits::{ClockEdge, ClockNode};
use crate::commands::clock::find_best_root::find_best_root::find_best_root;
use crate::commands::clock::find_best_root::find_best_split::FindRootResult;
use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::graph::edge::{GraphEdge, GraphEdgeKey, invert_edge};
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, GraphNodeKey};
use approx::ulps_eq;
use eyre::Report;

pub fn reroot_in_place<N, E, D>(
  graph: &mut Graph<N, E, D>,
  options: &ClockOptions,
  params: &BranchPointOptimizationParams,
) -> Result<GraphNodeKey, Report>
where
  N: GraphNode + ClockNode + Default,
  E: GraphEdge + ClockEdge + Default,
  D: Send + Sync,
{
  let FindRootResult { edge, split, total, .. } = find_best_root(graph, options, params)?;

  let old_root_key = { graph.get_exactly_one_root()?.read_arc().key() };
  let Some(edge_key) = edge else {
    return Ok(old_root_key); // Already at the best root
  };

  // Extract edge endpoints before potentially removing the edge
  let (source_key, target_key) = {
    let edge = graph.get_edge(edge_key).expect("Edge not found");
    (edge.read_arc().source(), edge.read_arc().target())
  };

  // Determine where to place the new root based on the split position along the edge
  // - split = 0.0 means root at target
  // - split = 1.0 means root at source
  // - 0.0 < split < 1.0 means somewhere in the middle - i.e. create a new node
  let new_root_key = if ulps_eq!(split, 0.0, max_ulps = 5) {
    target_key
  } else if ulps_eq!(split, 1.0, max_ulps = 5) {
    source_key
  } else {
    create_new_root_node(graph, edge_key, split, total)?
  };

  if new_root_key != old_root_key {
    apply_reroot(graph, old_root_key, new_root_key, options)?;

    // Clean up old root if it is now a trivial node, i.e. has exactly one parent and one child
    // Nb: the attributes of the old root node and branches (other than branch lengths) are discarded
    remove_node_if_trivial(graph, old_root_key)?;
  }

  Ok(new_root_key)
}

/// Create new root node by splitting the edge into two
fn create_new_root_node<N, E, D>(
  graph: &mut Graph<N, E, D>,
  edge_key: GraphEdgeKey,
  split: f64,
  total: ClockSet,
) -> Result<GraphNodeKey, Report>
where
  N: GraphNode + ClockNode + Default,
  E: GraphEdge + ClockEdge + Default,
  D: Send + Sync,
{
  let mut new_node = N::default();
  new_node.set_clock_set(total);
  let new_root_key = graph.add_node(new_node);

  // Extract edge data before removing the edge to avoid Arc reference issues
  let (source_key, target_key, branch_length) = {
    let edge = graph.get_edge(edge_key).expect("Edge not found");
    let source_key = edge.read_arc().source();
    let target_key = edge.read_arc().target();
    let branch_length = edge.read_arc().payload().read_arc().branch_length().unwrap_or_default();
    (source_key, target_key, branch_length)
  };

  let mut edge1 = E::default();
  edge1.set_branch_length(Some(split * branch_length));
  graph.add_edge(source_key, new_root_key, edge1)?;

  let mut edge2 = E::default();
  edge2.set_branch_length(Some((1.0 - split) * branch_length));
  graph.add_edge(new_root_key, target_key, edge2)?;

  graph.remove_edge(edge_key)?;

  Ok(new_root_key)
}

/// Modify graph topology to make the newly identified root the actual root.
fn apply_reroot<N, E, D>(
  graph: &mut Graph<N, E, D>,
  old_root_key: GraphNodeKey,
  new_root_key: GraphNodeKey,
  options: &ClockOptions,
) -> Result<(), Report>
where
  N: GraphNode + ClockNode,
  E: GraphEdge + ClockEdge,
  D: Send + Sync,
{
  // Find paths from the old root to the new desired root
  let paths = graph.path_from_node_to_node(new_root_key, old_root_key)?;

  // Invert every edge on the path from old to new root.
  // This will make the desired new root into an actual root. The old root might no longer be a root.
  for (_, edge) in &paths {
    if let Some(edge) = edge {
      invert_edge(graph, edge);
      let mut edge_payload = edge.read_arc().payload().write_arc();
      let edge_len = edge_payload.branch_length().unwrap();
      let branch_variance = options.variance_factor * edge_len + options.variance_offset;
      let tmp_to_parent = edge_payload.to_parent().clone();
      *edge_payload.to_parent_mut() = edge_payload.to_child().clone();
      *edge_payload.to_child_mut() = tmp_to_parent;
      *edge_payload.from_child_mut() = edge_payload.to_parent().propagate_averages(edge_len, branch_variance);
    }
  }

  // Some bookkeeping
  graph.build()?;
  Ok(())
}

/// Remove a node if it is trivial (i.e. has exactly one parent and one child), merging the edges
/// lengths.
/// Nb: when applied this function removes properties of the removed node and its connecting edges.
fn remove_node_if_trivial<N, E, D>(graph: &mut Graph<N, E, D>, node_key: GraphNodeKey) -> Result<(), Report>
where
  N: GraphNode + ClockNode,
  E: GraphEdge + ClockEdge + Default,
  D: Send + Sync,
{
  let (parent_edge_key, child_edge_key) = {
    let node = graph.get_node(node_key).expect("Node not found");
    let node = node.read_arc();
    if node.inbound().len() != 1 || node.outbound().len() != 1 {
      return Ok(());
    }
    (node.inbound()[0], node.outbound()[0])
  };

  let (parent_key, parent_branch) = {
    let parent_edge = graph.get_edge(parent_edge_key).expect("Parent edge not found");
    let parent_edge = parent_edge.read_arc();
    (
      parent_edge.source(),
      parent_edge.payload().read_arc().branch_length(),
    )
  };

  let (child_key, child_branch) = {
    let child_edge = graph.get_edge(child_edge_key).expect("Child edge not found");
    let child_edge = child_edge.read_arc();
    (
      child_edge.target(),
      child_edge.payload().read_arc().branch_length(),
    )
  };

  let merged_branch_length = match (parent_branch, child_branch) {
    (Some(parent_len), Some(child_len)) => Some(parent_len + child_len),
    (Some(len), None) | (None, Some(len)) => Some(len),
    (None, None) => None,
  };

  let mut merged_payload = E::default();
  merged_payload.set_branch_length(merged_branch_length);

  graph.remove_edge(parent_edge_key)?;
  graph.remove_edge(child_edge_key)?;
  graph.remove_node(node_key)?;

  graph.add_edge(parent_key, child_key, merged_payload)?;

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::commands::clock::clock_graph::ClockGraph;
  use crate::io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};

  #[test]
  fn test_remove_node_if_trivial_simple() -> Result<(), Report> {
    // define tree with trivial node:
    //        root
    //        /  \
    //      mid  B
    //      /
    //     A
    let mut graph: ClockGraph = nwk_read_str("((A:0.5)mid:0.3,B:0.2)root;")?;

    let mid_key = graph
      .get_nodes()
      .iter()
      .find_map(|node| {
        let node_guard = node.read_arc();
        let payload = node_guard.payload().read_arc();
        (payload.name.as_deref() == Some("mid")).then_some(node_guard.key())
      })
      .expect("Expected node named 'mid'");

    remove_node_if_trivial(&mut graph, mid_key)?;

    assert!(graph.get_node(mid_key).is_none(), "Expected node to be removed");

    let expected_newick = "(B:0.2[&mutations=\"\"],A:0.8[&mutations=\"\"])root[&mutations=\"\"];";
    let resulting_newick = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(
      resulting_newick, expected_newick,
      "Unexpected resulting tree structure after removing trivial node"
    );

    Ok(())
  }
}
