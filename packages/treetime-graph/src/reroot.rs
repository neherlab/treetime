use crate::edge::{GraphEdge, GraphEdgeKey, HasBranchLength, invert_edge};
use crate::graph::Graph;
use crate::node::{GraphNode, GraphNodeKey};
use eyre::Report;
use serde::{Deserialize, Serialize};

/// Information about an edge split during reroot.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EdgeSplitInfo {
  /// The original edge that was split.
  pub old_edge_key: GraphEdgeKey,
  /// The new node created at the split point.
  pub new_node_key: GraphNodeKey,
  /// The edge from the original source to the new node.
  pub parent_side_edge_key: GraphEdgeKey,
  /// The edge from the new node to the original target.
  pub child_side_edge_key: GraphEdgeKey,
  /// Position along the edge where the split occurred (0.0 = target, 1.0 = source).
  pub split_position: f64,
}

/// Information about an edge merge when removing a trivial node.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EdgeMergeInfo {
  /// The node that was removed.
  pub removed_node_key: GraphNodeKey,
  /// The edge from parent to the removed node.
  pub parent_edge_key: GraphEdgeKey,
  /// The edge from the removed node to its child.
  pub child_edge_key: GraphEdgeKey,
  /// The new merged edge from parent to child.
  pub merged_edge_key: GraphEdgeKey,
}

/// Result of a reroot operation.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RerootResult {
  /// The key of the new root node.
  pub new_root_key: GraphNodeKey,
  /// Information about edge split, if one occurred.
  pub edge_split: Option<EdgeSplitInfo>,
  /// Information about edge merge, if a trivial node was removed.
  pub edge_merge: Option<EdgeMergeInfo>,
  /// Keys of edges whose direction was inverted during rerooting.
  /// Empty when root did not change.
  pub inverted_edge_keys: Vec<GraphEdgeKey>,
}

/// Bundles all topology changes from a reroot operation for partition updates.
///
/// Passed to `PartitionRerootOps::apply_reroot` to update partition state in a single call.
#[derive(Clone, Debug, Default)]
pub struct RerootChanges {
  /// Edge split info if a new node was created at the reroot point.
  pub edge_split: Option<EdgeSplitInfo>,
  /// Edge merge info if the old root was removed as a trivial node.
  pub edge_merge: Option<EdgeMergeInfo>,
  /// Keys of edges on the path from old root to new root (post-inversion direction).
  /// Empty if root did not change or old root was removed.
  pub inverted_edge_keys: Vec<GraphEdgeKey>,
}

/// Split an edge by inserting a new node at `split_position` along it.
///
/// The original edge is removed and replaced with two new edges:
/// - parent-side: from original source to new node (length = `split_position * original_length`)
/// - child-side: from new node to original target (length = `(1 - split_position) * original_length`)
///
/// The new node is created with `N::default()`. Callers needing domain-specific
/// initialization (e.g. clock data) should configure the returned node afterward.
pub fn split_edge<N, E, D>(
  graph: &mut Graph<N, E, D>,
  edge_key: GraphEdgeKey,
  split_position: f64,
) -> Result<EdgeSplitInfo, Report>
where
  N: GraphNode + Default,
  E: GraphEdge + HasBranchLength + Default,
  D: Send + Sync,
{
  let new_node_key = graph.add_node(N::default());

  let (source_key, target_key, branch_length) = {
    let edge = graph.get_edge(edge_key).expect("Edge not found");
    let source_key = edge.read_arc().source();
    let target_key = edge.read_arc().target();
    let branch_length = edge.read_arc().payload().read_arc().branch_length().unwrap_or_default();
    (source_key, target_key, branch_length)
  };

  let mut parent_edge = E::default();
  parent_edge.set_branch_length(Some(split_position * branch_length));
  let parent_side_edge_key = graph.add_edge(source_key, new_node_key, parent_edge)?;

  let mut child_edge = E::default();
  child_edge.set_branch_length(Some((1.0 - split_position) * branch_length));
  let child_side_edge_key = graph.add_edge(new_node_key, target_key, child_edge)?;

  graph.remove_edge(edge_key)?;

  Ok(EdgeSplitInfo {
    old_edge_key: edge_key,
    new_node_key,
    parent_side_edge_key,
    child_side_edge_key,
    split_position,
  })
}

/// Invert edges along the path from old root to new root.
///
/// Returns the keys of all inverted edges (in old-root-to-new-root order).
/// Only inverts graph topology (edge direction). Domain-specific edge data
/// (clock messages, partition state) must be updated by the caller.
pub fn apply_reroot_topology<N, E, D>(
  graph: &mut Graph<N, E, D>,
  old_root_key: GraphNodeKey,
  new_root_key: GraphNodeKey,
) -> Result<Vec<GraphEdgeKey>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Send + Sync,
{
  let paths = graph.path_from_node_to_node(new_root_key, old_root_key)?;

  let mut inverted_edge_keys = Vec::new();
  for (_, edge) in &paths {
    if let Some(edge) = edge {
      inverted_edge_keys.push(edge.read_arc().key());
      invert_edge(graph, edge);
    }
  }

  graph.build()?;
  Ok(inverted_edge_keys)
}

/// Remove a node if it is trivial (exactly one parent and one child), merging the edges.
///
/// Returns `Some(EdgeMergeInfo)` if the node was removed, `None` if the node was not trivial.
/// The merged edge receives the sum of branch lengths. Domain-specific properties
/// of the removed node and its connecting edges are discarded.
pub fn remove_node_if_trivial<N, E, D>(
  graph: &mut Graph<N, E, D>,
  node_key: GraphNodeKey,
) -> Result<Option<EdgeMergeInfo>, Report>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength + Default,
  D: Send + Sync,
{
  let (parent_edge_key, child_edge_key) = {
    let node = graph.get_node(node_key).expect("Node not found");
    let node = node.read_arc();
    if node.inbound().len() != 1 || node.outbound().len() != 1 {
      return Ok(None);
    }
    (node.inbound()[0], node.outbound()[0])
  };

  let (parent_key, parent_branch) = {
    let parent_edge = graph.get_edge(parent_edge_key).expect("Parent edge not found");
    let parent_edge = parent_edge.read_arc();
    (parent_edge.source(), parent_edge.payload().read_arc().branch_length())
  };

  let (child_key, child_branch) = {
    let child_edge = graph.get_edge(child_edge_key).expect("Child edge not found");
    let child_edge = child_edge.read_arc();
    (child_edge.target(), child_edge.payload().read_arc().branch_length())
  };

  let merged_branch_length = match (parent_branch, child_branch) {
    (Some(a), Some(b)) => Some(a + b),
    (a, b) => a.or(b),
  };

  graph.remove_node(node_key)?;

  let mut merged_payload = E::default();
  merged_payload.set_branch_length(merged_branch_length);
  let merged_edge_key = graph.add_edge(parent_key, child_key, merged_payload)?;

  graph.build()?;

  Ok(Some(EdgeMergeInfo {
    removed_node_key: node_key,
    parent_edge_key,
    child_edge_key,
    merged_edge_key,
  }))
}
