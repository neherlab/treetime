use crate::edge::{GraphEdgeKey, invert_edge};
use crate::graph::Graph;
use crate::node::GraphNodeKey;
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

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
  /// Branch length of the parent-side edge (`split_position * original_length`).
  pub parent_side_length: Option<f64>,
  /// Branch length of the child-side edge (`(1 - split_position) * original_length`).
  pub child_side_length: Option<f64>,
  /// Position along the edge where the split occurred (0.0 = source, 1.0 = target).
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
  /// Branch length of the merged edge: the sum when both sides carry a length, otherwise whichever
  /// side has one, and `None` when neither does.
  pub merged_branch_length: Option<f64>,
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
/// Passed to each partition's reroot pass to update partition state in a single call.
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
/// The original branch length is supplied by the caller (`branch_length`, taken from its
/// branch-length value map). The two new edges carry no branch length; their lengths are returned in
/// [`EdgeSplitInfo`] for the caller to record in its map. A missing input length resolves to `0.0`
/// for the split, matching the input-tree derivation.
pub fn split_edge(
  graph: &mut Graph,
  edge_key: GraphEdgeKey,
  split_position: f64,
  branch_length: Option<f64>,
) -> Result<EdgeSplitInfo, Report> {
  let new_node_key = graph.add_node();

  let (source_key, target_key) = {
    let edge = graph.get_edge(edge_key).expect("Edge not found");
    let source_key = edge.read_arc().source();
    let target_key = edge.read_arc().target();
    (source_key, target_key)
  };

  let length = branch_length.unwrap_or_default();
  let parent_side_length = Some(split_position * length);
  let child_side_length = Some((1.0 - split_position) * length);

  let parent_side_edge_key = graph.add_edge(source_key, new_node_key)?;
  let child_side_edge_key = graph.add_edge(new_node_key, target_key)?;

  graph.remove_edge(edge_key)?;

  Ok(EdgeSplitInfo {
    old_edge_key: edge_key,
    new_node_key,
    parent_side_edge_key,
    child_side_edge_key,
    parent_side_length,
    child_side_length,
    split_position,
  })
}

/// Invert edges along the path from old root to new root.
///
/// Returns the keys of all inverted edges (in old-root-to-new-root order).
/// Only inverts graph topology (edge direction). Domain-specific edge data
/// (clock messages, partition state) must be updated by the caller.
pub fn apply_reroot_topology(
  graph: &mut Graph,
  old_root_key: GraphNodeKey,
  new_root_key: GraphNodeKey,
) -> Result<Vec<GraphEdgeKey>, Report> {
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
///
/// The parent-side and child-side branch lengths are supplied by the caller (from its branch-length
/// value map). The merged edge carries no branch length; the merged length is returned in
/// [`EdgeMergeInfo`] for the caller to record in its map. Merge semantics: the sum when both sides
/// carry a length, otherwise whichever side has one (never coerce `None` to `0.0` and sum).
pub fn remove_node_if_trivial(
  graph: &mut Graph,
  node_key: GraphNodeKey,
  parent_branch: Option<f64>,
  child_branch: Option<f64>,
) -> Result<Option<EdgeMergeInfo>, Report> {
  let (parent_edge_key, child_edge_key) = {
    let node = graph.get_node(node_key).expect("Node not found");
    let node = node.read_arc();
    if node.inbound().len() != 1 || node.outbound().len() != 1 {
      return Ok(None);
    }
    (node.inbound()[0], node.outbound()[0])
  };

  let parent_key = {
    let parent_edge = graph.get_edge(parent_edge_key).expect("Parent edge not found");
    parent_edge.read_arc().source()
  };

  let child_key = {
    let child_edge = graph.get_edge(child_edge_key).expect("Child edge not found");
    child_edge.read_arc().target()
  };

  let merged_branch_length = match (parent_branch, child_branch) {
    (Some(a), Some(b)) => Some(a + b),
    (a, b) => a.or(b),
  };

  graph.remove_node(node_key)?;

  let merged_edge_key = graph.add_edge(parent_key, child_key)?;

  graph.build()?;

  Ok(Some(EdgeMergeInfo {
    removed_node_key: node_key,
    parent_edge_key,
    child_edge_key,
    merged_edge_key,
    merged_branch_length,
  }))
}

/// The branch lengths of a trivial node's parent-side and child-side edges, read from a value map.
///
/// Returns `(None, None)` when the node is not trivial (not exactly one inbound and one outbound
/// edge), matching the guard in [`remove_node_if_trivial`]. A caller can therefore compute the merge
/// inputs unconditionally and leave the triviality decision to the removal.
pub fn trivial_node_branch_lengths(
  graph: &Graph,
  node_key: GraphNodeKey,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> (Option<f64>, Option<f64>) {
  let node = graph.get_node(node_key).expect("Node not found");
  let node = node.read_arc();
  if node.inbound().len() != 1 || node.outbound().len() != 1 {
    return (None, None);
  }
  let parent = branch_lengths.get(&node.inbound()[0]).copied().flatten();
  let child = branch_lengths.get(&node.outbound()[0]).copied().flatten();
  (parent, child)
}

/// Record the two edges produced by [`split_edge`] into a branch-length value map, dropping the
/// original edge's entry.
pub fn record_split(branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>, info: &EdgeSplitInfo) {
  branch_lengths.remove(&info.old_edge_key);
  branch_lengths.insert(info.parent_side_edge_key, info.parent_side_length);
  branch_lengths.insert(info.child_side_edge_key, info.child_side_length);
}

/// Record the merged edge produced by [`remove_node_if_trivial`] into a branch-length value map,
/// dropping the two consumed edges' entries.
pub fn record_merge(branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>, info: &EdgeMergeInfo) {
  branch_lengths.remove(&info.parent_edge_key);
  branch_lengths.remove(&info.child_edge_key);
  branch_lengths.insert(info.merged_edge_key, info.merged_branch_length);
}
