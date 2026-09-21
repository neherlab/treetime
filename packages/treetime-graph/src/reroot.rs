use crate::edge::{GraphEdgeKey, invert_edge};
use crate::graph::Graph;
use crate::node::GraphNodeKey;
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EdgeSplitInfo {
  pub old_edge_key: GraphEdgeKey,
  pub new_node_key: GraphNodeKey,
  pub parent_side_edge_key: GraphEdgeKey,
  pub child_side_edge_key: GraphEdgeKey,
  pub parent_side_length: Option<f64>,
  pub child_side_length: Option<f64>,
  pub split_position: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EdgeMergeInfo {
  pub removed_node_key: GraphNodeKey,
  pub parent_edge_key: GraphEdgeKey,
  pub child_edge_key: GraphEdgeKey,
  pub merged_edge_key: GraphEdgeKey,
  pub merged_branch_length: Option<f64>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RerootResult {
  pub new_root_key: GraphNodeKey,
  pub edge_split: Option<EdgeSplitInfo>,
  pub edge_merge: Option<EdgeMergeInfo>,
  pub inverted_edge_keys: Vec<GraphEdgeKey>,
}

#[derive(Clone, Debug, Default)]
pub struct RerootChanges {
  pub edge_split: Option<EdgeSplitInfo>,
  pub edge_merge: Option<EdgeMergeInfo>,
  pub inverted_edge_keys: Vec<GraphEdgeKey>,
}

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
pub fn split_edge(
  graph: &mut Graph,
  edge_key: GraphEdgeKey,
  split_position: f64,
  branch_length: Option<f64>,
) -> Result<EdgeSplitInfo, Report> {
  let new_node_key = graph.add_node();

  let (source_key, target_key) = {
    let edge = graph.get_edge(edge_key).expect("Edge not found");
    (edge.source(), edge.target())
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

pub fn apply_reroot_topology(
  graph: &mut Graph,
  old_root_key: GraphNodeKey,
  new_root_key: GraphNodeKey,
) -> Result<Vec<GraphEdgeKey>, Report> {
  let paths = graph.path_from_node_to_node(new_root_key, old_root_key)?;

  let mut inverted_edge_keys = Vec::new();
  for (_, edge) in &paths {
    if let Some(edge_key) = edge {
      inverted_edge_keys.push(*edge_key);
      invert_edge(graph, *edge_key);
    }
  }

  graph.build()?;
  Ok(inverted_edge_keys)
}

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
pub fn remove_node_if_trivial(
  graph: &mut Graph,
  node_key: GraphNodeKey,
  parent_branch: Option<f64>,
  child_branch: Option<f64>,
) -> Result<Option<EdgeMergeInfo>, Report> {
  let (parent_edge_key, child_edge_key) = {
    let node = graph.get_node(node_key).expect("Node not found");
    if node.inbound().len() != 1 || node.outbound().len() != 1 {
      return Ok(None);
    }
    (node.inbound()[0], node.outbound()[0])
  };

  let parent_key = graph.get_edge(parent_edge_key).expect("Parent edge not found").source();
  let child_key = graph.get_edge(child_edge_key).expect("Child edge not found").target();

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

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
pub fn trivial_node_branch_lengths(
  graph: &Graph,
  node_key: GraphNodeKey,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> (Option<f64>, Option<f64>) {
  let node = graph.get_node(node_key).expect("Node not found");
  if node.inbound().len() != 1 || node.outbound().len() != 1 {
    return (None, None);
  }
  let parent = branch_lengths.get(&node.inbound()[0]).copied().flatten();
  let child = branch_lengths.get(&node.outbound()[0]).copied().flatten();
  (parent, child)
}

pub fn record_split(branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>, info: &EdgeSplitInfo) {
  branch_lengths.remove(&info.old_edge_key);
  branch_lengths.insert(info.parent_side_edge_key, info.parent_side_length);
  branch_lengths.insert(info.child_side_edge_key, info.child_side_length);
}

pub fn record_merge(branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>, info: &EdgeMergeInfo) {
  branch_lengths.remove(&info.parent_edge_key);
  branch_lengths.remove(&info.child_edge_key);
  branch_lengths.insert(info.merged_edge_key, info.merged_branch_length);
}
