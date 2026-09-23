use crate::edge::GraphEdgeKey;
use crate::graph::Graph;
use crate::node::{GraphNodeKey, Node};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use std::collections::{BTreeSet, VecDeque};
use treetime_utils::make_internal_report;

#[allow(
  clippy::multiple_inherent_impl,
  reason = "split across files by concern; see graph.rs for the primary impl"
)]
impl Graph {
  pub fn iter_depth_first_preorder_forward<F>(&self, mut explorer: F) -> Result<(), Report>
  where
    F: FnMut(GraphNodeForward) -> Result<(), Report>,
  {
    let root_key = self
      .get_exactly_one_root()
      .wrap_err("Graph must have exactly one root")?
      .key();
    let mut stack = vec![root_key];
    while let Some(node_key) = stack.pop() {
      let node = self
        .get_node(node_key)
        .ok_or_else(|| make_internal_report!("Node not found on the graph: {node_key}"))?;
      explorer(GraphNodeForward::new(self, node))?;
      for (child_key, _) in self.children_keys_of(node).rev() {
        stack.push(child_key);
      }
    }
    Ok(())
  }

  pub fn iter_depth_first_postorder_forward<F>(&self, mut explorer: F) -> Result<(), Report>
  where
    F: FnMut(GraphNodeBackward) -> Result<(), Report>,
  {
    let root_key = self
      .get_exactly_one_root()
      .wrap_err("Graph must have exactly one root")?
      .key();
    let mut stack = vec![root_key];
    let mut visited = BTreeSet::new();
    while let Some(node_key) = stack.pop() {
      if visited.insert(node_key) {
        stack.push(node_key);
        let node = self
          .get_node(node_key)
          .ok_or_else(|| make_internal_report!("Node not found on the graph: {node_key}"))?;
        for (child_key, _) in self.children_keys_of(node).rev() {
          if !visited.contains(&child_key) {
            stack.push(child_key);
          }
        }
      } else {
        let node = self
          .get_node(node_key)
          .ok_or_else(|| make_internal_report!("Node not found on the graph: {node_key}"))?;
        explorer(GraphNodeBackward::new(self, node))?;
      }
    }
    Ok(())
  }

  pub fn iter_breadth_first_forward<F>(&self, mut explorer: F) -> Result<(), Report>
  where
    F: FnMut(GraphNodeForward) -> Result<(), Report>,
  {
    let root_key = self
      .get_exactly_one_root()
      .wrap_err("Graph must have exactly one root")?
      .key();
    let mut queue = VecDeque::from([root_key]);

    while let Some(node_key) = queue.pop_front() {
      let node = self
        .get_node(node_key)
        .ok_or_else(|| make_internal_report!("Node not found on the graph: {node_key}"))?;
      explorer(GraphNodeForward::new(self, node))?;
      for (child_key, _) in self.children_keys_of(node) {
        queue.push_back(child_key);
      }
    }
    Ok(())
  }

  pub fn iter_breadth_first_backward<F>(&self, mut explorer: F) -> Result<(), Report>
  where
    F: FnMut(GraphNodeBackward) -> Result<(), Report>,
  {
    let root_key = self
      .get_exactly_one_root()
      .wrap_err("Graph must have exactly one root")?
      .key();
    let mut queue = VecDeque::from([root_key]);
    let mut order = Vec::new();
    while let Some(node_key) = queue.pop_front() {
      order.push(node_key);
      let node = self
        .get_node(node_key)
        .ok_or_else(|| make_internal_report!("Node not found on the graph: {node_key}"))?;
      let mut child_keys = self
        .children_keys_of(node)
        .map(|(child_key, _)| child_key)
        .collect_vec();
      child_keys.sort_unstable();
      queue.extend(child_keys);
    }
    for node_key in order.into_iter().rev() {
      let node = self
        .get_node(node_key)
        .ok_or_else(|| make_internal_report!("Node not found on the graph: {node_key}"))?;
      explorer(GraphNodeBackward::new(self, node))?;
    }
    Ok(())
  }
}

#[must_use]
#[derive(Debug)]
pub struct GraphNodeForward {
  pub is_root: bool,
  pub is_leaf: bool,
  pub key: GraphNodeKey,
  pub parent_keys: Vec<(GraphNodeKey, GraphEdgeKey)>,
  pub child_edge_keys: Vec<GraphEdgeKey>,
}

impl GraphNodeForward {
  fn new(graph: &Graph, node: &Node) -> Self {
    let is_leaf = node.is_leaf();
    let is_root = node.is_root();
    let key = node.key();

    let parent_keys = graph.parents_keys_of(node).collect_vec();
    let child_edge_keys = graph.children_keys_of(node).map(|(_, edge_key)| edge_key).collect_vec();

    Self {
      is_root,
      is_leaf,
      key,
      parent_keys,
      child_edge_keys,
    }
  }
}

#[must_use]
#[derive(Debug)]
pub struct GraphNodeBackward {
  pub is_root: bool,
  pub is_leaf: bool,
  pub key: GraphNodeKey,
  pub child_keys: Vec<(GraphNodeKey, GraphEdgeKey)>,
  pub parent_edge_keys: Vec<GraphEdgeKey>,
}

impl GraphNodeBackward {
  fn new(graph: &Graph, node: &Node) -> Self {
    let is_leaf = node.is_leaf();
    let is_root = node.is_root();
    let key = node.key();

    let child_keys = graph.children_keys_of(node).collect_vec();
    let parent_edge_keys = graph.parents_keys_of(node).map(|(_, edge_key)| edge_key).collect_vec();

    Self {
      is_root,
      is_leaf,
      key,
      child_keys,
      parent_edge_keys,
    }
  }
}
