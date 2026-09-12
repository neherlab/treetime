use crate::edge::GraphEdgeKey;
use crate::graph::Graph;
use crate::node::{GraphNodeKey, Node};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use parking_lot::RwLock;
use std::collections::{BTreeSet, VecDeque};
use std::sync::Arc;
use traversal::Bft;

/// Represents graph node during forward traversal
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
  pub fn new(graph: &Graph, node: &Node) -> Self {
    let is_leaf = node.is_leaf();
    let is_root = node.is_root();
    let key = node.key();

    let parent_keys = graph
      .parents_of(node)
      .iter()
      .map(|(node, edge)| (node.read_arc().key(), edge.read_arc().key()))
      .collect_vec();

    let child_edge_keys = graph
      .children_of(node)
      .iter()
      .map(|(_, edge)| edge.read_arc().key())
      .collect_vec();

    Self {
      is_root,
      is_leaf,
      key,
      parent_keys,
      child_edge_keys,
    }
  }
}

/// Represents graph node during backwards traversal
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
  pub fn new(graph: &Graph, node: &Node) -> Self {
    let is_leaf = node.is_leaf();
    let is_root = node.is_root();
    let key = node.key();

    let child_keys = graph
      .children_of(node)
      .iter()
      .map(|(node, edge)| (node.read_arc().key(), edge.read_arc().key()))
      .collect_vec();

    let parent_edge_keys = graph
      .parents_of(node)
      .iter()
      .map(|(_, edge)| edge.read_arc().key())
      .collect_vec();

    Self {
      is_root,
      is_leaf,
      key,
      child_keys,
      parent_edge_keys,
    }
  }
}

/// Represents graph node during safe traversal
#[derive(Debug)]
pub struct GraphNodeSafe {
  pub is_root: bool,
  pub is_leaf: bool,
  pub key: GraphNodeKey,
}

impl GraphNodeSafe {
  pub fn from_node(node: &Arc<RwLock<Node>>) -> Self {
    let node = node.read();
    let is_leaf = node.is_leaf();
    let is_root = node.is_root();
    let key = node.key();
    Self { is_root, is_leaf, key }
  }
}

#[allow(
  clippy::multiple_inherent_impl,
  reason = "split across files by concern; see graph.rs for the primary impl"
)]
impl Graph {
  /// Serial depth-first preorder forward traversal (roots to leaves, parents before children).
  pub fn iter_depth_first_preorder_forward<F>(&self, mut explorer: F) -> Result<(), Report>
  where
    F: FnMut(GraphNodeForward) -> Result<(), Report>,
  {
    let root = self
      .get_exactly_one_root()
      .wrap_err("Graph must have exactly one root")?;
    let mut stack = Vec::from([(Arc::clone(&root), None)]);
    while let Some((current_node, _current_edge)) = stack.pop() {
      let current_node = current_node.read_arc();
      explorer(GraphNodeForward::new(self, &current_node))?;
      for (child, edge) in self.children_of(&current_node).into_iter().rev() {
        stack.push((child, Some(edge)));
      }
    }
    Ok(())
  }

  /// Serial depth-first postorder forward traversal (children before parents).
  pub fn iter_depth_first_postorder_forward<F>(&self, mut explorer: F) -> Result<(), Report>
  where
    F: FnMut(GraphNodeBackward) -> Result<(), Report>,
  {
    let root = self
      .get_exactly_one_root()
      .wrap_err("Graph must have exactly one root")?;
    let mut stack = Vec::new();
    let mut visited = BTreeSet::new();
    stack.push((Arc::clone(&root), None));
    while let Some((current_node, _)) = stack.pop() {
      let node_key = current_node.read_arc().key();
      if visited.insert(node_key) {
        stack.push((Arc::clone(&current_node), None));
        let children = self.children_of(&current_node.read_arc()).into_iter().rev();
        for (child, edge) in children {
          let child_key = child.read_arc().key();
          if !visited.contains(&child_key) {
            stack.push((child, Some(edge)));
          }
        }
      } else {
        explorer(GraphNodeBackward::new(self, &current_node.read_arc()))?;
      }
    }
    Ok(())
  }

  /// Serial breadth-first forward traversal (roots to leaves, along edge directions).
  ///
  /// Use this (rather than the parallel pass engine [`GraphPass::map_forward`](crate::pass::GraphPass::map_forward)) when the
  /// per-node work must capture mutable outer state, which a parallel callback cannot.
  pub fn iter_breadth_first_forward<F>(&self, mut explorer: F) -> Result<(), Report>
  where
    F: FnMut(GraphNodeForward) -> Result<(), Report>,
  {
    let root = self
      .get_exactly_one_root()
      .wrap_err("Graph must have exactly one root")?;
    let mut queue = VecDeque::new();
    queue.push_back(Arc::clone(&root));

    while let Some(current_node) = queue.pop_front() {
      explorer(GraphNodeForward::new(self, &current_node.read_arc()))?;
      let children = self.children_of(&current_node.read_arc());
      for (child, _) in children {
        queue.push_back(child);
      }
    }
    Ok(())
  }

  /// Serial breadth-first backward traversal (leaves to roots, against edge directions).
  pub fn iter_breadth_first_backward<F>(&self, mut explorer: F) -> Result<(), Report>
  where
    F: FnMut(GraphNodeBackward) -> Result<(), Report>,
  {
    let root = self
      .get_exactly_one_root()
      .wrap_err("Graph must have exactly one root")?;
    let nodes = Bft::new(&root, |node| self.iter_children_arc(node)).collect_vec();
    for (_, node) in nodes.into_iter().rev() {
      explorer(GraphNodeBackward::new(self, &node.write()))?;
    }
    Ok(())
  }

  fn iter_children_arc(&self, node: &Arc<RwLock<Node>>) -> impl Iterator<Item = &Arc<RwLock<Node>>> {
    let child_keys = self.child_keys_of(&node.read());
    self.nodes.iter().filter_map(move |node| {
      node
        .as_ref()
        .and_then(|node| child_keys.contains(&node.read_arc().key()).then_some(node))
    })
  }
}
