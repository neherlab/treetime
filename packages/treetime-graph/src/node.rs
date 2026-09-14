use crate::edge::GraphEdgeKey;
use derive_more::Display;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use std::hash::Hash;

#[derive(Copy, Clone, Debug, Display, Eq, PartialEq, Ord, PartialOrd, Hash, Serialize, Deserialize)]
pub struct GraphNodeKey(pub usize);

impl GraphNodeKey {
  #[inline]
  pub const fn as_usize(self) -> usize {
    self.0
  }
}

/// Internal representation of a node in a graph
#[derive(Debug, Serialize, Deserialize)]
pub struct Node {
  key: GraphNodeKey,
  outbound_edges: Vec<GraphEdgeKey>,
  inbound_edges: Vec<GraphEdgeKey>,
}

impl PartialEq<Self> for Node {
  fn eq(&self, other: &Self) -> bool {
    self.key == other.key
  }
}

impl Node {
  /// Create a new node.
  #[inline]
  pub fn new(key: GraphNodeKey) -> Node {
    Self {
      key,
      outbound_edges: Vec::new(),
      inbound_edges: Vec::new(),
    }
  }

  /// Get node key.
  #[inline]
  pub const fn key(&self) -> GraphNodeKey {
    self.key
  }

  /// Get node out-degree i.e. number of outbound edges.
  #[inline]
  pub fn degree_out(&self) -> usize {
    self.outbound().len()
  }

  /// Get node in-degree i.e. number of inbound edges.
  #[inline]
  pub fn degree_in(&self) -> usize {
    self.inbound().len()
  }

  /// Check if node is a leaf node, i.e. has no outbound edges.
  #[inline]
  pub fn is_leaf(&self) -> bool {
    self.outbound().is_empty()
  }

  /// Check if node is a root node, i.e. has no inbound edges.
  #[inline]
  pub fn is_root(&self) -> bool {
    self.inbound().is_empty()
  }

  /// Check if node is an internal node, i.e. has both inbound and outbound edges.
  #[inline]
  pub fn is_internal(&self) -> bool {
    !self.is_leaf() && !self.is_root()
  }

  #[inline]
  pub fn has_parents(&self) -> bool {
    self.degree_in() > 0
  }

  #[inline]
  pub fn has_one_parent(&self) -> bool {
    self.degree_in() == 1
  }

  #[inline]
  pub fn has_at_most_one_parent(&self) -> bool {
    self.degree_in() <= 1
  }

  #[inline]
  pub fn has_children(&self) -> bool {
    self.degree_out() > 0
  }

  #[inline]
  pub fn has_one_child(&self) -> bool {
    self.degree_out() == 1
  }

  #[inline]
  pub fn has_at_most_one_child(&self) -> bool {
    self.degree_out() <= 1
  }

  /// Get read access to outbound edges of the node.
  #[inline]
  pub fn outbound(&self) -> &[GraphEdgeKey] {
    self.outbound_edges.as_slice()
  }

  /// Get mutable access to the outbound edges of the node.
  #[inline]
  pub fn outbound_mut(&mut self) -> &mut Vec<GraphEdgeKey> {
    &mut self.outbound_edges
  }

  /// Get read access to inbound edges of the node.
  #[inline]
  pub fn inbound(&self) -> &[GraphEdgeKey] {
    self.inbound_edges.as_slice()
  }

  /// Get mutable access to the inbound edges of the node.
  #[inline]
  pub fn inbound_mut(&mut self) -> &mut Vec<GraphEdgeKey> {
    &mut self.inbound_edges
  }
}
