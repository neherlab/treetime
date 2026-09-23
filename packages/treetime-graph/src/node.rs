use crate::edge::GraphEdgeKey;
use derive_more::Display;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use std::hash::Hash;

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
  #[inline]
  pub(crate) fn new(key: GraphNodeKey) -> Node {
    Self {
      key,
      outbound_edges: Vec::new(),
      inbound_edges: Vec::new(),
    }
  }

  #[inline]
  pub const fn key(&self) -> GraphNodeKey {
    self.key
  }

  #[inline]
  pub fn degree_out(&self) -> usize {
    self.outbound().len()
  }

  #[inline]
  pub(crate) fn degree_in(&self) -> usize {
    self.inbound().len()
  }

  #[inline]
  pub fn is_leaf(&self) -> bool {
    self.outbound().is_empty()
  }

  #[inline]
  pub fn is_root(&self) -> bool {
    self.inbound().is_empty()
  }

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
  pub(crate) fn has_at_most_one_child(&self) -> bool {
    self.degree_out() <= 1
  }

  #[inline]
  pub fn outbound(&self) -> &[GraphEdgeKey] {
    self.outbound_edges.as_slice()
  }

  #[inline]
  pub(crate) fn outbound_mut(&mut self) -> &mut Vec<GraphEdgeKey> {
    &mut self.outbound_edges
  }

  #[inline]
  pub fn inbound(&self) -> &[GraphEdgeKey] {
    self.inbound_edges.as_slice()
  }

  #[inline]
  pub(crate) fn inbound_mut(&mut self) -> &mut Vec<GraphEdgeKey> {
    &mut self.inbound_edges
  }
}

#[derive(Copy, Clone, Debug, Display, Eq, PartialEq, Ord, PartialOrd, Hash, Serialize, Deserialize)]
pub struct GraphNodeKey(pub usize);

impl GraphNodeKey {
  #[inline]
  pub const fn as_usize(self) -> usize {
    self.0
  }
}
