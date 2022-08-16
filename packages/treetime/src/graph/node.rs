use crate::graph::edge::{Edge, GraphEdge};
use itertools::Itertools;
use parking_lot::{RwLock, RwLockReadGuard, RwLockWriteGuard};
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Weak};

type Outbound<N, E> = RwLock<Vec<Arc<Edge<N, E>>>>;
type Inbound<N, E> = RwLock<Vec<Weak<Edge<N, E>>>>;

pub trait Named {
  fn name(&self) -> &str;
}

pub trait GraphNode: Clone + Debug + Display + Sync + Send + Named {}

/// Represents a node in the graph. Data can be stored in and loaded from the
/// node in a thread safe manner.
#[derive(Debug)]
pub struct Node<N, E>
where
  N: GraphNode,
  E: GraphEdge,
{
  key: usize,
  data: Arc<RwLock<N>>,
  outbound: Outbound<N, E>,
  inbound: Inbound<N, E>,
  is_visited: AtomicBool,
}

impl<N, E> PartialEq<Self> for Node<N, E>
where
  E: GraphEdge,
  N: GraphNode,
{
  fn eq(&self, other: &Self) -> bool {
    self.key == other.key
  }
}

impl<N, E> Node<N, E>
where
  N: GraphNode,
  E: GraphEdge,
{
  /// Create a new node.
  #[inline]
  pub fn new(key: usize, data: N) -> Node<N, E> {
    Self {
      key,
      data: Arc::new(RwLock::new(data)),
      outbound: Outbound::new(Vec::new()),
      inbound: Inbound::new(Vec::new()),
      is_visited: AtomicBool::new(false),
    }
  }

  #[inline]
  pub fn payload(&self) -> Arc<RwLock<N>> {
    Arc::clone(&self.data)
  }

  /// Get node key.
  #[inline]
  pub const fn key(&self) -> usize {
    self.key
  }

  /// Get node degree ie. amount of outbound edges.
  #[inline]
  pub fn degree(&self) -> usize {
    self.outbound().len()
  }

  /// Check if node is a leaf node, i.e. has no outbound edges.
  #[inline]
  pub fn is_leaf(&self) -> bool {
    self.outbound().len() == 0
  }

  /// Check if node is a root node, i.e. has no inbound edges.
  #[inline]
  pub fn is_root(&self) -> bool {
    self.inbound().len() == 0
  }

  /// Get read access to outbound edges of the node.
  #[inline]
  pub fn outbound(&self) -> RwLockReadGuard<Vec<Arc<Edge<N, E>>>> {
    self.outbound.read()
  }

  /// Get read and write access to the outbound edges of the node. Will block other threads.
  #[inline]
  pub fn outbound_mut(&self) -> RwLockWriteGuard<Vec<Arc<Edge<N, E>>>> {
    self.outbound.write()
  }

  /// Get read access to inbound edges of the node.
  #[inline]
  pub fn inbound(&self) -> RwLockReadGuard<Vec<Weak<Edge<N, E>>>> {
    self.inbound.read()
  }

  /// Get read and write access to the outbound edges of the node. Will block other threads.
  #[inline]
  pub fn inbound_mut(&self) -> RwLockWriteGuard<Vec<Weak<Edge<N, E>>>> {
    self.inbound.write()
  }

  #[allow(clippy::type_complexity)]
  pub fn parents(&self) -> Vec<(Arc<RwLock<Node<N, E>>>, Arc<Edge<N, E>>)> {
    let parent_edges = self.inbound.read();
    parent_edges
      .iter()
      .cloned()
      .map(|parent_edge| {
        let parent_edge = parent_edge.upgrade().unwrap();
        let parent_node = parent_edge.source();
        (parent_node, parent_edge)
      })
      .collect_vec()
  }

  #[allow(clippy::type_complexity)]
  pub fn children(&self) -> Vec<(Arc<RwLock<Node<N, E>>>, Arc<Edge<N, E>>)> {
    let child_edges = self.outbound.read();
    child_edges
      .iter()
      .cloned()
      .map(|child_edge| {
        let child_node = child_edge.target();
        (child_node, child_edge)
      })
      .collect_vec()
  }

  #[inline]
  pub fn is_visited(&self) -> bool {
    self.is_visited.load(Ordering::Relaxed)
  }

  #[inline]
  pub fn mark_as_visited(&self) {
    self.is_visited.store(true, Ordering::Relaxed);
  }

  #[inline]
  pub fn mark_as_not_visited(&self) {
    self.is_visited.store(false, Ordering::Relaxed);
  }
}
