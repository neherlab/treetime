use crate::graph::node::{GraphNode, Node};
use parking_lot::RwLock;
use std::borrow::BorrowMut;
use std::fmt::{Debug, Display, Formatter};
use std::hash::Hash;
use std::mem::swap;
use std::sync::{Arc, Weak};

pub trait Weighted {
  fn weight(&self) -> f64;
}

pub trait GraphEdge: Clone + Debug + Display + Sync + Send + Weighted {}

/// Edge representing a connection between two nodes. Relevant data can be
/// stored in the edge atomically. Edge's target and source node's are
/// weak references and can't outlive the nodes they represent.
#[derive(Debug)]
pub struct Edge<N, E>
where
  N: GraphNode,
  E: GraphEdge,
{
  key: usize,
  source: Weak<RwLock<Node<N, E>>>,
  target: Weak<RwLock<Node<N, E>>>,
  data: Arc<RwLock<E>>,
}

impl<N, E> Edge<N, E>
where
  N: GraphNode,
  E: GraphEdge,
{
  /// Creates a new edge.
  pub fn new(key: usize, source: Weak<RwLock<Node<N, E>>>, target: Weak<RwLock<Node<N, E>>>, data: E) -> Edge<N, E> {
    Edge {
      key,
      source,
      target,
      data: Arc::new(RwLock::new(data)),
    }
  }

  #[inline]
  pub fn source(&self) -> Arc<RwLock<Node<N, E>>> {
    self.source.upgrade().unwrap()
  }

  #[inline]
  pub fn target(&self) -> Arc<RwLock<Node<N, E>>> {
    self.target.upgrade().unwrap()
  }

  #[inline]
  pub fn payload(&self) -> Arc<RwLock<E>> {
    Arc::clone(&self.data)
  }

  /// Invert direction of this edge.
  pub fn invert(&self) {
    let source = self.source();
    let target = self.target();
    swap(&mut source.write(), &mut target.write());
  }
}

impl<N, E> Display for Edge<N, E>
where
  N: GraphNode,
  E: GraphEdge,
{
  fn fmt(&self, fmt: &mut Formatter<'_>) -> std::fmt::Result {
    write!(fmt, "{} -> {}", self.source().read().key(), self.target().read().key())
  }
}
