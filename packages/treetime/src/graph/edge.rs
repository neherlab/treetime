use crate::graph::core::{CLOSED, OPEN};
use crate::graph::node::Node;
use parking_lot::Mutex;
use std::fmt::{Debug, Display, Formatter};
use std::hash::Hash;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Weak};

/// Edge representing a connection between two nodes. Relevant data can be
/// stored in the edge atomically. Edge's target and source node's are
/// weak references and can't outlive the nodes they represent.
#[derive(Debug)]
pub struct Edge<N, E>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  pub source: Weak<Mutex<Node<N, E>>>,
  pub target: Weak<Mutex<Node<N, E>>>,
  pub data: Mutex<E>,
  lock: AtomicBool,
}

impl<N, E> Edge<N, E>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  /// Creates a new edge.
  pub fn new(source: Weak<Mutex<Node<N, E>>>, target: Weak<Mutex<Node<N, E>>>, data: E) -> Edge<N, E> {
    Edge {
      source,
      target,
      data: Mutex::new(data),
      lock: AtomicBool::new(OPEN),
    }
  }

  /// Edge's source node.
  #[inline]
  pub fn source(&self) -> Arc<Mutex<Node<N, E>>> {
    self.source.upgrade().unwrap()
  }

  /// Edge's target node.
  #[inline]
  pub fn target(&self) -> Arc<Mutex<Node<N, E>>> {
    self.target.upgrade().unwrap()
  }

  /// Load data from the edge.
  #[inline]
  pub fn load(&self) -> E {
    self.data.lock().clone()
  }

  /// Store data into the edge.
  #[inline]
  pub fn store(&self, data: E) {
    let mut x = self.data.lock();
    *x = data;
  }

  #[inline]
  pub fn try_lock(&self) -> bool {
    self.lock.load(Ordering::Relaxed)
  }

  #[inline]
  pub fn close(&self) {
    self.lock.store(CLOSED, Ordering::Relaxed);
  }

  #[inline]
  pub fn open(&self) {
    self.lock.store(OPEN, Ordering::Relaxed);
  }
}

impl<N, E> Display for Edge<N, E>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  fn fmt(&self, fmt: &mut Formatter<'_>) -> std::fmt::Result {
    write!(fmt, "{} -> {}", self.source().lock().key(), self.target().lock().key())
  }
}
