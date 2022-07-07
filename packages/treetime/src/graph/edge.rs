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
pub struct Edge<K, N, E>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  pub source: Weak<Node<K, N, E>>,
  pub target: Weak<Node<K, N, E>>,
  pub data: Mutex<E>,
  pub lock: AtomicBool,
}

impl<K, N, E> Edge<K, N, E>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  /// Creates a new edge.
  pub fn new(source: &Arc<Node<K, N, E>>, target: &Arc<Node<K, N, E>>, data: E) -> Edge<K, N, E> {
    Edge {
      source: Arc::downgrade(source),
      target: Arc::downgrade(target),
      data: Mutex::new(data),
      lock: AtomicBool::new(OPEN),
    }
  }

  /// Edge's source node.
  #[inline]
  pub fn source(&self) -> Arc<Node<K, N, E>> {
    self.source.upgrade().unwrap()
  }

  /// Edge's target node.
  #[inline]
  pub fn target(&self) -> Arc<Node<K, N, E>> {
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

impl<K, N, E> Display for Edge<K, N, E>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  fn fmt(&self, fmt: &mut Formatter<'_>) -> std::fmt::Result {
    write!(fmt, "{} -> {}", self.source().key(), self.target().key())
  }
}
