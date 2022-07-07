use crate::graph::edge::Edge;
use crate::graph::node::Node;
use parking_lot::RwLock;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::{
  fmt::{Debug, Display},
  hash::Hash,
  sync::{
    atomic::{AtomicBool, Ordering},
    Arc, Weak,
  },
};

/// The Traverse enum is used when exploring edges in a graph. It's
/// states signify if an edge should be included in the search, skipped
/// or if the search should stop because a terminating condition has
/// been met such as finding a sink node.
pub enum Traverse {
  Include,
  Skip,
  Finish,
}

pub const OPEN: bool = false;
pub const CLOSED: bool = true;

pub enum Continue<T> {
  Yes(T),
  No(T),
}

/// Backtrack edges in an edge list starting from the last edge in the list.
/// Used for example to find the shortest path from the results of a breadth
/// first traversal.
pub fn backtrack_edges<K, N, E>(edges: &Vec<Weak<Edge<K, N, E>>>) -> Vec<Weak<Edge<K, N, E>>>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  let mut res = Vec::new();
  if edges.is_empty() {
    return res;
  }
  let w = &edges[edges.len() - 1];
  res.push(Weak::clone(w));
  let mut i = 0;
  for edge in edges.iter().rev() {
    let source = res[i].upgrade().unwrap().source();
    if edge.upgrade().unwrap().target().key == source.key {
      res.push(Weak::clone(edge));
      i += 1;
    }
  }
  res.reverse();
  res
}

// Opens all locks in
pub fn open_locks<K, N, E>(edges: &Vec<Weak<Edge<K, N, E>>>)
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  if edges.is_empty() {
    return;
  }
  edges[0].upgrade().unwrap().source().open();
  for weak in edges.iter() {
    let edge = weak.upgrade().unwrap();
    edge.open();
    edge.target().open();
  }
}

pub type Outbound<K, N, E> = RwLock<Vec<Arc<Edge<K, N, E>>>>;
pub type Inbound<K, N, E> = RwLock<Vec<Weak<Edge<K, N, E>>>>;

/// Connect two nodes if no previous connection exists.
pub fn connect<K, N, E>(source: &Arc<Node<K, N, E>>, target: &Arc<Node<K, N, E>>, data: E) -> bool
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  let connected = source
    .outbound()
    .iter()
    .any(|edge| edge.target.upgrade().unwrap().key == target.key);

  if !connected {
    let new_edge = Arc::new(Edge::new(source, target, data));
    source.outbound_mut().push(Arc::clone(&new_edge));
    target.inbound_mut().push(Arc::downgrade(&new_edge));
    true
  } else {
    false
  }
}

/// Disconnect two nodes from each other if they share an edge.
pub fn disconnect<K, N, E>(source: &Arc<Node<K, N, E>>, target: &Arc<Node<K, N, E>>) -> bool
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  let mut idx: (usize, usize) = (0, 0);
  let mut flag = false;
  for (i, edge) in source.outbound().iter().enumerate() {
    if edge.target().key == target.key {
      idx.0 = i;
      flag = true;
    }
  }
  for (i, edge) in target.inbound().iter().enumerate() {
    if edge.upgrade().unwrap().source().key == source.key {
      idx.1 = i;
    }
  }
  if flag {
    source.outbound_mut().remove(idx.0);
    source.inbound_mut().remove(idx.1);
  }
  flag
}

pub fn directed_breadth_first_traversal<K, N, E, F>(
  source: &Arc<Node<K, N, E>>,
  explorer: F,
) -> Option<Vec<Weak<Edge<K, N, E>>>>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
  F: Fn(&Arc<Edge<K, N, E>>) -> Traverse + Sync + Send + Copy,
{
  let mut frontiers: Vec<Weak<Edge<K, N, E>>>;
  let mut bounds: (usize, usize) = (0, 0);
  source.close();
  let initial = source.map_adjacent_dir(&explorer);
  match initial {
    Continue::No(segment) => {
      open_locks(&segment);
      return Some(segment);
    }
    Continue::Yes(segment) => {
      frontiers = segment;
    }
  }
  loop {
    bounds.1 = frontiers.len();
    if bounds.0 == bounds.1 {
      break;
    }
    let current_frontier = &frontiers[bounds.0..bounds.1];
    bounds.0 = bounds.1;
    let mut new_segments = Vec::new();
    for edge in current_frontier.iter() {
      let node = edge.upgrade().unwrap().target();
      let haystack = node.map_adjacent_dir(&explorer);
      match haystack {
        Continue::No(mut segment) => {
          new_segments.append(&mut segment);
          frontiers.append(&mut new_segments);
          open_locks(&frontiers);
          return Some(frontiers);
        }
        Continue::Yes(mut segment) => {
          new_segments.append(&mut segment);
        }
      }
    }
    frontiers.append(&mut new_segments);
  }
  open_locks(&frontiers);
  None
}

pub fn parallel_directed_breadth_first_traversal<K, N, E, F>(
  source: &Arc<Node<K, N, E>>,
  explorer: F,
) -> Option<Vec<Weak<Edge<K, N, E>>>>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
  F: Fn(&Arc<Edge<K, N, E>>) -> Traverse + Sync + Send + Copy,
{
  let mut frontiers: Vec<Weak<Edge<K, N, E>>>;
  let mut bounds: (usize, usize) = (0, 0);
  let terminate: Arc<AtomicBool> = Arc::new(AtomicBool::new(false));
  source.close();
  match source.map_adjacent_dir(&explorer) {
    Continue::No(segment) => {
      open_locks(&segment);
      return Some(segment);
    }
    Continue::Yes(segment) => {
      frontiers = segment;
    }
  }
  loop {
    bounds.1 = frontiers.len();
    if bounds.0 == bounds.1 {
      break;
    }
    let current_frontier = &frontiers[bounds.0..bounds.1];
    bounds.0 = bounds.1;
    let frontier_segments: Vec<_> = current_frontier
      .into_par_iter()
      .map(|edge| {
        if terminate.load(Ordering::Relaxed) {
          None
        } else {
          let node = edge.upgrade().unwrap().target();
          match node.map_adjacent_dir(&explorer) {
            Continue::No(segment) => {
              terminate.store(true, Ordering::Relaxed);
              Some(segment)
            }
            Continue::Yes(segment) => Some(segment),
          }
        }
      })
      .while_some()
      .collect();
    for mut segment in frontier_segments {
      frontiers.append(&mut segment);
    }
    if terminate.load(Ordering::Relaxed) {
      break;
    }
  }
  open_locks(&frontiers);
  terminate.load(Ordering::Relaxed).then(|| frontiers)
}
