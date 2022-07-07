use crate::graph::core::{Continue, Traverse};
use crate::graph::edge::Edge;
use crate::graph::node::Node;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Weak};

pub fn directed_breadth_first_traversal_parallel<K, N, E, F>(
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
