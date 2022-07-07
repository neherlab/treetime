use crate::graph::edge::Edge;
use crate::graph::node::Node;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::sync::{Arc, Weak};

pub fn directed_breadth_first_traversal_parallel<N, E, F>(
  source: &Arc<Node<N, E>>,
  explorer: F,
) -> Vec<Weak<Edge<N, E>>>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
  F: Fn(&Arc<Edge<N, E>>) + Sync + Send,
{
  let mut bounds = (0, 0);

  source.close();
  let mut frontiers = source.map_adjacent_dir(&explorer);

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
        let node = edge.upgrade().unwrap().target();
        let segment = node.map_adjacent_dir(&explorer);
        Some(segment)
      })
      .while_some()
      .collect();

    for mut segment in frontier_segments {
      frontiers.append(&mut segment);
    }
  }

  open_locks(&frontiers);

  frontiers
}

pub fn open_locks<N, E>(edges: &Vec<Weak<Edge<N, E>>>)
where
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
