use crate::graph::node::Node;
use itertools::Itertools;
use parking_lot::RwLock;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::sync::Arc;

pub fn directed_breadth_first_traversal_parallel<N, E, F>(sources: &[Arc<RwLock<Node<N, E>>>], explorer: F)
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
  F: Fn(&RwLock<Node<N, E>>) + Sync + Send,
{
  // Walk the graph one "frontier" at a time. Frontier is a set of nodes of a "level" in the graph, where each node
  // has its dependencies already resolved. Frontiers allow parallelism.
  let mut frontier = sources.to_vec();

  while !frontier.is_empty() {
    // Process each node in the current frontier concurrently
    frontier = frontier
      .into_par_iter()
      .map(|node| {

        // The actual visit. Here we call the user-provided function.
        explorer(&node);

        // We mark the node as visited so that it's not visited twice and telling following loop iterations
        // that its children can potentially be processed now
        {
          node.write().mark_as_visited();
        }

        // Here we iterate over node's children and decide which children to add to the next frontier:
        // We only queue the children which have ALL of its parents already visited.
        // This ensures "lazy", "exploration" traversal algorithm, where each node is visited exactly once
        // (as opposed to "greedy", "search-like" or "shortest path" algorithm where some of the nodes can be skipped
        // if they are not in the shortest path)
        let node = node.read();
        let children = node.outbound()
          .iter()
          .map(|edge| edge.target())
          .filter(|child| {
            child
              .read()
              .inbound()
              .iter()
              .filter_map(|edge| edge.upgrade().map(|edge| Arc::clone(&edge.source())))
              .all(|childs_parent| childs_parent.read().is_visited())
          })
          .collect_vec();

        children
      })
        // For each node's iteration, we receive a list of its children, so overall a list of lists.
        // We flatten anc collect it into a new frontier.
      .flatten()
      .collect();
  }
}
