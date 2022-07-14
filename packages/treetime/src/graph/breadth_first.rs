use crate::graph::node::Node;
use parking_lot::{RwLock, RwLockWriteGuard};
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::sync::Arc;

pub fn directed_breadth_first_traversal_parallel<N, E, F>(sources: &[Arc<RwLock<Node<N, E>>>], explorer: F)
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
  F: Fn(&RwLockWriteGuard<Node<N, E>>) + Sync + Send,
{
  // Walk the graph one "frontier" at a time. Frontier is a set of nodes of a "layer" in the graph, where each node
  // has its dependencies already resolved. Frontiers allow parallelism.
  let mut frontier = sources.to_vec();

  // We traverse the graph, gathering frontiers. The last frontier will be empty (nodes of the previous to last
  // frontier will have no unvisited children), ending this loop.
  while !frontier.is_empty() {
    // Process each node in the current frontier concurrently
    let frontier_candidates: Vec<Arc<RwLock<Node<N, E>>>> = frontier
      .into_par_iter()
      .map(|node| {
        {
          let mut node = node.write();
          if !node.is_visited() {
            // The actual visit. Here we call the user-provided function.
            explorer(&node);

            // We mark the node as visited so that it's not visited twice and telling following loop iterations
            // that its children can potentially be processed now
            node.mark_as_visited();
          }
        }

        // Gather node's children
        let node = node.read();
        let children: Vec<Arc<RwLock<Node<N, E>>>> = node.outbound()
          .par_iter()
          .map(|edge| edge.target())
          .collect();

        children
      })
        // For each node, we receive a list of its children, so overall a list of lists. We flatten it here into a
        // flat list.
      .flatten()
      .collect();

    // NOTE: this is a barrier. The separation of the two loops is necessary for synchronization.

    // Decide which children to add to the next frontier. We only add the children which have ALL of its parents
    // already visited. This ensures exploration where each node is visited exactly once and only when all of its
    // dependencies are resolved.
    let next_frontier = frontier_candidates
      .into_par_iter()
      .filter(|child| {
        child
          .read()
          .inbound()
          .iter()
          .filter_map(|edge| edge.upgrade().map(|edge| Arc::clone(&edge.source())))
          .all(|childs_parent| childs_parent.read().is_visited())
      })
      .collect();

    // The newly gathered frontier becomes the new current frontier
    frontier = next_frontier;
  }
}
