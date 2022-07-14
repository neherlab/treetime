use crate::graph::edge::Edge;
use crate::graph::node::Node;
use parking_lot::{RwLock, RwLockWriteGuard};
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use std::borrow::Borrow;
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::sync::Arc;

/// Policy trait, which defines how descendants and ascendants of graph nodes and edges are resolved
/// during a particular type of breadth-first traversal
pub trait BfsTraversalPolicy<N, E>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  /// Obtains descendants of a node during traversal
  fn node_descendants(node: &Arc<RwLock<Node<N, E>>>) -> Vec<Arc<RwLock<Node<N, E>>>>;

  /// Obtains ascendants of a node during traversal
  fn node_ascendants(node: &Arc<RwLock<Node<N, E>>>) -> Vec<Arc<RwLock<Node<N, E>>>>;

  /// Obtains ascendant node of an edge during traversal
  fn edge_ascendant(edge: &Arc<Edge<N, E>>) -> Arc<RwLock<Node<N, E>>>;
}

/// Policy trait implementation, which defines how descendants and ascendants of graph nodes and edges are resolved
/// during forward breadth-first traversal (from roots to leaves, along edge directions)
pub struct BfsTraversalPolicyForward;

impl<N, E> BfsTraversalPolicy<N, E> for BfsTraversalPolicyForward
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  /// Obtains descendants of a node during forward traversal
  fn node_descendants(node: &Arc<RwLock<Node<N, E>>>) -> Vec<Arc<RwLock<Node<N, E>>>> {
    // During forward traversal, descendants are the children (targets of outbound edges)
    let node = node.read();
    let edges = node.outbound();
    edges.par_iter().map(|edge| edge.target()).collect()
  }

  /// Obtains descendants of a node during forward traversal
  fn node_ascendants(node: &Arc<RwLock<Node<N, E>>>) -> Vec<Arc<RwLock<Node<N, E>>>> {
    // During forward traversal, ascendants are the parents (sources of inbound edges)
    let node = node.read();
    let edges = node.inbound();
    edges
      .par_iter()
      .filter_map(|edge| edge.upgrade().map(|edge| edge.source()))
      .collect()
  }

  /// Obtains ascendant node of an edge during forward traversal
  fn edge_ascendant(edge: &Arc<Edge<N, E>>) -> Arc<RwLock<Node<N, E>>> {
    // During backward traversal, ascendant node of an edge is the source edge
    edge.source()
  }
}

/// Policy trait implementation, which defines how descendants and ascendants of graph nodes and edges are resolved
/// during backward breadth-first traversal (from leaves to roots, against edge directions)
pub struct BfsTraversalPolicyBackward;

impl<N, E> BfsTraversalPolicy<N, E> for BfsTraversalPolicyBackward
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  /// Obtains descendants of a node during backward traversal
  fn node_descendants(node: &Arc<RwLock<Node<N, E>>>) -> Vec<Arc<RwLock<Node<N, E>>>> {
    // During backward traversal, descendants are the parents (sources of inbound edges)
    let node = node.read();
    let edges = node.inbound();
    edges
      .par_iter()
      .filter_map(|edge| edge.upgrade().map(|edge| edge.source()))
      .collect()
  }

  /// Obtains ascendants of a node during forward backward
  fn node_ascendants(node: &Arc<RwLock<Node<N, E>>>) -> Vec<Arc<RwLock<Node<N, E>>>> {
    // During backward traversal, ascendants are the children (targets of outbound edges)
    let node = node.read();
    let edges = node.outbound();
    edges.par_iter().map(|edge| edge.target()).collect()
  }

  /// Obtains ascendant node of an edge during backward traversal
  fn edge_ascendant(edge: &Arc<Edge<N, E>>) -> Arc<RwLock<Node<N, E>>> {
    // During backward traversal, ascendant node of an edge is the target edge
    edge.target()
  }
}

/// Performs parallel forward breadth-first traversal (from roots to leaves, along edge directions)
pub fn directed_breadth_first_traversal_forward<N, E, F>(sources: &[Arc<RwLock<Node<N, E>>>], explorer: F)
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
  F: Fn(&RwLockWriteGuard<Node<N, E>>) + Sync + Send,
{
  directed_breadth_first_traversal::<N, E, F, BfsTraversalPolicyForward>(sources, explorer);
}

/// Performs parallel backward breadth-first traversal (from leaves to roots, against edge directions)
pub fn directed_breadth_first_traversal_backward<N, E, F>(sources: &[Arc<RwLock<Node<N, E>>>], explorer: F)
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
  F: Fn(&RwLockWriteGuard<Node<N, E>>) + Sync + Send,
{
  directed_breadth_first_traversal::<N, E, F, BfsTraversalPolicyBackward>(sources, explorer);
}

/// Implements parallel breadth-first traversal of a directed graph, given source nodes, exploration function and
/// a traversal polity type.
///
/// TraversalPolicy here is a generic type, that defines how to access ascendants and descendants during a
/// concrete type of traversal.
fn directed_breadth_first_traversal<N, E, F, TraversalPolicy>(sources: &[Arc<RwLock<Node<N, E>>>], explorer: F)
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
  F: Fn(&RwLockWriteGuard<Node<N, E>>) + Sync + Send,
  TraversalPolicy: BfsTraversalPolicy<N, E>,
{
  // Walk the graph one "frontier" at a time. Frontier is a set of nodes of a "layer" in the graph, where each node
  // has its dependencies already resolved. Frontiers allow parallelism.
  let mut frontier = sources.to_vec();

  // We traverse the graph, gathering frontiers. The last frontier will be empty (nodes of the previous to last
  // frontier will have no unvisited descendants), ending this loop.
  while !frontier.is_empty() {
    // Process each node in the current frontier concurrently
    let frontier_candidate_nodes: Vec<Arc<RwLock<Node<N, E>>>> = frontier
      .into_par_iter()
      .map(|node| {
        {
          let mut node = node.write();
          if !node.is_visited() {
            // The actual visit. Here we call the user-provided function.
            explorer(&node);

            // We mark the node as visited so that it's not visited twice and telling following loop iterations
            // that its descendants can potentially be processed now
            node.mark_as_visited();
          }
        }

        // Gather node's descendants:
        //  - for forward traversal: children
        //  - for backwards traversal: parents
        TraversalPolicy::node_descendants(&node)
      })
        // For each node, we receive a list of its descendants, so overall a list of lists. We flatten it here into a
        // flat list.
      .flatten()
      .collect();

    // NOTE: this is a barrier. The separation of the two loops is necessary for synchronization.

    // Decide which descendants to add to the next frontier. We only add the descendants which have ALL of its
    // ascendants already visited. This ensures exploration where each node is visited exactly once and only when all
    // of its dependencies are resolved.
    let next_frontier = frontier_candidate_nodes
      .into_par_iter()
      .filter(|candidate_node| {
        TraversalPolicy::node_ascendants(candidate_node)
          .iter()
          .all(|asc| asc.read().is_visited())
      })
      .collect();

    // The newly gathered frontier becomes the new current frontier
    frontier = next_frontier;
  }
}
