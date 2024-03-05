use crate::graph::edge::{Edge, GraphEdge};
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, Node};
use itertools::Itertools;
use parking_lot::{RwLock, RwLockWriteGuard};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::sync::Arc;

pub enum GraphTraversalContinuation {
  Stop,
  Continue,
}

/// Policy trait, which defines how successors and predecessors of graph nodes and edges are resolved
/// during a particular type of breadth-first traversal
pub trait BfsTraversalPolicy<N, E>
where
  N: GraphNode,
  E: GraphEdge,
{
  /// Obtains successors of a node during traversal
  fn node_successors(graph: &Graph<N, E>, node: &Arc<RwLock<Node<N>>>) -> Vec<Arc<RwLock<Node<N>>>>;

  /// Obtains predecessors of a node during traversal
  fn node_predecessors(graph: &Graph<N, E>, node: &Arc<RwLock<Node<N>>>) -> Vec<Arc<RwLock<Node<N>>>>;

  /// Obtains predecessor node of an edge during traversal
  fn edge_predecessor(graph: &Graph<N, E>, edge: &Arc<RwLock<Edge<E>>>) -> Arc<RwLock<Node<N>>>;
}

/// Policy trait implementation, which defines how successors and predecessors of graph nodes and edges are resolved
/// during forward breadth-first traversal (from roots to leaves, along edge directions)
pub struct BfsTraversalPolicyForward;

impl<N, E> BfsTraversalPolicy<N, E> for BfsTraversalPolicyForward
where
  N: GraphNode,
  E: GraphEdge,
{
  /// Obtains successors of a node during forward traversal
  fn node_successors(graph: &Graph<N, E>, node: &Arc<RwLock<Node<N>>>) -> Vec<Arc<RwLock<Node<N>>>> {
    // During forward traversal, successors are the children (targets of outbound edges)
    graph
      .children_of(&node.read())
      .into_iter()
      .map(|(child, _)| child)
      .collect_vec()
  }

  /// Obtains successors of a node during forward traversal
  fn node_predecessors(graph: &Graph<N, E>, node: &Arc<RwLock<Node<N>>>) -> Vec<Arc<RwLock<Node<N>>>> {
    // During forward traversal, predecessors are the parents (sources of inbound edges)
    graph
      .parents_of(&node.read())
      .into_iter()
      .map(|(child, _)| child)
      .collect_vec()
  }

  /// Obtains predecessor node of an edge during forward traversal
  fn edge_predecessor(graph: &Graph<N, E>, edge: &Arc<RwLock<Edge<E>>>) -> Arc<RwLock<Node<N>>> {
    // During backward traversal, predecessor node of an edge is the source edge
    graph.get_node(edge.read().source()).unwrap()
  }
}

/// Policy trait implementation, which defines how successors and predecessors of graph nodes and edges are resolved
/// during backward breadth-first traversal (from leaves to roots, against edge directions)
pub struct BfsTraversalPolicyBackward;

impl<N, E> BfsTraversalPolicy<N, E> for BfsTraversalPolicyBackward
where
  N: GraphNode,
  E: GraphEdge,
{
  /// Obtains successors of a node during backward traversal
  fn node_successors(graph: &Graph<N, E>, node: &Arc<RwLock<Node<N>>>) -> Vec<Arc<RwLock<Node<N>>>> {
    // During backward traversal, successors are the parents (sources of inbound edges)
    graph
      .parents_of(&node.read())
      .into_iter()
      .map(|(child, _)| child)
      .collect_vec()
  }

  /// Obtains predecessors of a node during forward traversal
  fn node_predecessors(graph: &Graph<N, E>, node: &Arc<RwLock<Node<N>>>) -> Vec<Arc<RwLock<Node<N>>>> {
    // During backward traversal, predecessors are the children (targets of outbound edges)
    graph
      .children_of(&node.read())
      .into_iter()
      .map(|(child, _)| child)
      .collect_vec()
  }

  /// Obtains predecessor node of an edge during backward traversal
  fn edge_predecessor(graph: &Graph<N, E>, edge: &Arc<RwLock<Edge<E>>>) -> Arc<RwLock<Node<N>>> {
    // During backward traversal, predecessor node of an edge is the target edge
    graph.get_node(edge.read().target()).unwrap()
  }
}

/// Performs parallel forward breadth-first traversal (from roots to leaves, along edge directions)
pub fn directed_breadth_first_traversal_forward<N, E, F>(
  graph: &Graph<N, E>,
  sources: &[Arc<RwLock<Node<N>>>],
  explorer: F,
) where
  N: GraphNode,
  E: GraphEdge,
  F: Fn(&RwLockWriteGuard<Node<N>>) -> GraphTraversalContinuation + Sync + Send,
{
  directed_breadth_first_traversal::<N, E, F, BfsTraversalPolicyForward>(graph, sources, explorer);
}

/// Performs parallel backward breadth-first traversal (from leaves to roots, against edge directions)
pub fn directed_breadth_first_traversal_backward<N, E, F>(
  graph: &Graph<N, E>,
  sources: &[Arc<RwLock<Node<N>>>],
  explorer: F,
) where
  N: GraphNode,
  E: GraphEdge,
  F: Fn(&RwLockWriteGuard<Node<N>>) -> GraphTraversalContinuation + Sync + Send,
{
  directed_breadth_first_traversal::<N, E, F, BfsTraversalPolicyBackward>(graph, sources, explorer);
}

/// Implements parallel breadth-first traversal of a directed graph, given source nodes, exploration function and
/// a traversal polity type.
///
/// TraversalPolicy here is a generic type, that defines how to access predecessors and successors during a
/// concrete type of traversal.
pub fn directed_breadth_first_traversal<N, E, F, TraversalPolicy>(
  graph: &Graph<N, E>,
  sources: &[Arc<RwLock<Node<N>>>],
  explorer: F,
) where
  N: GraphNode,
  E: GraphEdge,
  F: Fn(&RwLockWriteGuard<Node<N>>) -> GraphTraversalContinuation + Sync + Send,
  TraversalPolicy: BfsTraversalPolicy<N, E>,
{
  // Walk the graph one "frontier" at a time. Frontier is a set of nodes of a "layer" in the graph, where each node
  // has its dependencies already resolved. Frontiers allow parallelism.
  let mut frontier = sources.to_vec();

  // We traverse the graph, gathering frontiers. The last frontier will be empty (nodes of the previous to last
  // frontier will have no unvisited successors), ending this loop.
  while !frontier.is_empty() {
    // Process each node in the current frontier concurrently
    let frontier_candidate_nodes: Vec<Arc<RwLock<Node<N>>>> = frontier
      .into_par_iter()
      .map(|node| {
        {
          let node = node.write();
          if !node.is_visited() {
            // The actual visit. Here we call the user-provided function.
            if let GraphTraversalContinuation::Stop = explorer(&node) {
                return vec![];
            }

            // We mark the node as visited so that it's not visited twice and telling following loop iterations
            // that its successors can potentially be processed now
            node.mark_as_visited();
          }
        }

        // Gather node's successors:
        //  - for forward traversal: children
        //  - for backwards traversal: parents
        TraversalPolicy::node_successors(graph, &node)
      })
        // For each node, we receive a list of its successors, so overall a list of lists. We flatten it here into a
        // flat list.
      .flatten()
      .collect();

    // NOTE: this is a barrier. The separation of the two loops is necessary for synchronization.

    // Decide which successors to add to the next frontier. We only add the successors which have ALL of its
    // predecessors already visited. This ensures exploration where each node is visited exactly once and only when all
    // of its dependencies are resolved.
    let next_frontier = frontier_candidate_nodes
      .into_par_iter()
      .filter(|candidate_node| {
        TraversalPolicy::node_predecessors(graph, candidate_node)
          .iter()
          .all(|asc| asc.read().is_visited())
      })
      .collect();

    // The newly gathered frontier becomes the new current frontier
    frontier = next_frontier;
  }
}
