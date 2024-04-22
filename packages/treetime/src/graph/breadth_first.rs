use crate::graph::edge::{Edge, GraphEdge};
use crate::graph::graph::{Graph, SafeNode, SafeNodeRefMut};
use crate::graph::node::GraphNode;
use itertools::Itertools;
use parking_lot::RwLock;
use rayon::prelude::*;
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
  fn node_successors(graph: &Graph<N, E>, node: &SafeNode<N>) -> Vec<SafeNode<N>>;

  /// Obtains predecessors of a node during traversal
  fn node_predecessors(graph: &Graph<N, E>, node: &SafeNode<N>) -> Vec<SafeNode<N>>;

  /// Obtains predecessor node of an edge during traversal
  fn edge_predecessor(graph: &Graph<N, E>, edge: &Arc<RwLock<Edge<E>>>) -> SafeNode<N>;
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
  fn node_successors(graph: &Graph<N, E>, node: &SafeNode<N>) -> Vec<SafeNode<N>> {
    // During forward traversal, successors are the children (targets of outbound edges)
    graph
      .children_of(&node.read())
      .into_iter()
      .map(|(child, _)| child)
      .collect_vec()
  }

  /// Obtains successors of a node during forward traversal
  fn node_predecessors(graph: &Graph<N, E>, node: &SafeNode<N>) -> Vec<SafeNode<N>> {
    // During forward traversal, predecessors are the parents (sources of inbound edges)
    graph
      .parents_of(&node.read())
      .into_iter()
      .map(|(child, _)| child)
      .collect_vec()
  }

  /// Obtains predecessor node of an edge during forward traversal
  fn edge_predecessor(graph: &Graph<N, E>, edge: &Arc<RwLock<Edge<E>>>) -> SafeNode<N> {
    // During backward traversal, predecessor node of an edge is the source edge
    graph.get_node(edge.read_arc().source()).unwrap()
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
  fn node_successors(graph: &Graph<N, E>, node: &SafeNode<N>) -> Vec<SafeNode<N>> {
    // During backward traversal, successors are the parents (sources of inbound edges)
    graph
      .parents_of(&node.read())
      .into_iter()
      .map(|(child, _)| child)
      .collect_vec()
  }

  /// Obtains predecessors of a node during forward traversal
  fn node_predecessors(graph: &Graph<N, E>, node: &SafeNode<N>) -> Vec<SafeNode<N>> {
    // During backward traversal, predecessors are the children (targets of outbound edges)
    graph
      .children_of(&node.read())
      .into_iter()
      .map(|(child, _)| child)
      .collect_vec()
  }

  /// Obtains predecessor node of an edge during backward traversal
  fn edge_predecessor(graph: &Graph<N, E>, edge: &Arc<RwLock<Edge<E>>>) -> SafeNode<N> {
    // During backward traversal, predecessor node of an edge is the target edge
    graph.get_node(edge.read().target()).unwrap()
  }
}

/// Implements parallel breadth-first traversal of a directed graph, given source nodes, visitor function and
/// a traversal policy type.
///
/// Any given node is guaranteed to be visited after all of its predecessors and before all of its successors,
/// however the exact order of visits is not specified. Each node is visited exactly once.
///
/// TraversalPolicy here is a generic type that defines how to access predecessors and successors during a
/// concrete type of traversal.
pub fn directed_breadth_first_traversal<N, E, F, TraversalPolicy>(
  graph: &Graph<N, E>,
  sources: &[SafeNode<N>],
  visitor: F,
) where
  N: GraphNode,
  E: GraphEdge,
  F: Fn(&SafeNodeRefMut<N>) -> GraphTraversalContinuation + Sync + Send,
  TraversalPolicy: BfsTraversalPolicy<N, E>,
{
  directed_breadth_first_traversal_biphasic::<_, _, _, TraversalPolicy>(graph, sources, |(node, phase)| match phase {
    GraphTraversalPhase::Pre => visitor(&node.write_arc()),
    GraphTraversalPhase::Post => GraphTraversalContinuation::Continue,
  });
}

/// Phase during biphasic traversal
#[derive(Debug)]
pub enum GraphTraversalPhase {
  Pre,  // Node is being visited for the first time, before its successors have been visited
  Post, // Node is being visited for the second time, after its successors have been visited
}

/// Implements biphasic parallel breadth-first traversal of a directed graph, given source nodes, visitor function and
/// a traversal policy type.
///
/// Any given node is guaranteed to be visited after all of its predecessors, however the exact order of visits is not
/// specified. Each node is visited exactly twice:
///  - before its successors (the `phase` is `Pre`)
///  - after its successors (the `phase` is `Post`)
/// This allows to execute different code before and after traversal of node's successors.
///
/// TraversalPolicy here is a generic type that defines how to access predecessors and successors during a
/// concrete type of traversal.
pub fn directed_breadth_first_traversal_biphasic<N, E, F, TraversalPolicy>(
  graph: &Graph<N, E>,
  sources: &[SafeNode<N>],
  visitor: F,
) where
  N: GraphNode,
  E: GraphEdge,
  F: Fn((SafeNode<N>, GraphTraversalPhase)) -> GraphTraversalContinuation + Sync + Send,
  TraversalPolicy: BfsTraversalPolicy<N, E>,
{
  // Walk the graph one "frontier" at a time. Frontier is a set of nodes of a "layer" in the graph, where each node
  // has its dependencies already resolved. Frontiers allow parallelism.
  // TODO: try to further split up the frontier into generations, such that independent parts of the frontier don't
  // wait until the remainder of the frontier is processed.
  let mut frontier = sources
    .iter()
    .map(|node| (Arc::clone(node), GraphTraversalPhase::Pre))
    .collect_vec();

  // We traverse the graph, gathering frontiers. The last frontier will be empty (nodes of the previous to last
  // frontier will have no unvisited successors), ending this loop.
  while !frontier.is_empty() {
    // Process each node in the current frontier concurrently
    frontier = frontier
      .par_iter()
      .map(|(node, phase)| {
        match phase {
          GraphTraversalPhase::Pre => {
            if !node.write_arc().mark_as_visited_pre() {
              // The first ("pre") visit
              if let GraphTraversalContinuation::Stop = visitor((Arc::clone(node), GraphTraversalPhase::Pre)) {
                return vec![];
              }

              // Schedule node's successors for first ("pre") traversal. Successors are:
              //  - for forward traversal: children
              //  - for backward traversal: parents
              let mut successors = TraversalPolicy::node_successors(graph, node)
                .into_iter()
                .map(|succ| (Arc::clone(&succ), GraphTraversalPhase::Pre)).collect_vec();

              // Schedule yourself for second traversal ("post"), right after children.
              successors.push((Arc::clone(node), GraphTraversalPhase::Post));

              successors
            } else {
              vec![]
            }
          }
          GraphTraversalPhase::Post => {
            if !node.write_arc().mark_as_not_visited_post() {
              // The second ("post") visit
              visitor((Arc::clone(node), GraphTraversalPhase::Post));
            }
            vec![]
          }
        }
      })
        // For each node, we receive a list of its successors, so overall a list of lists. We flatten it here into a
        // flat list.
      .flatten()
      .collect();
  }
}
