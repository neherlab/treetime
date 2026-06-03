use crate::edge::{Edge, GraphEdge};
use crate::graph::{Graph, SafeNode, SafeNodeRefMut};
use crate::node::GraphNode;
use eyre::Report;
use itertools::Itertools;
use parking_lot::{Mutex, RwLock};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::sync::Arc;
use treetime_utils::sync::mutex::extract_parallel_error;

pub enum GraphTraversalContinuation {
  Stop,
  Continue,
}

/// Policy trait, which defines how successors and predecessors of graph nodes and edges are resolved
/// during a particular type of breadth-first traversal
pub trait BfsTraversalPolicy<N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  /// Obtains successors of a node during traversal
  fn node_successors(graph: &Graph<N, E, D>, node: &SafeNode<N>) -> Vec<SafeNode<N>>;

  /// Obtains predecessors of a node during traversal
  fn node_predecessors(graph: &Graph<N, E, D>, node: &SafeNode<N>) -> Vec<SafeNode<N>>;

  /// Obtains predecessor node of an edge during traversal
  fn edge_predecessor(graph: &Graph<N, E, D>, edge: &Arc<RwLock<Edge<E>>>) -> SafeNode<N>;
}

/// Policy trait implementation, which defines how successors and predecessors of graph nodes and edges are resolved
/// during forward breadth-first traversal (from roots to leaves, along edge directions)
pub struct BfsTraversalPolicyForward;

impl<N, E, D> BfsTraversalPolicy<N, E, D> for BfsTraversalPolicyForward
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  /// Obtains successors of a node during forward traversal
  fn node_successors(graph: &Graph<N, E, D>, node: &SafeNode<N>) -> Vec<SafeNode<N>> {
    // During forward traversal, successors are the children (targets of outbound edges)
    graph
      .children_of(&node.read())
      .into_iter()
      .map(|(child, _)| child)
      .collect_vec()
  }

  /// Obtains successors of a node during forward traversal
  fn node_predecessors(graph: &Graph<N, E, D>, node: &SafeNode<N>) -> Vec<SafeNode<N>> {
    // During forward traversal, predecessors are the parents (sources of inbound edges)
    graph
      .parents_of(&node.read())
      .into_iter()
      .map(|(child, _)| child)
      .collect_vec()
  }

  /// Obtains predecessor node of an edge during forward traversal
  fn edge_predecessor(graph: &Graph<N, E, D>, edge: &Arc<RwLock<Edge<E>>>) -> SafeNode<N> {
    // During backward traversal, predecessor node of an edge is the source edge
    graph.get_node(edge.read_arc().source()).unwrap()
  }
}

/// Policy trait implementation, which defines how successors and predecessors of graph nodes and edges are resolved
/// during backward breadth-first traversal (from leaves to roots, against edge directions)
pub struct BfsTraversalPolicyBackward;

impl<N, E, D> BfsTraversalPolicy<N, E, D> for BfsTraversalPolicyBackward
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  /// Obtains successors of a node during backward traversal
  fn node_successors(graph: &Graph<N, E, D>, node: &SafeNode<N>) -> Vec<SafeNode<N>> {
    // During backward traversal, successors are the parents (sources of inbound edges)
    graph
      .parents_of(&node.read())
      .into_iter()
      .map(|(child, _)| child)
      .collect_vec()
  }

  /// Obtains predecessors of a node during forward traversal
  fn node_predecessors(graph: &Graph<N, E, D>, node: &SafeNode<N>) -> Vec<SafeNode<N>> {
    // During backward traversal, predecessors are the children (targets of outbound edges)
    graph
      .children_of(&node.read())
      .into_iter()
      .map(|(child, _)| child)
      .collect_vec()
  }

  /// Obtains predecessor node of an edge during backward traversal
  fn edge_predecessor(graph: &Graph<N, E, D>, edge: &Arc<RwLock<Edge<E>>>) -> SafeNode<N> {
    // During backward traversal, predecessor node of an edge is the target edge
    graph.get_node(edge.read().target()).unwrap()
  }
}

/// Performs parallel forward breadth-first traversal (from roots to leaves, along edge directions)
pub fn directed_breadth_first_traversal_forward<N, E, D, F>(
  graph: &Graph<N, E, D>,
  sources: &[SafeNode<N>],
  explorer: F,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
  F: Fn(&SafeNodeRefMut<N>) -> Result<GraphTraversalContinuation, Report> + Sync + Send,
{
  directed_breadth_first_traversal::<N, E, D, F, BfsTraversalPolicyForward>(graph, sources, explorer)
}

/// Performs parallel backward breadth-first traversal (from leaves to roots, against edge directions)
pub fn directed_breadth_first_traversal_backward<N, E, D, F>(
  graph: &Graph<N, E, D>,
  sources: &[SafeNode<N>],
  explorer: F,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
  F: Fn(&SafeNodeRefMut<N>) -> Result<GraphTraversalContinuation, Report> + Sync + Send,
{
  directed_breadth_first_traversal::<N, E, D, F, BfsTraversalPolicyBackward>(graph, sources, explorer)
}

/// Implements parallel breadth-first traversal of a directed graph, given source nodes, exploration function and
/// a traversal polity type.
///
/// TraversalPolicy here is a generic type, that defines how to access predecessors and successors during a
/// concrete type of traversal.
pub fn directed_breadth_first_traversal<N, E, D, F, TraversalPolicy>(
  graph: &Graph<N, E, D>,
  sources: &[SafeNode<N>],
  explorer: F,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
  F: Fn(&SafeNodeRefMut<N>) -> Result<GraphTraversalContinuation, Report> + Sync + Send,
  TraversalPolicy: BfsTraversalPolicy<N, E, D>,
{
  let error: Arc<Mutex<Option<Report>>> = Arc::new(Mutex::new(None));

  let mut frontier = sources.to_vec();

  while !frontier.is_empty() {
    let frontier_candidate_nodes: Vec<SafeNode<N>> = frontier
      .into_par_iter()
      .map(|node| {
        {
          let node = node.write_arc();
          if !node.is_visited() {
            match explorer(&node) {
              Ok(GraphTraversalContinuation::Continue) => {},
              Ok(GraphTraversalContinuation::Stop) => return vec![],
              Err(e) => {
                let mut guard = error.lock();
                if guard.is_none() {
                  *guard = Some(e);
                }
                return vec![];
              },
            }
            node.mark_as_visited();
          }
        }

        TraversalPolicy::node_successors(graph, &node)
      })
      .flatten()
      .collect();

    let next_frontier = frontier_candidate_nodes
      .into_par_iter()
      .filter(|candidate_node| {
        TraversalPolicy::node_predecessors(graph, candidate_node)
          .iter()
          .all(|asc| asc.read().is_visited())
      })
      .collect();

    frontier = next_frontier;
  }

  extract_parallel_error(error)
}
