use crate::graph::breadth_first::{
  directed_breadth_first_traversal, BfsTraversalPolicy, BfsTraversalPolicyBackward, BfsTraversalPolicyForward,
  GraphTraversalContinuation,
};
use crate::graph::edge::{Edge, GraphEdge};
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, Node};
use crate::make_internal_report;
use eyre::Report;
use itertools::Itertools;
use parking_lot::RwLock;
use std::borrow::Borrow;
use std::collections::HashSet;
use std::sync::atomic::AtomicBool;
use std::sync::atomic::Ordering::Relaxed;
use std::sync::Arc;

/// Finds edges on all paths between two nodes
pub fn find_paths<N: GraphNode, E: GraphEdge>(
  graph: &Graph<N, E>,
  start: &Arc<RwLock<Node<N, E>>>,
  finish: &Arc<RwLock<Node<N, E>>>,
) -> Result<Vec<Arc<Edge<N, E>>>, Report> {
  let mut edge_indices = HashSet::<usize>::new();

  for (i, edge) in graph.get_edges().iter().enumerate() {
    // An edge (connecting nodes `source` and `target`) is on a path between nodes `start` and `end` iff:
    //  - there exists a forward path from `start` to `source`
    //  - there exists a backward path from `finish` to `target`
    let has_forward_path = exists_forward_path_between(start, &edge.source());
    graph.reset_nodes();
    let has_backward_path = exists_backward_path_between(finish, &edge.target());
    graph.reset_nodes();

    if has_forward_path && has_backward_path {
      edge_indices.insert(i);
    }
  }

  edge_indices
    .into_iter()
    .map(|i| {
      graph.get_edge(i).ok_or_else(|| {
        let start = start.read().key();
        let finish = finish.read().key();
        make_internal_report!(
          "When searching for paths between node #{start} and node #{finish} in the graph: requested edge with #{i} not found",
        )
      })
    })
    .collect()
}

/// Checks whether a path exists between two nodes of a graph in forward direction (leaves to roots)
pub fn exists_forward_path_between<N, E>(start: &Arc<RwLock<Node<N, E>>>, finish: &Arc<RwLock<Node<N, E>>>) -> bool
where
  N: GraphNode,
  E: GraphEdge,
{
  exists_path_between::<N, E, BfsTraversalPolicyForward>(start, finish)
}

/// Checks whether a path exists between two nodes of a graph in backward direction (leaves to roots)
pub fn exists_backward_path_between<N, E>(start: &Arc<RwLock<Node<N, E>>>, finish: &Arc<RwLock<Node<N, E>>>) -> bool
where
  N: GraphNode,
  E: GraphEdge,
{
  exists_path_between::<N, E, BfsTraversalPolicyBackward>(start, finish)
}

/// Checks whether a path exists between two nodes of a graph, given a traversal policy
pub fn exists_path_between<N, E, P>(start: &Arc<RwLock<Node<N, E>>>, finish: &Arc<RwLock<Node<N, E>>>) -> bool
where
  N: GraphNode,
  E: GraphEdge,
  P: BfsTraversalPolicy<N, E>,
{
  let path_exists = AtomicBool::new(false);
  let finish_node_key = finish.read().key();
  directed_breadth_first_traversal::<N, E, _, P>(&[Arc::clone(start)], |node| {
    if node.key() == finish_node_key {
      path_exists.store(true, Relaxed);
      return GraphTraversalContinuation::Stop;
    }
    GraphTraversalContinuation::Continue
  });
  path_exists.into_inner()
}
