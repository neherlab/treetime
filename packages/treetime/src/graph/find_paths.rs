use crate::graph::breadth_first::{
  directed_breadth_first_traversal, BfsTraversalPolicy, BfsTraversalPolicyBackward, BfsTraversalPolicyForward,
  GraphTraversalContinuation,
};
use crate::graph::edge::{Edge, GraphEdge, GraphEdgeKey};
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, Node};
use crate::make_internal_report;
use eyre::Report;
use parking_lot::RwLock;
use std::collections::HashSet;
use std::sync::atomic::AtomicBool;
use std::sync::atomic::Ordering::Relaxed;
use std::sync::Arc;

/// Finds edges on all paths between two nodes
pub fn find_paths<N, E, D>(
  graph: &Graph<N, E, D>,
  start: &Arc<RwLock<Node<N>>>,
  finish: &Arc<RwLock<Node<N>>>,
) -> Result<Vec<Arc<RwLock<Edge<E>>>>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  let mut edge_keys = HashSet::<GraphEdgeKey>::new();

  for edge in graph.get_edges() {
    // An edge (connecting nodes `source` and `target`) is on a path between nodes `start` and `end` iff:
    //  - there exists a forward path from `start` to `source`
    //  - there exists a backward path from `finish` to `target`
    let source = graph.get_node(edge.read().source()).unwrap();
    let has_forward_path = exists_forward_path_between(graph, start, &source);
    graph.reset_nodes();

    let target = graph.get_node(edge.read().target()).unwrap();
    let has_backward_path = exists_backward_path_between(graph, finish, &target);
    graph.reset_nodes();

    if has_forward_path && has_backward_path {
      edge_keys.insert(edge.read().key());
    }
  }

  edge_keys
    .into_iter()
    .map(|edge_key| {
      graph.get_edge(edge_key).ok_or_else(|| {
        let start = start.read().key();
        let finish = finish.read().key();
        make_internal_report!(
          "When searching for paths between node #{start} and node #{finish} in the graph: requested edge with #{edge_key} not found",
        )
      })
    })
    .collect()
}

/// Checks whether a path exists between two nodes of a graph in forward direction (leaves to roots)
pub fn exists_forward_path_between<N, E, D>(
  graph: &Graph<N, E, D>,
  start: &Arc<RwLock<Node<N>>>,
  finish: &Arc<RwLock<Node<N>>>,
) -> bool
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  exists_path_between::<N, E, D, BfsTraversalPolicyForward>(graph, start, finish)
}

/// Checks whether a path exists between two nodes of a graph in backward direction (leaves to roots)
pub fn exists_backward_path_between<N, E, D>(
  graph: &Graph<N, E, D>,
  start: &Arc<RwLock<Node<N>>>,
  finish: &Arc<RwLock<Node<N>>>,
) -> bool
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  exists_path_between::<N, E, D, BfsTraversalPolicyBackward>(graph, start, finish)
}

/// Checks whether a path exists between two nodes of a graph, given a traversal policy
pub fn exists_path_between<N, E, D, P>(
  graph: &Graph<N, E, D>,
  start: &Arc<RwLock<Node<N>>>,
  finish: &Arc<RwLock<Node<N>>>,
) -> bool
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
  P: BfsTraversalPolicy<N, E, D>,
{
  let path_exists = AtomicBool::new(false);
  let finish_node_key = finish.read().key();
  directed_breadth_first_traversal::<N, E, D, _, P>(graph, &[Arc::clone(start)], |node| {
    if node.key() == finish_node_key {
      path_exists.store(true, Relaxed);
      return GraphTraversalContinuation::Stop;
    }
    GraphTraversalContinuation::Continue
  });
  path_exists.into_inner()
}
