use crate::breadth_first::{
  GraphTraversalContinuation, directed_breadth_first_traversal_backward, directed_breadth_first_traversal_forward,
};
use crate::edge::{GraphEdge, GraphEdgeKey};
use crate::graph::{Graph, NodeEdgePair, NodeEdgePayloadPair, SafeEdgePayloadRefMut, SafeNodePayloadRefMut};
use crate::node::{GraphNode, GraphNodeKey, Node};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use parking_lot::RwLock;
use std::collections::{BTreeSet, VecDeque};
use std::sync::Arc;
use traversal::Bft;
use treetime_utils::collections::container::{get_exactly_one, get_exactly_one_mut};

/// Represents graph node during forward traversal
#[must_use]
#[derive(Debug)]
pub struct GraphNodeForward<N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
{
  pub is_root: bool,
  pub is_leaf: bool,
  pub key: GraphNodeKey,
  pub payload: SafeNodePayloadRefMut<N>,
  pub parents: Vec<NodeEdgePayloadPair<N, E>>,
  pub parent_keys: Vec<(GraphNodeKey, GraphEdgeKey)>,
  pub child_edges: Vec<SafeEdgePayloadRefMut<E>>,
  pub child_edge_keys: Vec<GraphEdgeKey>,
  pub data: Arc<RwLock<D>>,
}

impl<N, E, D> GraphNodeForward<N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  pub fn new(graph: &Graph<N, E, D>, node: &Node<N>) -> Self {
    let is_leaf = node.is_leaf();
    let is_root = node.is_root();
    let key = node.key();

    let parent_keys = graph
      .parents_of(node)
      .iter()
      .map(|(node, edge)| (node.read_arc().key(), edge.read_arc().key()))
      .collect_vec();

    let child_edge_keys = graph
      .children_of(node)
      .iter()
      .map(|(_, edge)| edge.read_arc().key())
      .collect_vec();

    let payload = node.payload().write_arc();

    let parents = graph
      .parents_of(node)
      .iter()
      .map(|(node, edge)| (node.read().payload(), edge.read().payload()))
      .collect_vec();

    let child_edges = graph
      .children_of(node)
      .iter()
      .map(|(_, edge)| edge.write_arc().payload().write_arc())
      .collect_vec();

    let data = Arc::clone(&graph.data);

    Self {
      is_root,
      is_leaf,
      key,
      payload,
      parents,
      parent_keys,
      child_edges,
      child_edge_keys,
      data,
    }
  }

  pub fn get_exactly_one_parent(&self) -> Result<NodeEdgePayloadPair<N, E>, Report> {
    get_exactly_one(&self.parents)
      .cloned()
      .wrap_err("Nodes with multiple parents are not yet supported")
  }
}

/// Represents graph node during backwards traversal
#[must_use]
#[derive(Debug)]
pub struct GraphNodeBackward<N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  pub is_root: bool,
  pub is_leaf: bool,
  pub key: GraphNodeKey,
  pub payload: SafeNodePayloadRefMut<N>,
  pub children: Vec<NodeEdgePayloadPair<N, E>>,
  pub child_keys: Vec<(GraphNodeKey, GraphEdgeKey)>,
  pub parent_edges: Vec<SafeEdgePayloadRefMut<E>>,
  pub parent_edge_keys: Vec<GraphEdgeKey>,
  pub data: Arc<RwLock<D>>,
}

impl<N, E, D> GraphNodeBackward<N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  pub fn new(graph: &Graph<N, E, D>, node: &Node<N>) -> Self {
    let is_leaf = node.is_leaf();
    let is_root = node.is_root();
    let key = node.key();

    let child_keys = graph
      .children_of(node)
      .iter()
      .map(|(node, edge)| (node.read_arc().key(), edge.read_arc().key()))
      .collect_vec();

    let parent_edge_keys = graph
      .parents_of(node)
      .iter()
      .map(|(_, edge)| edge.read_arc().key())
      .collect_vec();

    let payload = node.payload().write_arc();

    let children = graph
      .children_of(node)
      .iter()
      .map(|(node, edge)| (node.read().payload(), edge.read().payload()))
      .collect_vec();

    let parent_edges = graph
      .parents_of(node)
      .iter()
      .map(|(_, edge)| edge.write_arc().payload().write_arc())
      .collect_vec();

    let data = Arc::clone(&graph.data);

    Self {
      is_root,
      is_leaf,
      key,
      payload,
      children,
      child_keys,
      parent_edges,
      parent_edge_keys,
      data,
    }
  }

  pub fn get_exactly_one_parent_edge(&mut self) -> Result<&mut SafeEdgePayloadRefMut<E>, Report> {
    get_exactly_one_mut(&mut self.parent_edges).wrap_err("Nodes with multiple parents are not yet supported")
  }
}

/// Represents graph node during safe traversal
#[derive(Debug)]
pub struct GraphNodeSafe<N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  pub is_root: bool,
  pub is_leaf: bool,
  pub key: GraphNodeKey,
  pub payload: Arc<RwLock<N>>,
  pub children: Vec<NodeEdgePair<N, E>>,
  pub parents: Vec<NodeEdgePair<N, E>>,
  pub data: Arc<RwLock<D>>,
}

impl<N, E, D> GraphNodeSafe<N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  pub fn from_node(graph: &Graph<N, E, D>, node: &Arc<RwLock<Node<N>>>) -> Self {
    let node = node.read();
    let is_leaf = node.is_leaf();
    let is_root = node.is_root();
    let key = node.key();
    let payload = node.payload();
    let parents = graph.parents_of(&node);
    let children = graph.children_of(&node);
    let data = Arc::clone(&graph.data);

    Self {
      is_root,
      is_leaf,
      key,
      payload,
      children,
      parents,
      data,
    }
  }
}

#[allow(
  clippy::multiple_inherent_impl,
  reason = "split across files by concern; see graph.rs for the primary impl"
)]
impl<N, E, D> Graph<N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  /// Return root-to-leaf frontiers whose nodes have only completed predecessors.
  pub fn breadth_first_frontiers_forward(&self) -> Result<Vec<Vec<GraphNodeKey>>, Report> {
    self.breadth_first_frontiers(self.roots.clone(), |node| node.outbound(), |node| node.inbound())
  }

  /// Return leaf-to-root frontiers whose nodes have only completed successors.
  pub fn breadth_first_frontiers_backward(&self) -> Result<Vec<Vec<GraphNodeKey>>, Report> {
    self.breadth_first_frontiers(self.leaves.clone(), |node| node.inbound(), |node| node.outbound())
  }

  /// Parallel breadth-first forward traversal (roots to leaves, along edge directions).
  ///
  /// The callback returns `Result<GraphTraversalContinuation, Report>`:
  /// - `Ok(Continue)`: mark node visited, expand successors
  /// - `Ok(Stop)`: stop expanding this node's branch (non-error early termination)
  /// - `Err(e)`: capture the first error, stop this branch, surface error to caller
  ///
  /// Other nodes already scheduled in the same frontier may still execute their callbacks
  /// before the traversal returns. Subsequent errors are discarded (first error wins).
  pub fn par_iter_breadth_first_forward<F>(&self, explorer: F) -> Result<(), Report>
  where
    F: Fn(GraphNodeForward<N, E, D>) -> Result<GraphTraversalContinuation, Report> + Sync + Send,
  {
    let roots = self.roots.iter().filter_map(|idx| self.get_node(*idx)).collect_vec();

    let result = directed_breadth_first_traversal_forward::<N, E, D, _>(self, roots.as_slice(), |node| {
      explorer(GraphNodeForward::new(self, node))
    });

    self.reset_nodes();
    result
  }

  /// Parallel breadth-first backward traversal (leaves to roots, against edge directions).
  ///
  /// See [`Self::par_iter_breadth_first_forward`] for callback semantics.
  pub fn par_iter_breadth_first_backward<F>(&self, explorer: F) -> Result<(), Report>
  where
    F: Fn(GraphNodeBackward<N, E, D>) -> Result<GraphTraversalContinuation, Report> + Sync + Send,
  {
    let leaves = self
      .leaves
      .iter()
      .filter_map(|idx| self.get_node(*idx))
      .rev()
      .collect_vec();

    let result = directed_breadth_first_traversal_backward::<N, E, D, _>(self, leaves.as_slice(), |node| {
      explorer(GraphNodeBackward::new(self, node))
    });

    self.reset_nodes();
    result
  }

  /// Serial depth-first preorder forward traversal (roots to leaves, parents before children).
  pub fn iter_depth_first_preorder_forward<F>(&self, mut explorer: F) -> Result<(), Report>
  where
    F: FnMut(GraphNodeForward<N, E, D>) -> Result<(), Report>,
  {
    let root = self
      .get_exactly_one_root()
      .wrap_err("Graph must have exactly one root")?;
    let mut stack = Vec::from([(Arc::clone(&root), None)]);
    while let Some((current_node, _current_edge)) = stack.pop() {
      let current_node = current_node.read_arc();
      explorer(GraphNodeForward::new(self, &current_node))?;
      for (child, edge) in self.children_of(&current_node).into_iter().rev() {
        stack.push((child, Some(edge)));
      }
    }
    Ok(())
  }

  /// Serial depth-first postorder forward traversal (children before parents).
  pub fn iter_depth_first_postorder_forward<F>(&self, mut explorer: F) -> Result<(), Report>
  where
    F: FnMut(GraphNodeBackward<N, E, D>) -> Result<(), Report>,
  {
    let root = self
      .get_exactly_one_root()
      .wrap_err("Graph must have exactly one root")?;
    let mut stack = Vec::new();
    let mut visited = BTreeSet::new();
    stack.push((Arc::clone(&root), None));
    while let Some((current_node, _)) = stack.pop() {
      let node_key = current_node.read_arc().key();
      if visited.insert(node_key) {
        stack.push((Arc::clone(&current_node), None));
        let children = self.children_of(&current_node.read_arc()).into_iter().rev();
        for (child, edge) in children {
          let child_key = child.read_arc().key();
          if !visited.contains(&child_key) {
            stack.push((child, Some(edge)));
          }
        }
      } else {
        explorer(GraphNodeBackward::new(self, &current_node.read_arc()))?;
      }
    }
    Ok(())
  }

  /// Serial breadth-first forward traversal (roots to leaves, along edge directions).
  ///
  /// Use this (rather than the parallel [`Self::par_iter_breadth_first_forward`]) when the
  /// per-node work must capture mutable outer state, which a parallel callback cannot.
  pub fn iter_breadth_first_forward<F>(&self, mut explorer: F) -> Result<(), Report>
  where
    F: FnMut(GraphNodeForward<N, E, D>) -> Result<(), Report>,
  {
    let root = self
      .get_exactly_one_root()
      .wrap_err("Graph must have exactly one root")?;
    let mut queue = VecDeque::new();
    queue.push_back(Arc::clone(&root));

    while let Some(current_node) = queue.pop_front() {
      explorer(GraphNodeForward::new(self, &current_node.read_arc()))?;
      let children = self.children_of(&current_node.read_arc());
      for (child, _) in children {
        queue.push_back(child);
      }
    }
    Ok(())
  }

  /// Serial breadth-first backward traversal (leaves to roots, against edge directions).
  pub fn iter_breadth_first_backward<F>(&self, mut explorer: F) -> Result<(), Report>
  where
    F: FnMut(GraphNodeBackward<N, E, D>) -> Result<(), Report>,
  {
    let root = self
      .get_exactly_one_root()
      .wrap_err("Graph must have exactly one root")?;
    let nodes = Bft::new(&root, |node| self.iter_children_arc(node)).collect_vec();
    for (_, node) in nodes.into_iter().rev() {
      explorer(GraphNodeBackward::new(self, &node.write()))?;
    }
    Ok(())
  }

  fn iter_children_arc(&self, node: &Arc<RwLock<Node<N>>>) -> impl Iterator<Item = &Arc<RwLock<Node<N>>>> {
    let child_keys = self.child_keys_of(&*node.read());
    self.nodes.iter().filter_map(move |node| {
      node
        .as_ref()
        .and_then(|node| child_keys.contains(&node.read_arc().key()).then_some(node))
    })
  }

  fn breadth_first_frontiers(
    &self,
    mut frontier: Vec<GraphNodeKey>,
    successors: impl Fn(&Node<N>) -> &[GraphEdgeKey],
    predecessors: impl Fn(&Node<N>) -> &[GraphEdgeKey],
  ) -> Result<Vec<Vec<GraphNodeKey>>, Report> {
    let mut visited = BTreeSet::new();
    let mut frontiers = Vec::new();

    while !frontier.is_empty() {
      frontier.sort();
      frontier.dedup();
      visited.extend(frontier.iter().copied());

      let candidates = frontier
        .iter()
        .map(|key| {
          let node = self.get_node(*key).ok_or_else(|| {
            treetime_utils::make_internal_report!("Node {key} not found while constructing frontiers")
          })?;
          successors(&node.read_arc())
            .iter()
            .map(|edge_key| {
              let edge = self.get_edge(*edge_key).ok_or_else(|| {
                treetime_utils::make_internal_report!("Edge {edge_key} not found while constructing frontiers")
              })?;
              let edge = edge.read_arc();
              Ok(if edge.source() == *key {
                edge.target()
              } else {
                edge.source()
              })
            })
            .collect::<Result<Vec<_>, Report>>()
        })
        .collect::<Result<Vec<_>, Report>>()?
        .into_iter()
        .flatten()
        .collect::<BTreeSet<_>>();

      frontiers.push(frontier);
      frontier = candidates
        .into_iter()
        .filter(|key| {
          self.get_node(*key).is_some_and(|node| {
            predecessors(&node.read_arc()).iter().all(|edge_key| {
              self.get_edge(*edge_key).is_some_and(|edge| {
                let edge = edge.read_arc();
                let predecessor = if edge.target() == *key {
                  edge.source()
                } else {
                  edge.target()
                };
                visited.contains(&predecessor)
              })
            })
          })
        })
        .collect();
    }

    if visited.len() != self.get_nodes().len() {
      return treetime_utils::make_internal_error!(
        "Cannot construct breadth-first frontiers: visited {} of {} nodes",
        visited.len(),
        self.get_nodes().len()
      );
    }

    Ok(frontiers)
  }
}
