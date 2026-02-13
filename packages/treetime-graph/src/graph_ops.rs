use crate::edge::{Edge, GraphEdge, GraphEdgeKey};
use crate::graph::{Graph, SafeEdge};
use crate::node::{GraphNode, GraphNodeKey, Node};
use eyre::{Report, WrapErr};
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_utils::mutex::unwrap_arc_rwlock;
use treetime_utils::{make_error, make_internal_report};

impl<N, E, D> Graph<N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  pub fn add_node(&mut self, node_payload: N) -> GraphNodeKey {
    let node_key = GraphNodeKey(self.nodes.len());
    let node = Arc::new(RwLock::new(Node::new(node_key, node_payload)));
    self.nodes.push(Some(node));
    node_key
  }

  #[allow(clippy::needless_collect)]
  pub fn remove_node(&mut self, node_key: GraphNodeKey) -> Result<(Node<N>, Vec<Edge<E>>), Report> {
    let edges_to_remove: Vec<GraphEdgeKey> = self
      .edges
      .iter()
      .filter_map(|e| {
        e.as_ref().and_then(|e| {
          let e = e.read();
          (e.source() == node_key || e.target() == node_key).then_some(e.key())
        })
      })
      .collect();

    let removed_edges = edges_to_remove
      .into_iter()
      .map(|edge_key| -> Result<Edge<E>, Report> { self.remove_edge(edge_key) })
      .collect::<Result<Vec<_>, _>>()?;

    let removed_node = self
      .nodes
      .get_mut(node_key.as_usize())
      .and_then(|node| node.take().map(unwrap_arc_rwlock))
      .transpose()
      .wrap_err_with(|| format!("When removing node: {node_key}"))?
      .ok_or_else(|| make_internal_report!("Attempted to remove non-existent node: {node_key}"))?;

    Ok((removed_node, removed_edges))
  }

  /// Add a new edge to the graph.
  pub fn add_edge(
    &mut self,
    source_key: GraphNodeKey,
    target_key: GraphNodeKey,
    edge_payload: E,
  ) -> Result<GraphEdgeKey, Report> {
    if source_key == target_key {
      return make_error!(
        "When adding a graph edge {source_key}->{target_key}: Attempted to connect node {source_key} to itself."
      );
    }

    let source_lock = self
      .get_node(source_key)
      .ok_or_else(|| format!("When adding a graph edge {source_key}->{target_key}: Node {source_key} not found."))
      .unwrap();

    let target_lock = self
      .get_node(target_key)
      .ok_or_else(|| format!("When adding a graph edge {source_key}->{target_key}: Node {target_key} not found."))
      .unwrap();

    let edge_key = GraphEdgeKey(self.edges.len());
    let new_edge = Arc::new(RwLock::new(Edge::new(edge_key, source_key, target_key, edge_payload)));

    {
      let (source, target) = (source_lock.read(), target_lock.read());

      let already_connected = source
        .outbound()
        .iter()
        .any(|edge| self.get_edge(*edge).unwrap().read().target() == target.key());

      if already_connected {
        return make_error!(
          "When adding a graph edge {source_key}->{target_key}: Nodes {source_key} and {target_key} are already connected."
        );
      }

      self.edges.push(Some(Arc::clone(&new_edge)));
    }

    {
      let (mut source, mut target) = (source_lock.write(), target_lock.write());
      source.outbound_mut().push(edge_key);
      target.inbound_mut().push(edge_key);
    }

    Ok(edge_key)
  }

  pub fn remove_edge(&mut self, edge_key: GraphEdgeKey) -> Result<Edge<E>, Report> {
    // Remove the edge key from inbound/outbound lists of nodes
    self.nodes.iter_mut().for_each(|node| {
      if let Some(node) = node {
        let mut node_locked = node.write_arc();
        node_locked.outbound_mut().retain(|&e| e != edge_key);
        node_locked.inbound_mut().retain(|&e| e != edge_key);
      }
    });

    // Remove the edge itself
    self
      .edges
      .get_mut(edge_key.as_usize())
      .and_then(|edge_slot| edge_slot.take().map(unwrap_arc_rwlock))
      .transpose()
      .wrap_err_with(|| format!("When removing edge: {edge_key}"))?
      .ok_or_else(|| make_internal_report!("Attempted to remove non-existent edge: {edge_key}"))
  }

  pub fn build(&mut self) -> Result<(), Report> {
    self.roots = self
      .nodes
      .iter()
      .filter_map(|node_option| {
        node_option
          .as_ref()
          .and_then(|node| node.read().is_root().then(|| node.read().key()))
      })
      .collect();

    self.leaves = self
      .nodes
      .iter()
      .filter_map(|node_option| {
        node_option
          .as_ref()
          .and_then(|node| node.read().is_leaf().then(|| node.read().key()))
      })
      .collect();

    Ok(())
  }
  #[allow(clippy::type_complexity)]
  pub fn collapse_edge(&mut self, edge_key: GraphEdgeKey) -> Result<(Node<N>, Edge<E>, Vec<SafeEdge<E>>), Report>
  where
    N: Clone,
    E: Clone,
  {
    let (source_key, target_key) = {
      let edge = self
        .get_edge(edge_key)
        .ok_or_else(|| make_internal_report!("Edge {} not found", edge_key))?;
      let edge = edge.read_arc();
      (edge.source(), edge.target())
    };

    let (target_inbound, target_outbound) = {
      let target_node = self
        .get_node(target_key)
        .ok_or_else(|| make_internal_report!("Target node {} not found", target_key))?;
      let target_node = target_node.read_arc();
      (target_node.inbound().to_vec(), target_node.outbound().to_vec())
    };

    for &inbound_edge_key in &target_inbound {
      if inbound_edge_key != edge_key {
        if let Some(inbound_edge) = self.get_edge(inbound_edge_key) {
          inbound_edge.write_arc().set_target(source_key);
          if let Some(source_node) = self.get_node(source_key) {
            let mut source_node = source_node.write_arc();
            if !source_node.inbound().contains(&inbound_edge_key) {
              source_node.inbound_mut().push(inbound_edge_key);
            }
          }
        }
      }
    }

    let mut new_edges = Vec::with_capacity(target_outbound.len());
    for &outbound_edge_key in &target_outbound {
      if outbound_edge_key != edge_key {
        if let Some(outbound_edge) = self.get_edge(outbound_edge_key) {
          new_edges.push(Arc::clone(&outbound_edge));
          outbound_edge.write_arc().set_source(source_key);
          if let Some(source_node) = self.get_node(source_key) {
            let mut source_node = source_node.write_arc();
            if !source_node.outbound().contains(&outbound_edge_key) {
              source_node.outbound_mut().push(outbound_edge_key);
            }
          }
        }
      }
    }

    if let Some(source_node) = self.get_node(source_key) {
      let mut source_node = source_node.write_arc();
      source_node.outbound_mut().retain(|&e| e != edge_key);
    }

    let removed_edge = self.remove_edge(edge_key)?;
    let (removed_node, _removed_edges) = self.remove_node(target_key)?;

    Ok((removed_node, removed_edge, new_edges))
  }
}
