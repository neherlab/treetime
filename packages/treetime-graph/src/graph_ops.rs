#[cfg(test)]
mod __tests__;

use crate::edge::{Edge, GraphEdgeKey};
use crate::graph::Graph;
use crate::node::{GraphNodeKey, Node};
use eyre::Report;
use treetime_utils::{make_error, make_internal_report, make_report};

#[allow(
  clippy::multiple_inherent_impl,
  reason = "split across files by concern; see graph.rs for the primary impl"
)]
impl Graph {
  pub fn add_node(&mut self) -> GraphNodeKey {
    let node_key = GraphNodeKey(self.nodes.len());
    self.nodes.push(Some(Node::new(node_key)));
    node_key
  }

  pub fn remove_node(&mut self, node_key: GraphNodeKey) -> Result<(Node, Vec<Edge>), Report> {
    let edges_to_remove: Vec<GraphEdgeKey> = self
      .edges
      .iter()
      .filter_map(|edge| {
        edge
          .as_ref()
          .and_then(|edge| (edge.source() == node_key || edge.target() == node_key).then_some(edge.key()))
      })
      .collect();

    let removed_edges = edges_to_remove
      .into_iter()
      .map(|edge_key| self.remove_edge(edge_key))
      .collect::<Result<Vec<_>, _>>()?;

    let removed_node = self
      .nodes
      .get_mut(node_key.as_usize())
      .and_then(Option::take)
      .ok_or_else(|| make_internal_report!("Attempted to remove non-existent node: {node_key}"))?;

    Ok((removed_node, removed_edges))
  }

  #[allow(
    clippy::expect_used,
    reason = "expect on a value an upstream invariant guarantees is present"
  )]
  pub fn add_edge(&mut self, source_key: GraphNodeKey, target_key: GraphNodeKey) -> Result<GraphEdgeKey, Report> {
    if source_key == target_key {
      return make_error!(
        "When adding a graph edge {source_key}->{target_key}: Attempted to connect node {source_key} to itself."
      );
    }

    let source = self.get_node(source_key).ok_or_else(|| {
      make_report!("When adding a graph edge {source_key}->{target_key}: Node {source_key} not found.")
    })?;

    if self.get_node(target_key).is_none() {
      return make_error!("When adding a graph edge {source_key}->{target_key}: Node {target_key} not found.");
    }

    let already_connected = source
      .outbound()
      .iter()
      .any(|edge| self.get_edge(*edge).is_some_and(|e| e.target() == target_key));

    if already_connected {
      return make_error!(
        "When adding a graph edge {source_key}->{target_key}: Nodes {source_key} and {target_key} are already connected."
      );
    }

    let edge_key = GraphEdgeKey(self.edges.len());
    self.edges.push(Some(Edge::new(edge_key, source_key, target_key)));
    self
      .get_node_mut(source_key)
      .expect("Edge source node must exist")
      .outbound_mut()
      .push(edge_key);
    self
      .get_node_mut(target_key)
      .expect("Edge target node must exist")
      .inbound_mut()
      .push(edge_key);

    Ok(edge_key)
  }

  #[allow(
    clippy::expect_used,
    reason = "expect on a value an upstream invariant guarantees is present"
  )]
  pub fn reparent_edge(&mut self, edge_key: GraphEdgeKey, new_source_key: GraphNodeKey) -> Result<(), Report> {
    let (old_source_key, target_key) = {
      let edge = self
        .get_edge(edge_key)
        .ok_or_else(|| make_internal_report!("When reparenting edge {edge_key}: edge not found"))?;
      (edge.source(), edge.target())
    };

    if old_source_key == new_source_key {
      return Ok(());
    }

    if new_source_key == target_key {
      return make_error!(
        "When reparenting edge {edge_key} to {new_source_key}: Attempted to connect node {new_source_key} to itself."
      );
    }

    {
      let new_source = self.get_node(new_source_key).ok_or_else(|| {
        make_report!("When reparenting edge {edge_key} to {new_source_key}: Node {new_source_key} not found.")
      })?;
      let already_connected = new_source
        .outbound()
        .iter()
        .any(|edge| self.get_edge(*edge).is_some_and(|e| e.target() == target_key));

      if already_connected {
        return make_error!(
          "When reparenting edge {edge_key} to {new_source_key}: Nodes {new_source_key} and {target_key} are already connected."
        );
      }
    }

    if let Some(old_source) = self.get_node_mut(old_source_key) {
      old_source.outbound_mut().retain(|&e| e != edge_key);
    }
    self
      .get_node_mut(new_source_key)
      .expect("New source node must exist")
      .outbound_mut()
      .push(edge_key);
    self
      .get_edge_mut(edge_key)
      .expect("Edge must exist")
      .set_source(new_source_key);

    Ok(())
  }

  pub fn remove_edge(&mut self, edge_key: GraphEdgeKey) -> Result<Edge, Report> {
    for node in self.nodes.iter_mut().flatten() {
      node.outbound_mut().retain(|&e| e != edge_key);
      node.inbound_mut().retain(|&e| e != edge_key);
    }

    self
      .edges
      .get_mut(edge_key.as_usize())
      .and_then(Option::take)
      .ok_or_else(|| make_internal_report!("Attempted to remove non-existent edge: {edge_key}"))
  }

  pub fn build(&mut self) -> Result<(), Report> {
    self.roots = self.get_nodes().filter(|node| node.is_root()).map(Node::key).collect();
    self.leaves = self.get_nodes().filter(|node| node.is_leaf()).map(Node::key).collect();
    Ok(())
  }

  #[allow(
    clippy::expect_used,
    reason = "expect on a value an upstream invariant guarantees is present"
  )]
  pub fn collapse_edge(&mut self, edge_key: GraphEdgeKey) -> Result<(Node, Edge, Vec<GraphEdgeKey>), Report> {
    let (source_key, target_key) = {
      let edge = self
        .get_edge(edge_key)
        .ok_or_else(|| make_internal_report!("Edge {edge_key} not found"))?;
      (edge.source(), edge.target())
    };

    let (target_inbound, target_outbound) = {
      let target_node = self
        .get_node(target_key)
        .ok_or_else(|| make_internal_report!("Target node {target_key} not found"))?;
      (target_node.inbound().to_vec(), target_node.outbound().to_vec())
    };

    for &inbound_edge_key in &target_inbound {
      if inbound_edge_key != edge_key && self.get_edge(inbound_edge_key).is_some() {
        self
          .get_edge_mut(inbound_edge_key)
          .expect("Inbound edge must exist")
          .set_target(source_key);
        if let Some(source_node) = self.get_node(source_key) {
          if !source_node.inbound().contains(&inbound_edge_key) {
            self
              .get_node_mut(source_key)
              .expect("Source node must exist")
              .inbound_mut()
              .push(inbound_edge_key);
          }
        }
      }
    }

    let mut new_edges = Vec::with_capacity(target_outbound.len());
    for &outbound_edge_key in &target_outbound {
      if outbound_edge_key != edge_key && self.get_edge(outbound_edge_key).is_some() {
        new_edges.push(outbound_edge_key);
        self
          .get_edge_mut(outbound_edge_key)
          .expect("Outbound edge must exist")
          .set_source(source_key);
        if let Some(source_node) = self.get_node(source_key) {
          if !source_node.outbound().contains(&outbound_edge_key) {
            self
              .get_node_mut(source_key)
              .expect("Source node must exist")
              .outbound_mut()
              .push(outbound_edge_key);
          }
        }
      }
    }

    if let Some(source_node) = self.get_node_mut(source_key) {
      source_node.outbound_mut().retain(|&e| e != edge_key);
    }

    let removed_edge = self.remove_edge(edge_key)?;
    let (removed_node, _removed_edges) = self.remove_node(target_key)?;

    Ok((removed_node, removed_edge, new_edges))
  }
}
