#[cfg(test)]
mod __tests__;

use crate::dependency_queue::{run_dependency_queue, validate_dependency_graph};
use crate::edge::GraphEdgeKey;
use crate::graph::Graph;
use crate::node::GraphNodeKey;
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::OnceLock;
use treetime_utils::make_internal_report;

pub struct GraphPass {
  nodes: Vec<GraphPassNode>,
  node_index: BTreeMap<GraphNodeKey, usize>,
  edge_keys: BTreeSet<GraphEdgeKey>,
  parents: Vec<Option<usize>>,
  children: Vec<Vec<usize>>,
}

impl GraphPass {
  pub fn new(graph: &Graph) -> Result<Self, Report> {
    let graph_nodes = graph.get_nodes().collect::<Vec<_>>();
    let mut nodes = Vec::with_capacity(graph_nodes.len());
    for &node in &graph_nodes {
      let key = node.key();
      let parent_edge = if let Some(edge_key) = graph.parent_inbound_edge(key)? {
        Some((graph.get_source_node_key(edge_key)?, edge_key))
      } else {
        None
      };
      nodes.push(GraphPassNode { key, parent_edge });
    }

    let node_index = nodes
      .iter()
      .enumerate()
      .map(|(index, node)| (node.key, index))
      .collect::<BTreeMap<_, _>>();

    let parents = nodes
      .iter()
      .map(|node| node.parent_edge.map(|(parent_key, _)| node_index[&parent_key]))
      .collect::<Vec<_>>();

    let mut children = vec![Vec::new(); nodes.len()];
    for &node in &graph_nodes {
      let parent_index = node_index[&node.key()];
      for (child, _edge) in graph.children_of(node) {
        let child_key = child.key();
        children[parent_index].push(node_index[&child_key]);
      }
    }

    let edge_keys = nodes
      .iter()
      .filter_map(|node| node.parent_edge.map(|(_, edge_key)| edge_key))
      .collect::<BTreeSet<_>>();

    let prerequisites = children.iter().map(Vec::len).collect::<Vec<_>>();
    let successors = parents
      .iter()
      .map(|parent| parent.iter().copied().collect::<Vec<_>>())
      .collect::<Vec<_>>();
    validate_dependency_graph(&prerequisites, &successors)?;

    Ok(Self {
      nodes,
      node_index,
      edge_keys,
      parents,
      children,
    })
  }

  #[allow(
    clippy::expect_used,
    reason = "expect on a value an upstream invariant guarantees is present"
  )]
  pub fn map_backward<N, E, NodeOut, EdgeOut>(
    &self,
    nodes: &BTreeMap<GraphNodeKey, N>,
    edges: &BTreeMap<GraphEdgeKey, E>,
    missing_node: impl FnMut(GraphNodeKey) -> Result<N, Report>,
    visit: impl Fn(
      GraphPassBackwardContext<'_, N, E, NodeOut, EdgeOut>,
    ) -> Result<GraphPassNodeOutput<NodeOut, EdgeOut>, Report>
    + Sync
    + Send,
  ) -> Result<GraphMapOutputs<NodeOut, EdgeOut>, Report>
  where
    N: Sync,
    E: Sync,
    NodeOut: Send + Sync,
    EdgeOut: Send + Sync,
  {
    let created = self.validate_and_create_missing(nodes, edges, missing_node)?;

    let prerequisites = self.children.iter().map(Vec::len).collect::<Vec<_>>();
    let successors = self
      .parents
      .iter()
      .map(|parent| parent.iter().copied().collect::<Vec<_>>())
      .collect::<Vec<_>>();

    let completed = std::iter::repeat_with(OnceLock::<GraphPassNodeOutput<NodeOut, EdgeOut>>::new)
      .take(self.nodes.len())
      .collect::<Vec<_>>();

    run_dependency_queue(&prerequisites, &successors, |index| {
      let node = &self.nodes[index];
      let input = self.resolve_node(nodes, &created, index);
      let parent_edge = node.parent_edge.map(|(_, edge_key)| (edge_key, &edges[&edge_key]));

      let children = self.children[index]
        .iter()
        .map(|&child_index| {
          let child = &self.nodes[child_index];
          let (_, edge_key) = child.parent_edge.expect("Backward child must have a parent edge");
          let output = completed[child_index]
            .get()
            .expect("Backward child must complete before its parent");
          GraphPassChildBackward {
            node_key: child.key,
            edge_key,
            node: &output.node,
            edge: output.parent_message.as_ref(),
          }
        })
        .collect::<Vec<_>>();

      let context = GraphPassBackwardContext {
        key: node.key,
        is_leaf: self.children[index].is_empty(),
        is_root: self.parents[index].is_none(),
        input,
        parent_edge,
        children: &children,
      };
      let output = visit(context)?;
      assert!(
        completed[index].set(output).is_ok(),
        "Dependency queue must publish each indexed slot once"
      );
      Ok(())
    })?;

    self.collect_map_outputs(completed)
  }

  #[allow(
    clippy::expect_used,
    reason = "expect on a value an upstream invariant guarantees is present"
  )]
  pub fn map_forward<N, E, NodeOut, EdgeOut>(
    &self,
    nodes: &BTreeMap<GraphNodeKey, N>,
    edges: &BTreeMap<GraphEdgeKey, E>,
    missing_node: impl FnMut(GraphNodeKey) -> Result<N, Report>,
    visit: impl Fn(GraphPassForwardContext<'_, N, E, NodeOut>) -> Result<GraphPassNodeOutput<NodeOut, EdgeOut>, Report>
    + Sync
    + Send,
  ) -> Result<GraphMapOutputs<NodeOut, EdgeOut>, Report>
  where
    N: Sync,
    E: Sync,
    NodeOut: Send + Sync,
    EdgeOut: Send + Sync,
  {
    let created = self.validate_and_create_missing(nodes, edges, missing_node)?;

    let prerequisites = self
      .parents
      .iter()
      .map(|parent| usize::from(parent.is_some()))
      .collect::<Vec<_>>();
    let successors = self.children.clone();

    let completed = std::iter::repeat_with(OnceLock::<GraphPassNodeOutput<NodeOut, EdgeOut>>::new)
      .take(self.nodes.len())
      .collect::<Vec<_>>();

    run_dependency_queue(&prerequisites, &successors, |index| {
      let node = &self.nodes[index];
      let input = self.resolve_node(nodes, &created, index);
      let parent_edge = node.parent_edge.map(|(_, edge_key)| (edge_key, &edges[&edge_key]));

      let parent = self.parents[index].map(|parent_index| {
        &completed[parent_index]
          .get()
          .expect("Forward parent must complete before its child")
          .node
      });

      let parent_key = self.parents[index].map(|parent_index| self.nodes[parent_index].key);
      let context = GraphPassForwardContext {
        key: node.key,
        is_leaf: self.children[index].is_empty(),
        is_root: self.parents[index].is_none(),
        input,
        parent_key,
        parent_edge,
        parent,
      };
      let output = visit(context)?;
      assert!(
        completed[index].set(output).is_ok(),
        "Dependency queue must publish each indexed slot once"
      );
      Ok(())
    })?;

    self.collect_map_outputs(completed)
  }

  fn validate_and_create_missing<N, E>(
    &self,
    nodes: &BTreeMap<GraphNodeKey, N>,
    edges: &BTreeMap<GraphEdgeKey, E>,
    mut missing_node: impl FnMut(GraphNodeKey) -> Result<N, Report>,
  ) -> Result<BTreeMap<GraphNodeKey, N>, Report> {
    if nodes.keys().any(|key| !self.node_index.contains_key(key))
      || edges.keys().any(|key| !self.edge_keys.contains(key))
    {
      return Err(make_internal_report!(
        "Partition contains stale topology entries while indexing a pass"
      ));
    }
    if edges.len() != self.edge_keys.len() || self.edge_keys.iter().any(|key| !edges.contains_key(key)) {
      return Err(make_internal_report!(
        "Partition edge set does not match the graph while indexing a pass"
      ));
    }
    let mut created = BTreeMap::new();
    for node in &self.nodes {
      if !nodes.contains_key(&node.key) {
        created.insert(node.key, missing_node(node.key)?);
      }
    }
    Ok(created)
  }

  #[allow(
    clippy::expect_used,
    reason = "expect on a value an upstream invariant guarantees is present"
  )]
  fn resolve_node<'a, N>(
    &self,
    nodes: &'a BTreeMap<GraphNodeKey, N>,
    created: &'a BTreeMap<GraphNodeKey, N>,
    index: usize,
  ) -> &'a N {
    let key = self.nodes[index].key;
    nodes
      .get(&key)
      .or_else(|| created.get(&key))
      .expect("Indexed node must have an input")
  }

  #[allow(
    clippy::expect_used,
    reason = "expect on a value an upstream invariant guarantees is present"
  )]
  fn collect_map_outputs<NodeOut, EdgeOut>(
    &self,
    completed: Vec<OnceLock<GraphPassNodeOutput<NodeOut, EdgeOut>>>,
  ) -> Result<GraphMapOutputs<NodeOut, EdgeOut>, Report> {
    let mut nodes = BTreeMap::new();
    let mut edges = BTreeMap::new();
    for (index, slot) in completed.into_iter().enumerate() {
      let output = slot
        .into_inner()
        .expect("Every indexed slot must publish an output after a successful map");
      let key = self.nodes[index].key;
      if nodes.insert(key, output.node).is_some() {
        return Err(make_internal_report!(
          "Duplicate node {key} while collecting a graph map"
        ));
      }
      if let Some(edge) = output.parent_message {
        let (_, edge_key) = self.nodes[index]
          .parent_edge
          .expect("Parent message must belong to a parent edge");
        if edges.insert(edge_key, edge).is_some() {
          return Err(make_internal_report!(
            "Duplicate edge {edge_key} while collecting a graph map"
          ));
        }
      }
    }
    Ok(GraphMapOutputs { nodes, edges })
  }
}

pub struct GraphPassChildBackward<'a, NodeOut, EdgeOut> {
  pub node_key: GraphNodeKey,
  pub edge_key: GraphEdgeKey,
  pub node: &'a NodeOut,
  pub edge: Option<&'a EdgeOut>,
}

pub struct GraphPassBackwardContext<'a, N, E, NodeOut, EdgeOut> {
  pub key: GraphNodeKey,
  pub is_leaf: bool,
  pub is_root: bool,
  pub input: &'a N,
  pub parent_edge: Option<(GraphEdgeKey, &'a E)>,
  pub children: &'a [GraphPassChildBackward<'a, NodeOut, EdgeOut>],
}

pub struct GraphPassForwardContext<'a, N, E, NodeOut> {
  pub key: GraphNodeKey,
  pub is_leaf: bool,
  pub is_root: bool,
  pub input: &'a N,
  pub parent_key: Option<GraphNodeKey>,
  pub parent_edge: Option<(GraphEdgeKey, &'a E)>,
  pub parent: Option<&'a NodeOut>,
}

pub struct GraphPassNodeOutput<NodeOut, EdgeOut> {
  pub node: NodeOut,
  pub parent_message: Option<EdgeOut>,
}

pub struct GraphMapOutputs<NodeOut, EdgeOut> {
  pub nodes: BTreeMap<GraphNodeKey, NodeOut>,
  pub edges: BTreeMap<GraphEdgeKey, EdgeOut>,
}

struct GraphPassNode {
  key: GraphNodeKey,
  parent_edge: Option<(GraphNodeKey, GraphEdgeKey)>,
}
