#[cfg(test)]
mod __tests__;

use crate::dependency_queue::{run_dependency_queue, validate_dependency_graph};
use crate::edge::{GraphEdge, GraphEdgeKey};
use crate::graph::Graph;
use crate::node::{GraphNode, GraphNodeKey};
use eyre::Report;
use parking_lot::Mutex;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::OnceLock;
use treetime_utils::make_internal_report;

pub fn with_indexed_graph_payloads<N, E, D, R>(
  graph: &Graph<N, E, D>,
  visit: impl FnOnce(&mut IndexedPass<N, E>) -> Result<R, Report>,
) -> Result<R, Report>
where
  N: GraphNode + Default,
  E: GraphEdge + Default,
  D: Send + Sync,
{
  let topology = IndexedPassTopology::new(graph)?;
  let mut nodes = graph
    .get_nodes()
    .iter()
    .map(|node| {
      let node = node.read_arc();
      let key = node.key();
      let data = std::mem::take(&mut *node.payload().write_arc());
      (key, data)
    })
    .collect::<BTreeMap<_, _>>();
  let mut edges = graph
    .get_edges()
    .iter()
    .map(|edge| {
      let edge = edge.read_arc();
      let key = edge.key();
      let data = std::mem::take(&mut *edge.payload().write_arc());
      (key, data)
    })
    .collect::<BTreeMap<_, _>>();

  let mut pass = IndexedPass::from_topology(&topology, &mut nodes, &mut edges, |_| {
    unreachable!("graph payload extraction includes every node")
  })?;
  let result = visit(&mut pass);
  let (mut nodes, mut edges) = pass.into_maps()?;

  for node in graph.get_nodes() {
    let node = node.read_arc();
    let key = node.key();
    *node.payload().write_arc() = nodes.remove(&key).expect("Indexed pass must restore every graph node");
  }
  for edge in graph.get_edges() {
    let edge = edge.read_arc();
    let key = edge.key();
    *edge.payload().write_arc() = edges.remove(&key).expect("Indexed pass must restore every graph edge");
  }

  result
}

pub struct IndexedPass<N, E> {
  slots: Vec<IndexedPassSlot<N, E>>,
  node_indices: Vec<Option<usize>>,
  edge_indices: Vec<Option<usize>>,
  parents: Vec<Option<usize>>,
  children: Vec<Vec<usize>>,
}

pub struct IndexedPassSlot<N, E> {
  pub key: GraphNodeKey,
  pub node: N,
  pub parent_key: Option<GraphNodeKey>,
  pub parent_edge: Option<(GraphEdgeKey, E)>,
}

pub struct IndexedPassDependencies<'a, N, E> {
  slots: &'a [OnceLock<IndexedPassSlot<N, E>>],
  node_indices: &'a [Option<usize>],
  edge_indices: &'a [Option<usize>],
}

impl<N, E> IndexedPass<N, E> {
  pub fn new<GN, GE>(
    graph: &Graph<GN, GE, impl Send + Sync>,
    nodes: &mut BTreeMap<GraphNodeKey, N>,
    edges: &mut BTreeMap<GraphEdgeKey, E>,
    missing_node: impl FnMut(GraphNodeKey) -> Result<N, Report>,
  ) -> Result<Self, Report>
  where
    GN: GraphNode,
    GE: GraphEdge,
  {
    let topology = IndexedPassTopology::new(graph)?;
    Self::from_topology(&topology, nodes, edges, missing_node)
  }

  pub fn try_for_each_backward(
    &mut self,
    visit: impl Fn(&IndexedPassDependencies<N, E>, &mut IndexedPassSlot<N, E>) -> Result<(), Report> + Sync + Send,
  ) -> Result<(), Report>
  where
    N: Send + Sync,
    E: Send + Sync,
  {
    let prerequisites = self.children.iter().map(Vec::len).collect::<Vec<_>>();
    let successors = self
      .parents
      .iter()
      .map(|parent| parent.iter().copied().collect::<Vec<_>>())
      .collect::<Vec<_>>();
    self.try_for_each_ready(&prerequisites, &successors, visit)
  }

  pub fn try_for_each_forward(
    &mut self,
    visit: impl Fn(&IndexedPassDependencies<N, E>, &mut IndexedPassSlot<N, E>) -> Result<(), Report> + Sync + Send,
  ) -> Result<(), Report>
  where
    N: Send + Sync,
    E: Send + Sync,
  {
    let prerequisites = self
      .parents
      .iter()
      .map(|parent| usize::from(parent.is_some()))
      .collect::<Vec<_>>();
    let successors = self.children.clone();
    self.try_for_each_ready(&prerequisites, &successors, visit)
  }

  pub fn into_maps(self) -> Result<(BTreeMap<GraphNodeKey, N>, BTreeMap<GraphEdgeKey, E>), Report> {
    let mut nodes = BTreeMap::new();
    let mut edges = BTreeMap::new();
    for slot in self.slots {
      if nodes.insert(slot.key, slot.node).is_some() {
        return Err(make_internal_report!(
          "Duplicate partition node {} while restoring a pass",
          slot.key
        ));
      }
      if let Some((edge_key, edge)) = slot.parent_edge
        && edges.insert(edge_key, edge).is_some()
      {
        return Err(make_internal_report!(
          "Duplicate partition edge {edge_key} while restoring a pass"
        ));
      }
    }
    Ok((nodes, edges))
  }

  fn from_topology(
    topology: &IndexedPassTopology,
    nodes: &mut BTreeMap<GraphNodeKey, N>,
    edges: &mut BTreeMap<GraphEdgeKey, E>,
    mut missing_node: impl FnMut(GraphNodeKey) -> Result<N, Report>,
  ) -> Result<Self, Report> {
    let graph_node_keys = topology.nodes.iter().map(|node| node.key).collect::<BTreeSet<_>>();
    let graph_edge_keys = topology
      .nodes
      .iter()
      .filter_map(|node| node.parent.map(|(_, edge_key)| edge_key))
      .collect::<BTreeSet<_>>();
    if nodes.keys().any(|key| !graph_node_keys.contains(key)) || edges.keys().any(|key| !graph_edge_keys.contains(key))
    {
      return Err(make_internal_report!(
        "Partition contains stale topology entries while indexing a pass"
      ));
    }
    let missing_keys = graph_node_keys
      .iter()
      .filter(|key| !nodes.contains_key(key))
      .copied()
      .collect::<Vec<_>>();
    for key in missing_keys {
      nodes.insert(key, missing_node(key)?);
    }
    if edges.len() != graph_edge_keys.len() || graph_edge_keys.iter().any(|key| !edges.contains_key(key)) {
      return Err(make_internal_report!(
        "Partition edge set does not match the graph while indexing a pass"
      ));
    }

    let mut nodes = std::mem::take(nodes);
    let mut edges = std::mem::take(edges);
    let node_capacity = topology
      .nodes
      .iter()
      .map(|node| node.key.as_usize())
      .max()
      .map_or(0, |index| index + 1);
    let edge_capacity = graph_edge_keys
      .iter()
      .map(|key| key.as_usize())
      .max()
      .map_or(0, |index| index + 1);
    let mut node_indices = vec![None; node_capacity];
    let mut edge_indices = vec![None; edge_capacity];
    let mut slots = Vec::with_capacity(topology.nodes.len());

    for node in &topology.nodes {
      let data = nodes.remove(&node.key).map_or_else(|| missing_node(node.key), Ok)?;
      let (parent_key, parent_edge) = if let Some((parent_key, edge_key)) = node.parent {
        let edge = edges.remove(&edge_key).expect("Validated partition edge must exist");
        (Some(parent_key), Some((edge_key, edge)))
      } else {
        (None, None)
      };
      let index = slots.len();
      node_indices[node.key.as_usize()] = Some(index);
      if let Some((edge_key, _)) = &parent_edge {
        edge_indices[edge_key.as_usize()] = Some(index);
      }
      slots.push(IndexedPassSlot {
        key: node.key,
        node: data,
        parent_key,
        parent_edge,
      });
    }

    let parents = topology
      .nodes
      .iter()
      .map(|node| {
        node
          .parent
          .map(|(key, _)| node_indices[key.as_usize()].expect("Indexed parent must exist"))
      })
      .collect::<Vec<_>>();
    let mut children = vec![Vec::new(); slots.len()];
    for (index, parent) in parents.iter().enumerate() {
      if let Some(parent) = parent {
        children[*parent].push(index);
      }
    }

    debug_assert!(nodes.is_empty() && edges.is_empty());

    Ok(Self {
      slots,
      node_indices,
      edge_indices,
      parents,
      children,
    })
  }

  fn try_for_each_ready(
    &mut self,
    prerequisites: &[usize],
    successors: &[Vec<usize>],
    visit: impl Fn(&IndexedPassDependencies<N, E>, &mut IndexedPassSlot<N, E>) -> Result<(), Report> + Sync + Send,
  ) -> Result<(), Report>
  where
    N: Send + Sync,
    E: Send + Sync,
  {
    let pending = std::mem::take(&mut self.slots)
      .into_iter()
      .map(|slot| Mutex::new(Some(slot)))
      .collect::<Vec<_>>();
    let completed = std::iter::repeat_with(OnceLock::new)
      .take(pending.len())
      .collect::<Vec<_>>();
    let dependencies = IndexedPassDependencies {
      slots: &completed,
      node_indices: &self.node_indices,
      edge_indices: &self.edge_indices,
    };

    let result = run_dependency_queue(prerequisites, successors, |index| {
      let mut slot = pending[index]
        .lock()
        .take()
        .expect("Dependency queue must schedule each indexed slot once");
      let result = visit(&dependencies, &mut slot);
      assert!(
        completed[index].set(slot).is_ok(),
        "Dependency queue must publish each indexed slot once"
      );
      result
    });

    self.slots = pending
      .into_iter()
      .zip(completed)
      .map(|(pending, completed)| {
        completed
          .into_inner()
          .or_else(|| pending.into_inner())
          .expect("Every indexed slot must remain available after traversal")
      })
      .collect();
    result
  }
}

impl<N, E> IndexedPassDependencies<'_, N, E> {
  pub fn slot(&self, key: GraphNodeKey) -> &IndexedPassSlot<N, E> {
    let index = self.node_indices[key.as_usize()].expect("Indexed dependency node must have a slot");
    self.slots[index]
      .get()
      .expect("Indexed dependency must complete before its successor")
  }

  pub fn node(&self, key: GraphNodeKey) -> &N {
    &self.slot(key).node
  }

  pub fn edge(&self, key: GraphEdgeKey) -> &E {
    let index = self.edge_indices[key.as_usize()].expect("Indexed dependency edge must have a slot");
    &self.slots[index]
      .get()
      .expect("Indexed dependency must complete before its successor")
      .parent_edge
      .as_ref()
      .expect("Indexed edge owner must have a parent edge")
      .1
  }
}

struct IndexedPassTopology {
  nodes: Vec<IndexedPassTopologyNode>,
}

impl IndexedPassTopology {
  fn new<N, E>(graph: &Graph<N, E, impl Send + Sync>) -> Result<Self, Report>
  where
    N: GraphNode,
    E: GraphEdge,
  {
    let nodes = graph
      .get_nodes()
      .iter()
      .map(|node| {
        let key = node.read_arc().key();
        let parent = if let Some(edge_key) = graph.parent_inbound_edge(key)? {
          Some((graph.get_source_node_key(edge_key)?, edge_key))
        } else {
          None
        };
        Ok(IndexedPassTopologyNode { key, parent })
      })
      .collect::<Result<Vec<_>, Report>>()?;
    let node_indices = nodes
      .iter()
      .enumerate()
      .map(|(index, node)| (node.key, index))
      .collect::<BTreeMap<_, _>>();
    let parents = nodes
      .iter()
      .map(|node| node.parent.map(|(key, _)| node_indices[&key]))
      .collect::<Vec<_>>();
    let mut children = vec![Vec::new(); nodes.len()];
    for (index, parent) in parents.iter().enumerate() {
      if let Some(parent) = parent {
        children[*parent].push(index);
      }
    }
    let prerequisites = children.iter().map(Vec::len).collect::<Vec<_>>();
    let successors = parents
      .iter()
      .map(|parent| parent.iter().copied().collect::<Vec<_>>())
      .collect::<Vec<_>>();
    validate_dependency_graph(&prerequisites, &successors)?;
    Ok(Self { nodes })
  }
}

struct IndexedPassTopologyNode {
  key: GraphNodeKey,
  parent: Option<(GraphNodeKey, GraphEdgeKey)>,
}
