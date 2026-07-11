use eyre::Report;
use parking_lot::RwLock;
use std::collections::BTreeMap;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey};
use treetime_utils::make_internal_report;

pub struct IndexedPass<N, E> {
  slots: Vec<IndexedPassSlot<N>>,
  edges: Vec<RwLock<Option<(GraphEdgeKey, E)>>>,
  node_indices: Vec<Option<usize>>,
  edge_indices: Vec<Option<usize>>,
  frontiers: Vec<std::ops::Range<usize>>,
}

pub struct IndexedPassSlot<N> {
  pub key: GraphNodeKey,
  pub node: N,
  pub parent_key: Option<GraphNodeKey>,
}

impl<N, E> IndexedPass<N, E> {
  pub fn new<GN, GE>(
    graph: &Graph<GN, GE, ()>,
    nodes: &mut BTreeMap<GraphNodeKey, N>,
    edges: &mut BTreeMap<GraphEdgeKey, E>,
    mut missing_node: impl FnMut(GraphNodeKey) -> Result<N, Report>,
  ) -> Result<Self, Report>
  where
    GN: GraphNode,
    GE: GraphEdge,
  {
    let frontiers = graph.breadth_first_frontiers_backward()?;
    let topology = frontiers
      .iter()
      .flatten()
      .map(|key| {
        let parent = if let Some(edge_key) = graph.parent_inbound_edge(*key)? {
          Some((graph.get_source_node_key(edge_key)?, edge_key))
        } else {
          None
        };
        Ok((*key, parent))
      })
      .collect::<Result<BTreeMap<_, _>, Report>>()?;
    let graph_node_keys = topology.keys().copied().collect::<std::collections::BTreeSet<_>>();
    let graph_edge_keys = topology
      .values()
      .filter_map(|parent| parent.map(|(_, edge_key)| edge_key))
      .collect::<std::collections::BTreeSet<_>>();
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
    let mut node_indices = vec![None; graph.num_nodes()];
    let edge_capacity = graph
      .get_edges()
      .iter()
      .map(|edge| edge.read_arc().key().as_usize())
      .max()
      .map_or(0, |index| index + 1);
    let mut edge_indices = vec![None; edge_capacity];
    let mut slots = Vec::with_capacity(nodes.len());
    let mut edge_slots = Vec::with_capacity(nodes.len());
    let mut ranges = Vec::with_capacity(frontiers.len());

    for frontier in frontiers {
      let start = slots.len();
      for key in frontier {
        let node = nodes.remove(&key).map_or_else(|| missing_node(key), Ok)?;
        let parent = topology[&key];
        let (parent_key, parent_edge) = if let Some((parent_key, edge_key)) = parent {
          let edge = edges.remove(&edge_key).expect("Validated partition edge must exist");
          (Some(parent_key), Some((edge_key, edge)))
        } else {
          (None, None)
        };
        let index = slots.len();
        node_indices[key.as_usize()] = Some(index);
        if let Some((edge_key, _)) = &parent_edge {
          edge_indices[edge_key.as_usize()] = Some(index);
        }
        slots.push(IndexedPassSlot { key, node, parent_key });
        edge_slots.push(RwLock::new(parent_edge));
      }
      ranges.push(start..slots.len());
    }

    debug_assert!(nodes.is_empty() && edges.is_empty());

    Ok(Self {
      slots,
      edges: edge_slots,
      node_indices,
      edge_indices,
      frontiers: ranges,
    })
  }

  pub fn try_for_each_backward_frontier(
    &mut self,
    mut visit: impl FnMut(
      &[Option<usize>],
      &[Option<usize>],
      &[RwLock<Option<(GraphEdgeKey, E)>>],
      usize,
      &[IndexedPassSlot<N>],
      &mut [IndexedPassSlot<N>],
    ) -> Result<(), Report>,
  ) -> Result<(), Report> {
    for range in self.frontiers.clone() {
      let (completed, remaining) = self.slots.split_at_mut(range.start);
      let (frontier, _) = remaining.split_at_mut(range.len());
      visit(
        &self.node_indices,
        &self.edge_indices,
        &self.edges,
        0,
        completed,
        frontier,
      )?;
    }
    Ok(())
  }

  pub fn try_for_each_forward_frontier(
    &mut self,
    mut visit: impl FnMut(
      &[Option<usize>],
      &[Option<usize>],
      &[RwLock<Option<(GraphEdgeKey, E)>>],
      &[IndexedPassSlot<N>],
      usize,
      &mut [IndexedPassSlot<N>],
      &[IndexedPassSlot<N>],
    ) -> Result<(), Report>,
  ) -> Result<(), Report> {
    for range in self.frontiers.clone().into_iter().rev() {
      let (remaining, completed) = self.slots.split_at_mut(range.end);
      let (future, frontier) = remaining.split_at_mut(range.start);
      visit(
        &self.node_indices,
        &self.edge_indices,
        &self.edges,
        future,
        range.end,
        frontier,
        completed,
      )?;
    }
    Ok(())
  }

  pub fn into_maps(self) -> Result<(BTreeMap<GraphNodeKey, N>, BTreeMap<GraphEdgeKey, E>), Report> {
    let mut nodes = BTreeMap::new();
    let mut edges = BTreeMap::new();
    for (slot, edge_slot) in self.slots.into_iter().zip(self.edges) {
      if nodes.insert(slot.key, slot.node).is_some() {
        return Err(make_internal_report!(
          "Duplicate partition node {} while restoring a pass",
          slot.key
        ));
      }
      if let Some((edge_key, edge)) = edge_slot.into_inner()
        && edges.insert(edge_key, edge).is_some()
      {
        return Err(make_internal_report!(
          "Duplicate partition edge {edge_key} while restoring a pass"
        ));
      }
    }
    Ok((nodes, edges))
  }
}
