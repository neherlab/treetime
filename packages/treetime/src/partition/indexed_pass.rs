use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey};
use treetime_utils::make_internal_report;

pub struct IndexedPass<N, E> {
  slots: Vec<IndexedPassSlot<N, E>>,
  node_indices: Vec<Option<usize>>,
  edge_indices: Vec<Option<usize>>,
  frontiers: Vec<std::ops::Range<usize>>,
}

pub struct IndexedPassSlot<N, E> {
  pub key: GraphNodeKey,
  pub node: N,
  pub parent_key: Option<GraphNodeKey>,
  pub parent_edge: Option<(GraphEdgeKey, E)>,
}

impl<N, E> IndexedPass<N, E> {
  pub fn new<GN, GE>(
    graph: &Graph<GN, GE, ()>,
    mut nodes: BTreeMap<GraphNodeKey, N>,
    mut edges: BTreeMap<GraphEdgeKey, E>,
    mut missing_node: impl FnMut(GraphNodeKey) -> Result<N, Report>,
  ) -> Result<Self, Report>
  where
    GN: GraphNode,
    GE: GraphEdge,
  {
    let frontiers = graph.breadth_first_frontiers_backward()?;
    let mut node_indices = vec![None; graph.num_nodes()];
    let edge_capacity = graph
      .get_edges()
      .iter()
      .map(|edge| edge.read_arc().key().as_usize())
      .max()
      .map_or(0, |index| index + 1);
    let mut edge_indices = vec![None; edge_capacity];
    let mut slots = Vec::with_capacity(nodes.len());
    let mut ranges = Vec::with_capacity(frontiers.len());

    for frontier in frontiers {
      let start = slots.len();
      for key in frontier {
        let node = nodes.remove(&key).map_or_else(|| missing_node(key), Ok)?;
        let (parent_key, parent_edge) = if let Some(edge_key) = graph.parent_inbound_edge(key)? {
          let parent_key = graph.get_source_node_key(edge_key)?;
          let edge = edges
            .remove(&edge_key)
            .ok_or_else(|| make_internal_report!("Partition edge {edge_key} is missing while indexing a pass"))?;
          (Some(parent_key), Some((edge_key, edge)))
        } else {
          (None, None)
        };
        let index = slots.len();
        node_indices[key.as_usize()] = Some(index);
        if let Some((edge_key, _)) = &parent_edge {
          edge_indices[edge_key.as_usize()] = Some(index);
        }
        slots.push(IndexedPassSlot {
          key,
          node,
          parent_key,
          parent_edge,
        });
      }
      ranges.push(start..slots.len());
    }

    if !nodes.is_empty() || !edges.is_empty() {
      return Err(make_internal_report!(
        "Partition contains {} stale nodes and {} stale edges while indexing a pass",
        nodes.len(),
        edges.len()
      ));
    }

    Ok(Self {
      slots,
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
      usize,
      &[IndexedPassSlot<N, E>],
      &mut [IndexedPassSlot<N, E>],
    ) -> Result<(), Report>,
  ) -> Result<(), Report> {
    for range in self.frontiers.clone() {
      let (completed, remaining) = self.slots.split_at_mut(range.start);
      let (frontier, _) = remaining.split_at_mut(range.len());
      visit(&self.node_indices, &self.edge_indices, 0, completed, frontier)?;
    }
    Ok(())
  }

  pub fn try_for_each_forward_frontier(
    &mut self,
    mut visit: impl FnMut(
      &[Option<usize>],
      &[Option<usize>],
      usize,
      &mut [IndexedPassSlot<N, E>],
      &[IndexedPassSlot<N, E>],
    ) -> Result<(), Report>,
  ) -> Result<(), Report> {
    for range in self.frontiers.clone().into_iter().rev() {
      let (remaining, completed) = self.slots.split_at_mut(range.end);
      let (_, frontier) = remaining.split_at_mut(range.start);
      visit(&self.node_indices, &self.edge_indices, range.end, frontier, completed)?;
    }
    Ok(())
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
}
