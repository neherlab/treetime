use crate::payload::traits::{TimetreeEdge, TimetreeNode};
use eyre::Report;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::{Distribution, NegLog};
use treetime_graph::edge::{GraphEdge, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey};
use treetime_graph::pass::{
  GraphMapOutputs, GraphPass, GraphPassBackwardContext, GraphPassForwardContext, GraphPassNodeOutput,
};

/// Per-node date-inference fields, held as a value keyed by [`GraphNodeKey`] instead of on the graph
/// node payload.
///
/// `time_distribution` is the node's posterior over its date, refined in place across the two date
/// passes; `time` the committed point estimate; `bad_branch` the exclusion flag a parent reads when
/// gathering child messages; `date_constraint` the fixed input date lifted back into the posterior on
/// every backward pass, carried read-only; `contradicted` a per-pass flag the forward pass raises when
/// the rest of the tree gives the given date no probability, folded into a diagnostic count and never
/// stored on a payload.
#[derive(Debug, Clone, Default)]
pub struct DateNodeState {
  pub time_distribution: Option<Arc<Distribution<NegLog>>>,
  pub time: Option<f64>,
  pub bad_branch: bool,
  pub date_constraint: Option<Arc<Distribution<NegLog>>>,
  pub contradicted: bool,
}

/// Per-edge date-inference fields, held as a value keyed by [`GraphEdgeKey`] instead of on the graph
/// edge payload.
///
/// `branch_length_distribution` is the branch's time-duration law the passes convolve across;
/// `msg_to_parent` the backward message the child sends up, divided back out as the cavity on the
/// forward pass; `time_length` the branch's committed duration (the Newick weight).
#[derive(Debug, Clone, Default)]
pub struct DateEdgeState {
  pub branch_length_distribution: Option<Arc<Distribution<NegLog>>>,
  pub msg_to_parent: Option<Arc<Distribution<NegLog>>>,
  pub time_length: Option<f64>,
}

/// The date-inference state for a whole tree, routed through the timetree date passes in place of the
/// `NodeTimetree`/`EdgeTimetree` payload fields.
///
/// Keyed by stable node and edge ids, so the maps stay valid across a reroot and polytomy resolution,
/// which add and drop nodes and edges while leaving ids stable (gaps, never renumbered).
#[derive(Debug, Clone, Default)]
pub struct TimetreeState {
  pub nodes: BTreeMap<GraphNodeKey, DateNodeState>,
  pub edges: BTreeMap<GraphEdgeKey, DateEdgeState>,
}

impl TimetreeState {
  /// Empty per-node/per-edge state for every node and edge of `graph`, all fields default.
  pub fn new<N, E, D>(graph: &Graph<N, E, D>) -> Self
  where
    N: GraphNode,
    E: GraphEdge,
    D: Send + Sync,
  {
    let nodes = graph
      .get_nodes()
      .iter()
      .map(|node| (node.read_arc().key(), DateNodeState::default()))
      .collect();
    let edges = graph
      .get_edges()
      .iter()
      .map(|edge| (edge.read_arc().key(), DateEdgeState::default()))
      .collect();
    Self { nodes, edges }
  }

  /// Seed date state from the graph payloads, reading the durable date fields the passes need.
  ///
  /// Used at the timetree date-pass call sites, where these fields live on `NodeTimetree`/`EdgeTimetree`
  /// between passes. `contradicted` has no payload counterpart and starts false.
  pub fn seed_from_payloads<N, E, D>(graph: &Graph<N, E, D>) -> Self
  where
    N: GraphNode + TimetreeNode,
    E: GraphEdge + TimetreeEdge,
    D: Send + Sync,
  {
    let nodes = graph
      .get_nodes()
      .iter()
      .map(|node| {
        let node = node.read_arc();
        let payload = node.payload().read_arc();
        let state = DateNodeState {
          time_distribution: payload.time_distribution().clone(),
          time: payload.time(),
          bad_branch: payload.bad_branch(),
          date_constraint: payload.date_constraint().clone(),
          contradicted: false,
        };
        (node.key(), state)
      })
      .collect();
    let edges = graph
      .get_edges()
      .iter()
      .map(|edge| {
        let edge = edge.read_arc();
        let payload = edge.payload().read_arc();
        let state = DateEdgeState {
          branch_length_distribution: payload.branch_length_distribution().clone(),
          msg_to_parent: payload.msg_to_parent().clone(),
          time_length: payload.time_length(),
        };
        (edge.key(), state)
      })
      .collect();
    Self { nodes, edges }
  }

  /// Write the date posterior, committed time, and backward messages back into the graph payloads.
  ///
  /// The transitional repopulation the timetree pipeline needs while the date passes run on the value
  /// but the refinement loop, coalescent statistics, confidence extraction, and tree writers still read
  /// `time`, `time_distribution`, and `msg_to_parent` off the payloads.
  pub fn write_to_payloads<N, E, D>(&self, graph: &Graph<N, E, D>)
  where
    N: GraphNode + TimetreeNode,
    E: GraphEdge + TimetreeEdge,
    D: Send + Sync,
  {
    for node_ref in graph.get_nodes() {
      let node = node_ref.read_arc();
      let state = self.node(node.key());
      let mut payload = node.payload().write_arc();
      payload.set_time_distribution(state.time_distribution.clone());
      payload.set_time(state.time);
    }
    for edge_ref in graph.get_edges() {
      let edge = edge_ref.read_arc();
      let state = self.edge(edge.key());
      let mut payload = edge.payload().write_arc();
      payload.set_msg_to_parent(state.msg_to_parent.clone());
    }
  }

  #[must_use]
  pub fn node(&self, key: GraphNodeKey) -> &DateNodeState {
    self
      .nodes
      .get(&key)
      .unwrap_or_else(|| panic!("Timetree state is missing node {key}"))
  }

  #[must_use]
  pub fn node_mut(&mut self, key: GraphNodeKey) -> &mut DateNodeState {
    self
      .nodes
      .get_mut(&key)
      .unwrap_or_else(|| panic!("Timetree state is missing node {key}"))
  }

  #[must_use]
  pub fn edge(&self, key: GraphEdgeKey) -> &DateEdgeState {
    self
      .edges
      .get(&key)
      .unwrap_or_else(|| panic!("Timetree state is missing edge {key}"))
  }

  /// Run a value-returning backward pass over the date state through the graph's dependency engine,
  /// replacing the per-node and per-edge maps with the visitor's outputs. Reproduces the same
  /// thread-independent, deterministic child fold order as the payload-based passes.
  pub fn map_backward<GN, GE, D, F>(&mut self, graph: &Graph<GN, GE, D>, visit: F) -> Result<(), Report>
  where
    GN: GraphNode,
    GE: GraphEdge,
    D: Send + Sync,
    F: Fn(
        GraphPassBackwardContext<'_, DateNodeState, DateEdgeState, DateNodeState, DateEdgeState>,
      ) -> Result<GraphPassNodeOutput<DateNodeState, DateEdgeState>, Report>
      + Sync
      + Send,
  {
    let pass = GraphPass::new(
      graph,
      &mut self.nodes,
      &mut self.edges,
      |_| Ok(DateNodeState::default()),
    )?;
    let GraphMapOutputs { nodes, edges } = pass.try_map_backward(visit)?;
    self.nodes = nodes;
    self.edges = edges;
    Ok(())
  }

  /// Run a value-returning forward pass over the date state through the graph's dependency engine,
  /// replacing the per-node and per-edge maps with the visitor's outputs. Each node reads its single
  /// parent's already-published output.
  pub fn map_forward<GN, GE, D, F>(&mut self, graph: &Graph<GN, GE, D>, visit: F) -> Result<(), Report>
  where
    GN: GraphNode,
    GE: GraphEdge,
    D: Send + Sync,
    F: Fn(
        GraphPassForwardContext<'_, DateNodeState, DateEdgeState, DateNodeState>,
      ) -> Result<GraphPassNodeOutput<DateNodeState, DateEdgeState>, Report>
      + Sync
      + Send,
  {
    let pass = GraphPass::new(
      graph,
      &mut self.nodes,
      &mut self.edges,
      |_| Ok(DateNodeState::default()),
    )?;
    let GraphMapOutputs { nodes, edges } = pass.try_map_forward(visit)?;
    self.nodes = nodes;
    self.edges = edges;
    Ok(())
  }
}
