use crate::payload::clock_set::ClockSet;
use crate::payload::traits::ClockNode;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdge;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey, Outlier};
use treetime_graph::pass::{
  GraphMapOutputs, GraphPass, GraphPassBackwardContext, GraphPassForwardContext, GraphPassNodeOutput,
};

/// Per-node clock inference fields, held as a value keyed by [`GraphNodeKey`] instead of on the
/// graph node payload.
///
/// `time` is the observed or estimated date the regression reads for a leaf; `div` the cumulative
/// divergence; `bad_branch` and `is_outlier` the two exclusion flags; `clock_set` the accumulated
/// root-to-tip moment sums. The backward pass recomputes `clock_set` from scratch, so its seeded
/// value is never read; `time`, `bad_branch`, and `is_outlier` are durable inputs carried between
/// passes.
#[derive(Debug, Clone, Default)]
pub struct ClockNodeState {
  pub clock_set: ClockSet,
  pub div: f64,
  pub time: Option<f64>,
  pub bad_branch: bool,
  pub is_outlier: bool,
}

/// Per-edge clock messages, held as a value keyed by [`GraphEdgeKey`] instead of on the graph edge
/// payload. All three are recomputed by the regression passes and re-oriented on reroot.
#[derive(Debug, Clone, Default)]
pub struct ClockEdgeState {
  pub clock_to_parent: ClockSet,
  pub clock_to_child: ClockSet,
  /// The propagated `to_parent` message, kept to avoid recomputing the propagated message.
  pub clock_from_child: ClockSet,
}

/// The clock inference state for a whole tree, routed through the clock pipeline in place of the
/// `NodeClock`/`EdgeClock` (and `NodeTimetree`/`EdgeTimetree`) payload fields.
///
/// Keyed by stable node and edge ids, so the maps stay valid across a reroot, which adds a split
/// node, drops a trivial node, and re-orients the inverted path while leaving ids stable (gaps,
/// never renumbered).
#[derive(Debug, Clone, Default)]
pub struct ClockState {
  pub nodes: BTreeMap<GraphNodeKey, ClockNodeState>,
  pub edges: BTreeMap<GraphEdgeKey, ClockEdgeState>,
}

impl ClockState {
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
      .map(|node| (node.read_arc().key(), ClockNodeState::default()))
      .collect();
    let edges = graph
      .get_edges()
      .iter()
      .map(|edge| (edge.read_arc().key(), ClockEdgeState::default()))
      .collect();
    Self { nodes, edges }
  }

  /// Seed clock state from the graph payloads, reading the durable clock inputs each pass needs:
  /// `time` (via [`ClockNode::likely_time`]), `is_outlier`, `div`, and `clock_set`.
  ///
  /// Used at the timetree call sites, where these fields live on `NodeTimetree` between shared
  /// clock calls. Reading `time` here matches reading `likely_time()` live at the pass, because no
  /// timetree pass runs between a shared clock call's entry and its regression.
  pub fn seed_from_payloads<N, E, D>(graph: &Graph<N, E, D>) -> Self
  where
    N: GraphNode + ClockNode,
    E: GraphEdge,
    D: Send + Sync,
  {
    let nodes = graph
      .get_nodes()
      .iter()
      .map(|node| {
        let node = node.read_arc();
        let payload = node.payload().read_arc();
        let state = ClockNodeState {
          clock_set: payload.clock_set().clone(),
          div: ClockNode::div(&*payload),
          time: payload.likely_time(),
          bad_branch: false,
          is_outlier: payload.is_outlier(),
        };
        (node.key(), state)
      })
      .collect();
    let edges = graph
      .get_edges()
      .iter()
      .map(|edge| (edge.read_arc().key(), ClockEdgeState::default()))
      .collect();
    Self { nodes, edges }
  }

  /// Write the divergence and outlier flag of every node back into the graph payloads.
  ///
  /// The transitional repopulation the timetree pipeline needs after [`clock_filter`], whose own
  /// downstream (outlier bad-branch propagation, confidence intervals, tree writers) reads `div`
  /// and `is_outlier` off `NodeTimetree`.
  ///
  /// [`clock_filter`]: crate::clock::clock_filter::clock_filter_inplace
  pub fn write_div_is_outlier_to_payloads<N, E, D>(&self, graph: &Graph<N, E, D>)
  where
    N: GraphNode + ClockNode + Outlier,
    E: GraphEdge,
    D: Send + Sync,
  {
    for node_ref in graph.get_nodes() {
      let node = node_ref.read_arc();
      let state = self.node(node.key());
      let mut payload = node.payload().write_arc();
      ClockNode::set_div(&mut *payload, state.div);
      payload.set_is_outlier(state.is_outlier);
    }
  }

  #[must_use]
  pub fn node(&self, key: GraphNodeKey) -> &ClockNodeState {
    self
      .nodes
      .get(&key)
      .unwrap_or_else(|| panic!("Clock state is missing node {key}"))
  }

  #[must_use]
  pub fn node_mut(&mut self, key: GraphNodeKey) -> &mut ClockNodeState {
    self
      .nodes
      .get_mut(&key)
      .unwrap_or_else(|| panic!("Clock state is missing node {key}"))
  }

  #[must_use]
  pub fn edge(&self, key: GraphEdgeKey) -> &ClockEdgeState {
    self
      .edges
      .get(&key)
      .unwrap_or_else(|| panic!("Clock state is missing edge {key}"))
  }

  /// Run a value-returning backward pass over the clock state through the graph's dependency
  /// engine, replacing the per-node and per-edge maps with the visitor's outputs.
  ///
  /// The engine reads the node/edge inputs from the current state and reproduces the same
  /// thread-independent, deterministic child fold order as the payload-based passes.
  pub fn map_backward<GN, GE, D, F>(&mut self, graph: &Graph<GN, GE, D>, visit: F) -> Result<(), Report>
  where
    GN: GraphNode,
    GE: GraphEdge,
    D: Send + Sync,
    F: Fn(
        GraphPassBackwardContext<'_, ClockNodeState, ClockEdgeState, ClockNodeState, ClockEdgeState>,
      ) -> Result<GraphPassNodeOutput<ClockNodeState, ClockEdgeState>, Report>
      + Sync
      + Send,
  {
    let pass = GraphPass::new(graph, &mut self.nodes, &mut self.edges, |_| {
      Ok(ClockNodeState::default())
    })?;
    let GraphMapOutputs { nodes, edges } = pass.try_map_backward(visit)?;
    self.nodes = nodes;
    self.edges = edges;
    Ok(())
  }

  /// Run a value-returning forward pass over the clock state through the graph's dependency engine,
  /// replacing the per-node and per-edge maps with the visitor's outputs. Each node reads its single
  /// parent's already-published output.
  pub fn map_forward<GN, GE, D, F>(&mut self, graph: &Graph<GN, GE, D>, visit: F) -> Result<(), Report>
  where
    GN: GraphNode,
    GE: GraphEdge,
    D: Send + Sync,
    F: Fn(
        GraphPassForwardContext<'_, ClockNodeState, ClockEdgeState, ClockNodeState>,
      ) -> Result<GraphPassNodeOutput<ClockNodeState, ClockEdgeState>, Report>
      + Sync
      + Send,
  {
    let pass = GraphPass::new(graph, &mut self.nodes, &mut self.edges, |_| {
      Ok(ClockNodeState::default())
    })?;
    let GraphMapOutputs { nodes, edges } = pass.try_map_forward(visit)?;
    self.nodes = nodes;
    self.edges = edges;
    Ok(())
  }
}

impl ClockNodeState {
  /// The date the regression reads for this node, matching `NodeClock::likely_time`.
  #[must_use]
  pub fn likely_time(&self) -> Option<f64> {
    self.time
  }
}
