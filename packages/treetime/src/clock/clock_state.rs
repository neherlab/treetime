use crate::clock::clock_set::ClockSet;
use eyre::Report;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
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
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ClockNodeState {
  pub clock_set: ClockSet,
  pub div: f64,
  pub time: Option<f64>,
  pub bad_branch: bool,
  pub is_outlier: bool,
}

/// Per-edge clock messages, held as a value keyed by [`GraphEdgeKey`] instead of on the graph edge
/// payload. The three clock messages are recomputed by the regression passes and re-oriented on
/// reroot. `time_length` and `gamma` are the branch's solver-updated duration and relaxed-clock rate
/// multiplier the re-estimation reads to convert time back to divergence; they are seeded from the
/// date state in the refinement loop and stay at their defaults elsewhere, where the regression reads
/// input branch lengths instead.
#[derive(Debug, Clone, SmartDefault, PartialEq)]
pub struct ClockEdgeState {
  pub clock_to_parent: ClockSet,
  pub clock_to_child: ClockSet,
  /// The propagated `to_parent` message, kept to avoid recomputing the propagated message.
  pub clock_from_child: ClockSet,
  pub time_length: Option<f64>,
  #[default = 1.0]
  pub gamma: f64,
}

/// The clock inference state for a whole tree, routed through the clock pipeline as the per-node
/// [`ClockNodeState`] and per-edge [`ClockEdgeState`] fields, seeded from the timetree state's dates.
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
  pub fn new<D>(graph: &Graph<D>) -> Self
  where
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

  /// Seed clock state from value inputs instead of the graph payload, for tests that drive the clock
  /// passes without populating node payloads.
  ///
  /// Each node's date comes from `times` (`None` for a missing key, matching a leaf without a date),
  /// reproducing the timetree state's likely-time selection at the clock seed point. The divergence and outlier flag
  /// start at their payload defaults (`0.0` and `false`), `bad_branch` starts false, the clock set
  /// starts default (the backward pass recomputes the root clock set before it is read), and every
  /// edge starts default.
  pub fn seed_from_values<D>(graph: &Graph<D>, times: &BTreeMap<GraphNodeKey, Option<f64>>) -> Self
  where
    D: Send + Sync,
  {
    let nodes = graph
      .get_nodes()
      .iter()
      .map(|node| {
        let key = node.read_arc().key();
        let state = ClockNodeState {
          clock_set: ClockSet::default(),
          div: 0.0,
          time: times.get(&key).copied().flatten(),
          bad_branch: false,
          is_outlier: false,
        };
        (key, state)
      })
      .collect();
    let edges = graph
      .get_edges()
      .iter()
      .map(|edge| (edge.read_arc().key(), ClockEdgeState::default()))
      .collect();
    Self { nodes, edges }
  }

  /// Rebuild the per-node and per-edge clock maps to match the current graph, sourcing each node's
  /// date from `times`.
  ///
  /// Rebuilds the maps to match the current graph, so it stays valid across a reroot or polytomy
  /// resolution that added or dropped nodes and edges. Each node's date comes from `times`, the clock
  /// set starts default, and `bad_branch` starts false. The divergence and outlier flag are the two
  /// fields that live only in the value: they are preserved from the previous state for nodes that
  /// survived, and default for nodes that are new. Every edge resets to default messages, which the
  /// following backward pass recomputes.
  ///
  /// Used in the refinement loop, where the date passes have refined the node times on the threaded
  /// [`TimetreeState`](crate::timetree::timetree_state::TimetreeState) value between clock calls, so
  /// the regression must read the refined dates from the value.
  ///
  /// `edge_inputs` carries each edge's solver-updated time length and relaxed-clock rate multiplier
  /// (also from the date state); the re-estimation reads them to convert time back to divergence.
  /// Every other seed path leaves these at their defaults, where the regression reads input branch
  /// lengths instead.
  ///
  /// [`reseed_transitional_from_payloads`]: ClockState::reseed_transitional_from_payloads
  pub fn reseed_transitional_from_times<D>(
    &mut self,
    graph: &Graph<D>,
    times: &BTreeMap<GraphNodeKey, Option<f64>>,
    edge_inputs: &BTreeMap<GraphEdgeKey, (Option<f64>, f64)>,
  ) where
    D: Send + Sync,
  {
    self.reseed_transitional(graph, |key| times.get(&key).copied().flatten());
    for (key, &(time_length, gamma)) in edge_inputs {
      if let Some(edge) = self.edges.get_mut(key) {
        edge.time_length = time_length;
        edge.gamma = gamma;
      }
    }
  }

  fn reseed_transitional<D, F>(&mut self, graph: &Graph<D>, time_of: F)
  where
    D: Send + Sync,
    F: Fn(GraphNodeKey) -> Option<f64>,
  {
    let nodes = graph
      .get_nodes()
      .iter()
      .map(|node| {
        let node = node.read_arc();
        let key = node.key();
        let (div, is_outlier) = self
          .nodes
          .get(&key)
          .map_or((0.0, false), |state| (state.div, state.is_outlier));
        let state = ClockNodeState {
          clock_set: ClockSet::default(),
          div,
          time: time_of(key),
          bad_branch: false,
          is_outlier,
        };
        (key, state)
      })
      .collect();
    let edges = graph
      .get_edges()
      .iter()
      .map(|edge| (edge.read_arc().key(), ClockEdgeState::default()))
      .collect();
    self.nodes = nodes;
    self.edges = edges;
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
  pub fn map_backward<D, F>(&mut self, graph: &Graph<D>, visit: F) -> Result<(), Report>
  where
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
  pub fn map_forward<D, F>(&mut self, graph: &Graph<D>, visit: F) -> Result<(), Report>
  where
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
  /// The date the regression reads for this node, taken from the clock state's node time.
  #[must_use]
  pub fn likely_time(&self) -> Option<f64> {
    self.time
  }
}
