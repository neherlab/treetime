#![allow(
  clippy::panic,
  reason = "mandated crash on missing graph node/edge access"
)]

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

/// Per-node durable clock inputs, held as a value keyed by [`GraphNodeKey`].
///
/// `time` is the observed or estimated date the regression reads for a leaf; `bad_branch` the
/// exclusion flag set from parsimony/date assignment. These are inputs the passes read but never
/// overwrite: the regression accumulators and fitted results live in [`ClockNodeState`].
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ClockNodeInput {
  pub time: Option<f64>,
  pub bad_branch: bool,
}

/// Per-edge durable clock inputs, held as a value keyed by [`GraphEdgeKey`].
///
/// `time_length` and `gamma` are the branch's solver-updated duration and relaxed-clock rate
/// multiplier the re-estimation reads to convert time back to divergence; they are seeded from the
/// date state in the refinement loop and stay at their defaults elsewhere, where the regression reads
/// input branch lengths instead.
#[derive(Debug, Clone, SmartDefault, PartialEq)]
pub struct ClockEdgeInput {
  pub time_length: Option<f64>,
  #[default = 1.0]
  pub gamma: f64,
}

/// Per-node clock regression accumulators and fitted results, keyed by [`GraphNodeKey`].
///
/// `clock_set` is the accumulated root-to-tip moment sums the backward pass recomputes from scratch;
/// `div` the cumulative divergence; `is_outlier` the fitted exclusion flag the clock filter writes.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ClockNodeState {
  pub clock_set: ClockSet,
  pub div: f64,
  pub is_outlier: bool,
}

/// Per-edge clock messages, held as a value keyed by [`GraphEdgeKey`]. The three clock messages are
/// recomputed by the regression passes and re-oriented on reroot.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ClockEdgeState {
  pub clock_to_parent: ClockSet,
  pub clock_to_child: ClockSet,
  /// The propagated `to_parent` message, kept to avoid recomputing the propagated message.
  pub clock_from_child: ClockSet,
}

/// Durable clock inputs for a whole tree: the observed dates, exclusion flags, and (in the refinement
/// loop) the solver-updated per-edge time lengths and relaxed-clock rate multipliers.
///
/// Kept as a separate owner from the produced [`ClockState`] regression results so the passes read
/// inputs they never overwrite. Keyed by stable node and edge ids, so the maps stay valid across a
/// reroot, which adds a split node, drops a trivial node, and re-orients the inverted path while
/// leaving ids stable (gaps, never renumbered).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ClockInputs {
  pub nodes: BTreeMap<GraphNodeKey, ClockNodeInput>,
  pub edges: BTreeMap<GraphEdgeKey, ClockEdgeInput>,
}

impl ClockInputs {
  /// Empty per-node/per-edge inputs for every node and edge of `graph`, all fields default.
  pub fn new(graph: &Graph) -> Self {
    let nodes = graph
      .get_nodes()
      .map(|node| (node.key(), ClockNodeInput::default()))
      .collect();
    let edges = graph
      .get_edges()
      .map(|edge| (edge.key(), ClockEdgeInput::default()))
      .collect();
    Self { nodes, edges }
  }

  /// Seed the per-node observed dates from `times` (`None` for a missing key, matching a leaf without
  /// a date), leaving every other input at its default and every edge default. Reproduces the
  /// timetree state's likely-time selection at the clock seed point.
  pub fn seed_from_times(graph: &Graph, times: &BTreeMap<GraphNodeKey, Option<f64>>) -> Self {
    let mut inputs = Self::new(graph);
    for node in graph.get_nodes() {
      let key = node.key();
      inputs.nodes.entry(key).or_default().time = times.get(&key).copied().flatten();
    }
    inputs
  }

  /// Rebuild the per-node and per-edge input maps to match the current graph, sourcing each node's
  /// date from `times` and each edge's solver-updated time length and relaxed-clock rate multiplier
  /// from `edge_inputs`.
  ///
  /// Stays valid across a reroot or polytomy resolution that added or dropped nodes and edges. Each
  /// node's date comes from `times`, `bad_branch` starts false. Every edge sources its time length and
  /// gamma from `edge_inputs`, defaulting where absent. Every other seed path leaves the edge inputs
  /// at their defaults, where the regression reads input branch lengths instead.
  pub fn reseed_from_times(
    &mut self,
    graph: &Graph,
    times: &BTreeMap<GraphNodeKey, Option<f64>>,
    edge_inputs: &BTreeMap<GraphEdgeKey, (Option<f64>, f64)>,
  ) {
    let nodes = graph
      .get_nodes()
      .map(|node| {
        let key = node.key();
        (
          key,
          ClockNodeInput {
            time: times.get(&key).copied().flatten(),
            bad_branch: false,
          },
        )
      })
      .collect();
    let edges = graph
      .get_edges()
      .map(|edge| {
        let key = edge.key();
        let (time_length, gamma) = edge_inputs.get(&key).copied().unwrap_or((None, 1.0));
        (key, ClockEdgeInput { time_length, gamma })
      })
      .collect();
    self.nodes = nodes;
    self.edges = edges;
  }

  #[must_use]
  pub fn node(&self, key: GraphNodeKey) -> &ClockNodeInput {
    self
      .nodes
      .get(&key)
      .unwrap_or_else(|| panic!("Clock inputs are missing node {key}"))
  }

  #[must_use]
  pub fn node_mut(&mut self, key: GraphNodeKey) -> &mut ClockNodeInput {
    self
      .nodes
      .get_mut(&key)
      .unwrap_or_else(|| panic!("Clock inputs are missing node {key}"))
  }

  #[must_use]
  pub fn edge(&self, key: GraphEdgeKey) -> &ClockEdgeInput {
    self
      .edges
      .get(&key)
      .unwrap_or_else(|| panic!("Clock inputs are missing edge {key}"))
  }

  /// The date the regression reads for a node, taken from the clock inputs.
  #[must_use]
  pub fn likely_time(&self, key: GraphNodeKey) -> Option<f64> {
    self.node(key).time
  }
}

/// The clock regression results for a whole tree, routed through the clock pipeline as the per-node
/// [`ClockNodeState`] and per-edge [`ClockEdgeState`] fields.
///
/// Produced by the regression passes and read back by the reroot search, the clock filter, and the
/// output gather. Keyed by stable node and edge ids, so the maps stay valid across a reroot.
#[derive(Debug, Clone, Default)]
pub struct ClockState {
  pub nodes: BTreeMap<GraphNodeKey, ClockNodeState>,
  pub edges: BTreeMap<GraphEdgeKey, ClockEdgeState>,
}

impl ClockState {
  /// Empty per-node/per-edge regression results for every node and edge of `graph`, all fields
  /// default.
  pub fn new(graph: &Graph) -> Self {
    let nodes = graph
      .get_nodes()
      .map(|node| (node.key(), ClockNodeState::default()))
      .collect();
    let edges = graph
      .get_edges()
      .map(|edge| (edge.key(), ClockEdgeState::default()))
      .collect();
    Self { nodes, edges }
  }

  /// Rebuild the per-node and per-edge result maps to match the current graph, preserving the
  /// value-resident divergence and outlier flag for surviving nodes.
  ///
  /// Rebuilds the maps to match the current graph, so it stays valid across a reroot or polytomy
  /// resolution that added or dropped nodes and edges. The clock set starts default (the following
  /// backward pass recomputes it). The divergence and outlier flag are preserved from the previous
  /// state for nodes that survived, and default for nodes that are new. Every edge resets to default
  /// messages, which the following backward pass recomputes.
  ///
  /// Used in the refinement loop, where the clock inputs have been reseeded from the refined date
  /// state between clock calls.
  pub fn reseed_transitional(&mut self, graph: &Graph) {
    let nodes = graph
      .get_nodes()
      .map(|node| {
        let key = node.key();
        let (div, is_outlier) = self
          .nodes
          .get(&key)
          .map_or((0.0, false), |state| (state.div, state.is_outlier));
        (
          key,
          ClockNodeState {
            clock_set: ClockSet::default(),
            div,
            is_outlier,
          },
        )
      })
      .collect();
    let edges = graph
      .get_edges()
      .map(|edge| (edge.key(), ClockEdgeState::default()))
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
  /// engine, replacing the per-node and per-edge result maps with the visitor's outputs.
  ///
  /// The engine reads the node/edge results from the current state and uses a
  /// thread-independent, deterministic child fold order.
  pub fn map_backward<F>(&mut self, graph: &Graph, visit: F) -> Result<(), Report>
  where
    F: Fn(
        GraphPassBackwardContext<'_, ClockNodeState, ClockEdgeState, ClockNodeState, ClockEdgeState>,
      ) -> Result<GraphPassNodeOutput<ClockNodeState, ClockEdgeState>, Report>
      + Sync
      + Send,
  {
    let pass = GraphPass::new(graph)?;
    let GraphMapOutputs { nodes, edges } =
      pass.map_backward(&self.nodes, &self.edges, |_| Ok(ClockNodeState::default()), visit)?;
    self.nodes = nodes;
    self.edges = edges;
    Ok(())
  }

  /// Run a value-returning forward pass over the clock state through the graph's dependency engine,
  /// replacing the per-node and per-edge result maps with the visitor's outputs. Each node reads its
  /// single parent's already-published output.
  pub fn map_forward<F>(&mut self, graph: &Graph, visit: F) -> Result<(), Report>
  where
    F: Fn(
        GraphPassForwardContext<'_, ClockNodeState, ClockEdgeState, ClockNodeState>,
      ) -> Result<GraphPassNodeOutput<ClockNodeState, ClockEdgeState>, Report>
      + Sync
      + Send,
  {
    let pass = GraphPass::new(graph)?;
    let GraphMapOutputs { nodes, edges } =
      pass.map_forward(&self.nodes, &self.edges, |_| Ok(ClockNodeState::default()), visit)?;
    self.nodes = nodes;
    self.edges = edges;
    Ok(())
  }
}
