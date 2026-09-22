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

#[derive(Debug, Clone, Default, PartialEq)]
pub struct ClockNodeInput {
  pub time: Option<f64>,
  pub bad_branch: bool,
}

#[derive(Debug, Clone, SmartDefault, PartialEq)]
pub struct ClockEdgeInput {
  pub time_length: Option<f64>,
  #[default = 1.0]
  pub gamma: f64,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct ClockNodeState {
  pub clock_set: ClockSet,
  pub div: f64,
  pub is_outlier: bool,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct ClockEdgeState {
  pub clock_to_parent: ClockSet,
  pub clock_to_child: ClockSet,
  pub clock_from_child: ClockSet,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct ClockInputs {
  pub nodes: BTreeMap<GraphNodeKey, ClockNodeInput>,
  pub edges: BTreeMap<GraphEdgeKey, ClockEdgeInput>,
}

impl ClockInputs {
  pub(crate) fn new(graph: &Graph) -> Self {
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

  pub(crate) fn seed_from_times(graph: &Graph, times: &BTreeMap<GraphNodeKey, Option<f64>>) -> Self {
    let mut inputs = Self::new(graph);
    for node in graph.get_nodes() {
      let key = node.key();
      inputs.nodes.entry(key).or_default().time = times.get(&key).copied().flatten();
    }
    inputs
  }

  pub(crate) fn reseed_from_times(
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

  #[allow(clippy::panic, reason = "panics on a violated internal invariant")]
  #[must_use]
  pub fn node(&self, key: GraphNodeKey) -> &ClockNodeInput {
    self
      .nodes
      .get(&key)
      .unwrap_or_else(|| panic!("Clock inputs are missing node {key}"))
  }

  #[allow(clippy::panic, reason = "panics on a violated internal invariant")]
  #[must_use]
  pub(crate) fn node_mut(&mut self, key: GraphNodeKey) -> &mut ClockNodeInput {
    self
      .nodes
      .get_mut(&key)
      .unwrap_or_else(|| panic!("Clock inputs are missing node {key}"))
  }

  #[allow(clippy::panic, reason = "panics on a violated internal invariant")]
  #[must_use]
  pub(crate) fn edge(&self, key: GraphEdgeKey) -> &ClockEdgeInput {
    self
      .edges
      .get(&key)
      .unwrap_or_else(|| panic!("Clock inputs are missing edge {key}"))
  }

  #[must_use]
  pub(crate) fn likely_time(&self, key: GraphNodeKey) -> Option<f64> {
    self.node(key).time
  }
}

#[derive(Debug, Clone, Default)]
pub struct ClockState {
  pub nodes: BTreeMap<GraphNodeKey, ClockNodeState>,
  pub edges: BTreeMap<GraphEdgeKey, ClockEdgeState>,
}

impl ClockState {
  pub(crate) fn new(graph: &Graph) -> Self {
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

  pub(crate) fn reseed_transitional(&mut self, graph: &Graph) {
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

  #[allow(clippy::panic, reason = "panics on a violated internal invariant")]
  #[must_use]
  pub fn node(&self, key: GraphNodeKey) -> &ClockNodeState {
    self
      .nodes
      .get(&key)
      .unwrap_or_else(|| panic!("Clock state is missing node {key}"))
  }

  #[allow(clippy::panic, reason = "panics on a violated internal invariant")]
  #[must_use]
  pub(crate) fn node_mut(&mut self, key: GraphNodeKey) -> &mut ClockNodeState {
    self
      .nodes
      .get_mut(&key)
      .unwrap_or_else(|| panic!("Clock state is missing node {key}"))
  }

  #[allow(clippy::panic, reason = "panics on a violated internal invariant")]
  #[must_use]
  pub(crate) fn edge(&self, key: GraphEdgeKey) -> &ClockEdgeState {
    self
      .edges
      .get(&key)
      .unwrap_or_else(|| panic!("Clock state is missing edge {key}"))
  }

  pub(crate) fn map_backward<F>(&mut self, graph: &Graph, visit: F) -> Result<(), Report>
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

  pub(crate) fn map_forward<F>(&mut self, graph: &Graph, visit: F) -> Result<(), Report>
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
