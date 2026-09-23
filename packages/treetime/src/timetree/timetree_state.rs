use crate::clock::date_constraints::DateConstraints;
use crate::coalescent::node_time::{CoalescentNodeTime, CoalescentNodeTimes};
use eyre::Report;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::{Distribution, NegLog};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pass::{
  GraphMapOutputs, GraphPass, GraphPassBackwardContext, GraphPassForwardContext, GraphPassNodeOutput,
};

#[derive(Debug, Clone, Default)]
pub struct TimetreeState {
  pub nodes: BTreeMap<GraphNodeKey, DateNodeState>,
  pub edges: BTreeMap<GraphEdgeKey, DateEdgeState>,
}

impl TimetreeState {
  pub fn new(graph: &Graph) -> Self {
    let nodes = graph
      .get_nodes()
      .map(|node| (node.key(), DateNodeState::default()))
      .collect();
    let edges = graph
      .get_edges()
      .map(|edge| (edge.key(), DateEdgeState::default()))
      .collect();
    Self { nodes, edges }
  }

  pub(crate) fn seed_from_values(graph: &Graph, constraints: &DateConstraints) -> Self {
    let nodes = graph
      .get_nodes()
      .map(|node| {
        let key = node.key();
        let state = DateNodeState {
          time_distribution: constraints.time_distributions.get(&key).cloned().flatten(),
          time: None,
          bad_branch: constraints.bad_branches.get(&key).copied().unwrap_or(false),
          contradicted: false,
        };
        (key, state)
      })
      .collect();
    let edges = graph
      .get_edges()
      .map(|edge| (edge.key(), DateEdgeState::default()))
      .collect();
    Self { nodes, edges }
  }

  pub(crate) fn reseed_from_values(&mut self, graph: &Graph) {
    let nodes = graph
      .get_nodes()
      .map(|node| {
        let key = node.key();
        let state = self
          .nodes
          .get(&key)
          .map_or_else(DateNodeState::default, |node| DateNodeState {
            time_distribution: node.time_distribution.clone(),
            time: node.time,
            bad_branch: node.bad_branch,
            contradicted: false,
          });
        (key, state)
      })
      .collect();
    let edges = graph
      .get_edges()
      .map(|edge| {
        let key = edge.key();
        let state = self
          .edges
          .get(&key)
          .map_or_else(DateEdgeState::default, |edge| DateEdgeState {
            branch_length_distribution: edge.branch_length_distribution.clone(),
            msg_to_parent: edge.msg_to_parent.clone(),
            time_length: edge.time_length,
            gamma: edge.gamma,
          });
        (key, state)
      })
      .collect();
    self.nodes = nodes;
    self.edges = edges;
  }

  pub(crate) fn reset_date_edges_for_topology_change(&mut self, graph: &Graph) {
    for edge_ref in graph.get_edges() {
      let key = edge_ref.key();
      let entry = self.edges.entry(key).or_default();
      entry.branch_length_distribution = None;
      entry.msg_to_parent = None;
      entry.gamma = 1.0;
    }
    for node_ref in graph.get_nodes() {
      let key = node_ref.key();
      self.nodes.entry(key).or_default();
    }
  }

  #[must_use]
  pub(crate) fn likely_times(&self, constraints: &DateConstraints) -> BTreeMap<GraphNodeKey, Option<f64>> {
    self
      .nodes
      .iter()
      .map(|(key, node)| {
        let time = constraints
          .date_constraints
          .get(key)
          .cloned()
          .flatten()
          .as_ref()
          .or(node.time_distribution.as_ref())
          .and_then(|dist| dist.likely_time());
        (*key, time)
      })
      .collect()
  }

  #[must_use]
  pub(crate) fn coalescent_node_times(&self) -> CoalescentNodeTimes {
    self
      .nodes
      .iter()
      .map(|(key, node)| {
        let entry = CoalescentNodeTime {
          time: node.time,
          time_dist_likely: node.time_distribution.as_ref().and_then(|dist| dist.likely_time()),
          bad_branch: node.bad_branch,
        };
        (*key, entry)
      })
      .collect()
  }

  #[allow(clippy::panic, reason = "panics on a violated internal invariant")]
  #[must_use]
  pub fn node(&self, key: GraphNodeKey) -> &DateNodeState {
    self
      .nodes
      .get(&key)
      .unwrap_or_else(|| panic!("Timetree state is missing node {key}"))
  }

  #[allow(clippy::panic, reason = "panics on a violated internal invariant")]
  #[must_use]
  pub fn node_mut(&mut self, key: GraphNodeKey) -> &mut DateNodeState {
    self
      .nodes
      .get_mut(&key)
      .unwrap_or_else(|| panic!("Timetree state is missing node {key}"))
  }

  #[allow(clippy::panic, reason = "panics on a violated internal invariant")]
  #[must_use]
  pub fn edge(&self, key: GraphEdgeKey) -> &DateEdgeState {
    self
      .edges
      .get(&key)
      .unwrap_or_else(|| panic!("Timetree state is missing edge {key}"))
  }

  #[allow(clippy::panic, reason = "panics on a violated internal invariant")]
  #[must_use]
  pub(crate) fn edge_mut(&mut self, key: GraphEdgeKey) -> &mut DateEdgeState {
    self
      .edges
      .get_mut(&key)
      .unwrap_or_else(|| panic!("Timetree state is missing edge {key}"))
  }

  pub(crate) fn map_backward<F>(&mut self, graph: &Graph, visit: F) -> Result<(), Report>
  where
    F: Fn(
        GraphPassBackwardContext<'_, DateNodeState, DateEdgeState, DateNodeState, DateEdgeState>,
      ) -> Result<GraphPassNodeOutput<DateNodeState, DateEdgeState>, Report>
      + Sync
      + Send,
  {
    let pass = GraphPass::new(graph)?;
    let GraphMapOutputs { nodes, edges } =
      pass.map_backward(&self.nodes, &self.edges, |_| Ok(DateNodeState::default()), visit)?;
    self.nodes = nodes;
    self.edges = edges;
    Ok(())
  }

  pub(crate) fn map_forward<F>(&mut self, graph: &Graph, visit: F) -> Result<(), Report>
  where
    F: Fn(
        GraphPassForwardContext<'_, DateNodeState, DateEdgeState, DateNodeState>,
      ) -> Result<GraphPassNodeOutput<DateNodeState, DateEdgeState>, Report>
      + Sync
      + Send,
  {
    let pass = GraphPass::new(graph)?;
    let GraphMapOutputs { nodes, edges } =
      pass.map_forward(&self.nodes, &self.edges, |_| Ok(DateNodeState::default()), visit)?;
    self.nodes = nodes;
    self.edges = edges;
    Ok(())
  }
}

#[derive(Debug, Clone, Default)]
pub struct DateNodeState {
  pub time_distribution: Option<Arc<Distribution<NegLog>>>,
  pub time: Option<f64>,
  pub bad_branch: bool,
  pub contradicted: bool,
}

#[derive(Debug, Clone, SmartDefault)]
pub struct DateEdgeState {
  pub branch_length_distribution: Option<Arc<Distribution<NegLog>>>,
  pub msg_to_parent: Option<Arc<Distribution<NegLog>>>,
  pub time_length: Option<f64>,
  #[default = 1.0]
  pub gamma: f64,
}
