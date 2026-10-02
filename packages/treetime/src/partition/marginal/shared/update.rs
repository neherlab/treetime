use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::MutationCounts;
use eyre::Report;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::LogLh;

pub trait MarginalPasses {
  type Node: MarginalNodeState;
  type BackwardInput;
  type Backward;
  type Forward;
  type Estimate;

  fn backward_input(node_states: &BTreeMap<GraphNodeKey, Self::Node>) -> &Self::BackwardInput;

  fn backward_input_with_reset_log_lh(node_states: &BTreeMap<GraphNodeKey, Self::Node>) -> Self::BackwardInput;

  fn marginal_backward(
    &self,
    gtr: &GTR,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    input: &Self::BackwardInput,
  ) -> Result<MarginalBackward<Self::Node, Self::Backward>, Report>;

  fn marginal_forward(
    &self,
    gtr: &GTR,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: &BTreeMap<GraphNodeKey, Self::Node>,
    backward: &BTreeMap<GraphEdgeKey, Self::Backward>,
  ) -> Result<MarginalForward<Self::Node, Self::Forward, Self::Estimate>, Report>;

  fn count_transitions(
    &self,
    gtr: &GTR,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: &BTreeMap<GraphNodeKey, Self::Node>,
    backward: &BTreeMap<GraphEdgeKey, Self::Backward>,
    forward: &BTreeMap<GraphEdgeKey, Self::Forward>,
  ) -> Result<MutationCounts, Report>;

  fn marginal_update(
    &self,
    gtr: &GTR,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    input: &Self::BackwardInput,
  ) -> Result<MarginalUpdate<Self::Node, Self::Backward, Self::Forward, Self::Estimate>, Report> {
    let MarginalBackward { node_states, backward } = self.marginal_backward(gtr, graph, branch_lengths, input)?;
    let log_lh = self.root_log_lh(graph, &node_states)?;
    let MarginalForward {
      node_states,
      forward,
      estimates,
    } = self.marginal_forward(gtr, graph, branch_lengths, &node_states, &backward)?;
    Ok(MarginalUpdate {
      node_states,
      edges: MarginalEdges {
        backward,
        forward,
        estimates,
      },
      log_lh,
    })
  }

  fn get_log_lh(&self, node_states: &BTreeMap<GraphNodeKey, Self::Node>, node_key: GraphNodeKey) -> LogLh {
    node_states
      .get(&node_key)
      .map_or(LogLh::ZERO, MarginalNodeState::log_lh)
  }

  fn root_log_lh(&self, graph: &Graph, node_states: &BTreeMap<GraphNodeKey, Self::Node>) -> Result<LogLh, Report> {
    let root_key = graph.get_exactly_one_root()?.key();
    Ok(self.get_log_lh(node_states, root_key))
  }
}

#[derive(Clone, Debug, Serialize)]
pub struct MarginalUpdate<Node, Backward, Forward, Estimate> {
  pub(crate) node_states: BTreeMap<GraphNodeKey, Node>,
  pub(crate) edges: MarginalEdges<Backward, Forward, Estimate>,
  pub(crate) log_lh: LogLh,
}

#[derive(Clone, Debug, Serialize)]
pub struct MarginalEdges<Backward, Forward, Estimate> {
  pub backward: BTreeMap<GraphEdgeKey, Backward>,
  pub forward: BTreeMap<GraphEdgeKey, Forward>,
  pub estimates: BTreeMap<GraphEdgeKey, Estimate>,
}

impl<Backward, Forward, Estimate> Default for MarginalEdges<Backward, Forward, Estimate> {
  fn default() -> Self {
    Self {
      backward: BTreeMap::new(),
      forward: BTreeMap::new(),
      estimates: BTreeMap::new(),
    }
  }
}

#[derive(Clone, Debug, Serialize)]
pub struct MarginalBackward<Node, Backward> {
  pub(crate) node_states: BTreeMap<GraphNodeKey, Node>,
  pub(crate) backward: BTreeMap<GraphEdgeKey, Backward>,
}

#[derive(Clone, Debug, Serialize)]
pub struct MarginalForward<Node, Forward, Estimate> {
  pub(crate) node_states: BTreeMap<GraphNodeKey, Node>,
  pub(crate) forward: BTreeMap<GraphEdgeKey, Forward>,
  pub(crate) estimates: BTreeMap<GraphEdgeKey, Estimate>,
}

pub trait MarginalNodeState {
  fn log_lh(&self) -> LogLh;
}
