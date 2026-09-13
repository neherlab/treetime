use crate::gtr::infer_gtr::common::MutationCounts;
use eyre::Report;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::LogLh;

/// The marginal pass surface every sequence representation provides, over its own role-typed per-node
/// and per-edge element types.
///
/// A representation supplies the two passes and the transition counts; the full update, the per-node
/// likelihood reads, and the log-likelihood reset are derived from them here, so the four-step update
/// exists once rather than once per representation. Stable inputs (`graph`, `branch_lengths`) are
/// borrowed and every pass returns new owned maps, so a failed pass leaves its inputs intact.
pub trait PartitionMarginalOps {
  type Node: MarginalNodeState;
  type Backward;
  type Forward;
  type Estimate;

  /// Run the marginal backward pass (children before parent).
  fn marginal_backward(
    &self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: &BTreeMap<GraphNodeKey, Self::Node>,
  ) -> Result<MarginalBackward<Self::Node, Self::Backward>, Report>;

  /// Run the marginal forward pass (parent before children) over the backward messages.
  fn marginal_forward(
    &self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: &BTreeMap<GraphNodeKey, Self::Node>,
    backward: &BTreeMap<GraphEdgeKey, Self::Backward>,
  ) -> Result<MarginalForward<Self::Node, Self::Forward, Self::Estimate>, Report>;

  /// Count posterior-weighted state transitions over the tree, the input GTR inference reads.
  fn count_transitions(
    &self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    node_states: &BTreeMap<GraphNodeKey, Self::Node>,
    backward: &BTreeMap<GraphEdgeKey, Self::Backward>,
    forward: &BTreeMap<GraphEdgeKey, Self::Forward>,
  ) -> Result<MutationCounts, Report>;

  /// Run a full marginal update (backward, then forward) over the given node states.
  ///
  /// The substitution log likelihood is read at the root between the two passes: after the backward
  /// pass the root profile holds the likelihood of the observed data under the model, while the
  /// forward pass overwrites every node profile with its posterior.
  fn marginal_update(
    &self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: BTreeMap<GraphNodeKey, Self::Node>,
  ) -> Result<MarginalUpdate<Self::Node, Self::Backward, Self::Forward, Self::Estimate>, Report> {
    let MarginalBackward { node_states, backward } = self.marginal_backward(graph, branch_lengths, &node_states)?;
    let log_lh = self.root_log_lh(graph, &node_states)?;
    let MarginalForward {
      node_states,
      forward,
      estimates,
    } = self.marginal_forward(graph, branch_lengths, &node_states, &backward)?;
    Ok(MarginalUpdate {
      node_states,
      backward,
      forward,
      estimates,
      log_lh,
    })
  }

  /// The profile log likelihood recorded for one node, or zero when the node has no state.
  fn get_log_lh(&self, node_states: &BTreeMap<GraphNodeKey, Self::Node>, node_key: GraphNodeKey) -> LogLh {
    node_states
      .get(&node_key)
      .map_or(LogLh::ZERO, MarginalNodeState::log_lh)
  }

  /// The profile log likelihood recorded at the tree root.
  fn root_log_lh(&self, graph: &Graph, node_states: &BTreeMap<GraphNodeKey, Self::Node>) -> Result<LogLh, Report> {
    let root_key = graph.get_exactly_one_root()?.read_arc().key();
    Ok(self.get_log_lh(node_states, root_key))
  }

  /// Zero every node's profile log likelihood before a backward pass whose result is read as a
  /// likelihood.
  ///
  /// The backward pass reads a leaf's profile log likelihood as its message to the parent and folds it
  /// up the tree into the root log likelihood. A leaf profile carried over from an earlier forward
  /// pass holds that pass's posterior log likelihood, which adds a large model-independent constant to
  /// the result. That constant does not move the optimum of a rate search in exact arithmetic, but its
  /// magnitude erodes precision in Brent's parabolic interpolation and shifts the selected rate.
  fn reset_node_log_lh(&self, node_states: &mut BTreeMap<GraphNodeKey, Self::Node>) {
    for node in node_states.values_mut() {
      node.set_log_lh(LogLh::ZERO);
    }
  }
}

/// The result of one full marginal update (backward pass then forward pass): the refreshed per-node
/// states, the per-edge backward and forward messages, the per-edge estimates, and the substitution
/// log likelihood (the root node likelihood after the backward pass), as distinct owned values.
///
/// Generic over the backend element types. Sparse, dense, and discrete representations fill the same
/// map shape (one node map, three edge maps, one scalar) with their own node-state, message, and
/// estimate types, so a single named result serves all three and maps field-for-field onto the
/// reconstruction structs that own these maps.
#[derive(Clone, Debug, Serialize)]
pub struct MarginalUpdate<Node, Backward, Forward, Estimate> {
  pub node_states: BTreeMap<GraphNodeKey, Node>,
  pub backward: BTreeMap<GraphEdgeKey, Backward>,
  pub forward: BTreeMap<GraphEdgeKey, Forward>,
  pub estimates: BTreeMap<GraphEdgeKey, Estimate>,
  pub log_lh: LogLh,
}

/// The result of one marginal backward pass: the refreshed per-node states and the per-edge messages
/// toward the parent, as distinct owned values.
#[derive(Clone, Debug, Serialize)]
pub struct MarginalBackward<Node, Backward> {
  pub node_states: BTreeMap<GraphNodeKey, Node>,
  pub backward: BTreeMap<GraphEdgeKey, Backward>,
}

/// The result of one marginal forward pass: the refreshed per-node states, the per-edge messages
/// toward the child, and the per-edge estimates, as distinct owned values.
#[derive(Clone, Debug, Serialize)]
pub struct MarginalForward<Node, Forward, Estimate> {
  pub node_states: BTreeMap<GraphNodeKey, Node>,
  pub forward: BTreeMap<GraphEdgeKey, Forward>,
  pub estimates: BTreeMap<GraphEdgeKey, Estimate>,
}

/// A per-node marginal state that carries the node's profile log likelihood, so the shared update can
/// read and clear it without knowing the profile representation.
pub trait MarginalNodeState {
  fn log_lh(&self) -> LogLh;

  fn set_log_lh(&mut self, log_lh: LogLh);
}
