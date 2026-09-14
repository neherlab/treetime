use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::MutationCounts;
use eyre::Report;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::LogLh;

/// The result of one full marginal update (backward pass then forward pass): the refreshed per-node
/// states, the per-edge results of that update, and the substitution log likelihood (the root node
/// likelihood after the backward pass), as distinct owned values.
///
/// Generic over the backend element types. Sparse, dense, and discrete representations fill the same
/// map shape (one node map, three edge maps, one scalar) with their own node-state, message, and
/// estimate types, so a single named result serves all three and maps field-for-field onto the
/// reconstruction structs that own these maps.
#[derive(Clone, Debug, Serialize)]
pub struct MarginalUpdate<Node, Backward, Forward, Estimate> {
  pub node_states: BTreeMap<GraphNodeKey, Node>,
  pub edges: MarginalEdges<Backward, Forward, Estimate>,
  pub log_lh: LogLh,
}

/// The per-edge results of one marginal update: the backward messages toward the parent, the forward
/// messages toward the child, and the estimates derived during the forward pass.
///
/// Grouped because they share one lifetime. All three are produced by a single update, read by the
/// consumers of that update, and invalidated together by any change of topology or branch length. The
/// node states, in contrast, survive a structural change and seed the next update, so a caller that
/// keeps inference state across stages keeps the node states and drops this value.
#[derive(Clone, Debug, Serialize)]
pub struct MarginalEdges<Backward, Forward, Estimate> {
  pub backward: BTreeMap<GraphEdgeKey, Backward>,
  pub forward: BTreeMap<GraphEdgeKey, Forward>,
  pub estimates: BTreeMap<GraphEdgeKey, Estimate>,
}

impl<Backward, Forward, Estimate> Default for MarginalEdges<Backward, Forward, Estimate> {
  /// The empty per-edge result set: what a representation holds before its first update, and after a
  /// structural change has invalidated the previous update's edge results.
  fn default() -> Self {
    Self {
      backward: BTreeMap::new(),
      forward: BTreeMap::new(),
      estimates: BTreeMap::new(),
    }
  }
}

/// The result of a marginal update read for its node states alone: the refreshed per-node states and
/// the substitution log likelihood.
#[derive(Clone, Debug, Serialize)]
pub struct MarginalStates<Node> {
  pub node_states: BTreeMap<GraphNodeKey, Node>,
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

/// The marginal pass surface every sequence representation provides, over its own role-typed per-node
/// and per-edge element types.
///
/// A representation supplies model access, the two marginal passes, and transition counting. The full
/// update, the per-node likelihood reads, and the log-likelihood reset derive from those here, so the
/// four-step update and the likelihood plumbing exist once rather than once per representation. Stable
/// inputs (`graph`, branch lengths) are borrowed and every pass returns new owned maps, so a failed pass
/// leaves its inputs intact.
pub trait MarginalPasses {
  type Node: MarginalNodeState;
  type Backward;
  type Forward;
  type Estimate;

  /// The current substitution model.
  fn gtr(&self) -> &GTR;

  /// Replace the substitution model.
  fn set_gtr(&mut self, gtr: GTR);

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
  /// pass the root profile holds the likelihood of the observed data under the model, while the forward
  /// pass overwrites every node profile with its posterior.
  #[allow(clippy::needless_pass_by_value)]
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
      edges: MarginalEdges {
        backward,
        forward,
        estimates,
      },
      log_lh,
    })
  }

  /// Run a full marginal update and return only the refreshed node states and the substitution log
  /// likelihood, dropping the per-edge messages and estimates.
  fn marginal_states(
    &self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: BTreeMap<GraphNodeKey, Self::Node>,
  ) -> Result<MarginalStates<Self::Node>, Report> {
    let MarginalUpdate {
      node_states, log_lh, ..
    } = self.marginal_update(graph, branch_lengths, node_states)?;
    Ok(MarginalStates { node_states, log_lh })
  }

  /// The profile log likelihood recorded for one node, or zero when the node has no state.
  fn get_log_lh(&self, node_states: &BTreeMap<GraphNodeKey, Self::Node>, node_key: GraphNodeKey) -> LogLh {
    node_states
      .get(&node_key)
      .map_or(LogLh::ZERO, MarginalNodeState::log_lh)
  }

  /// The profile log likelihood recorded at the tree root, or zero when the root has no state.
  fn root_log_lh(&self, graph: &Graph, node_states: &BTreeMap<GraphNodeKey, Self::Node>) -> Result<LogLh, Report> {
    let root_key = graph.get_exactly_one_root()?.read_arc().key();
    Ok(self.get_log_lh(node_states, root_key))
  }

  /// Zero every node's profile log likelihood before a backward pass whose result is read as a
  /// likelihood.
  ///
  /// The backward pass reads a leaf's profile log likelihood as its message to the parent and folds it
  /// up the tree into the root log likelihood. A leaf profile carried over from an earlier forward pass
  /// holds that pass's posterior log likelihood, which adds a large model-independent constant to the
  /// result. That constant does not move the optimum of a rate search in exact arithmetic, but its
  /// magnitude erodes precision in Brent's parabolic interpolation and shifts the selected rate.
  fn reset_node_log_lh(&self, node_states: &mut BTreeMap<GraphNodeKey, Self::Node>) {
    for node in node_states.values_mut() {
      node.set_log_lh(LogLh::ZERO);
    }
  }
}
