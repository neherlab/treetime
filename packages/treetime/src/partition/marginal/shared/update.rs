use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
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
