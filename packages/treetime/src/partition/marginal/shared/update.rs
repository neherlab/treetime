use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::LogLh;

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
