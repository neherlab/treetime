use crate::payload::traits::TimetreeNode;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey};

/// Per-node inferred times routed into the coalescent collectors as a value, in place of the
/// `NodeTimetree` payload's `time`/`time_distribution` fields.
///
/// Both sources are carried because the two collectors read the node's date differently:
/// [`collect_coalescent_edges`](crate::coalescent::edge_data::collect_coalescent_edges) prefers the
/// committed point estimate and falls back to the distribution peak, while
/// [`collect_tree_events`](crate::coalescent::events::collect_tree_events) reads only the
/// distribution peak.
#[derive(Clone, Copy, Debug, Default)]
pub struct CoalescentNodeTime {
  /// Committed point-estimate time, when one exists.
  pub time: Option<f64>,
  /// Likely time of the node's time distribution, when one exists.
  pub time_dist_likely: Option<f64>,
}

/// Per-node inferred times keyed by node, threaded into the coalescent collectors instead of read
/// off the graph payload.
pub type CoalescentNodeTimes = BTreeMap<GraphNodeKey, CoalescentNodeTime>;

/// Build the coalescent node-time map from the graph payloads.
///
/// Used by callers that still hold the node times on the payload (tests). The timetree pipeline
/// builds the map from its value state instead.
pub fn coalescent_node_times_from_payloads<N, E, D>(graph: &Graph<N, E, D>) -> CoalescentNodeTimes
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge,
  D: Sync + Send,
{
  graph
    .get_nodes()
    .iter()
    .map(|node| {
      let node = node.read_arc();
      let payload = node.payload().read_arc();
      let entry = CoalescentNodeTime {
        time: payload.time(),
        time_dist_likely: payload.time_distribution().as_ref().and_then(|dist| dist.likely_time()),
      };
      (node.key(), entry)
    })
    .collect()
}
