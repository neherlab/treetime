use std::collections::BTreeMap;
use treetime_graph::node::GraphNodeKey;

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
  /// Whether the node is excluded as a bad branch. The collectors skip a node with this set, reading
  /// it from here instead of off the graph payload.
  pub bad_branch: bool,
}

/// Per-node inferred times keyed by node, threaded into the coalescent collectors instead of read
/// off the graph payload.
pub type CoalescentNodeTimes = BTreeMap<GraphNodeKey, CoalescentNodeTime>;
