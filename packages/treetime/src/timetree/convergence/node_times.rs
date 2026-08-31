use crate::partition::timetree::partition::GraphTimetree;
use crate::payload::traits::TimetreeNode;
use std::collections::BTreeMap;
use treetime_graph::node::GraphNodeKey;

/// Inferred time of every dated node, keyed by node.
pub type NodeTimeSnapshot = BTreeMap<GraphNodeKey, f64>;

/// How far node times moved between two snapshots, in years.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NodeTimeChange {
  /// Largest absolute movement of any single node. `None` when nothing is comparable.
  pub max: Option<f64>,
  /// Root mean square movement over the comparable nodes.
  pub rms: Option<f64>,
}

pub fn capture_node_times(graph: &GraphTimetree) -> NodeTimeSnapshot {
  graph
    .get_nodes()
    .iter()
    .filter_map(|node| {
      let node = node.read_arc();
      let time = node.payload().read_arc().time()?;
      time.is_finite().then(|| (node.key(), time))
    })
    .collect()
}

/// Compare two snapshots over the nodes present and dated in both.
///
/// Nodes appearing in only one snapshot are skipped: polytomy resolution introduces nodes that
/// have no earlier position to be compared against. A round that changed the topology is not
/// judged converged regardless, since `n_resolved` must also be zero.
///
/// Both statistics are reported because their ratio is diagnostic: `max` far above `rms` means a
/// single node is oscillating while the rest of the tree is still, which is a local bistability
/// rather than a loop that has failed to settle.
pub fn measure_node_time_change(previous: &NodeTimeSnapshot, current: &NodeTimeSnapshot) -> NodeTimeChange {
  let changes: Vec<f64> = previous
    .iter()
    .filter_map(|(key, prev)| current.get(key).map(|curr| (curr - prev).abs()))
    .collect();

  if changes.is_empty() {
    return NodeTimeChange::default();
  }

  let max = changes.iter().copied().fold(f64::NEG_INFINITY, f64::max);
  let rms = (changes.iter().map(|change| change * change).sum::<f64>() / changes.len() as f64).sqrt();
  NodeTimeChange {
    max: Some(max),
    rms: Some(rms),
  }
}
