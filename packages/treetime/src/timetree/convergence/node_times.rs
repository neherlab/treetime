use crate::timetree::timetree_state::TimetreeState;
use std::collections::BTreeMap;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

pub(crate) fn capture_node_times(graph: &Graph, state: &TimetreeState) -> NodeTimeSnapshot {
  graph
    .get_nodes()
    .filter_map(|node| {
      let key = node.key();
      let time = state.node(key).time?;
      time.is_finite().then_some((key, time))
    })
    .collect()
}

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub(crate) fn measure_node_time_change(previous: &NodeTimeSnapshot, current: &NodeTimeSnapshot) -> NodeTimeChange {
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

pub(crate) type NodeTimeSnapshot = BTreeMap<GraphNodeKey, f64>;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NodeTimeChange {
  pub max: Option<f64>,
  pub rms: Option<f64>,
}
