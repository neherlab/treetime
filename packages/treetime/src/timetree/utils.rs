use crate::clock::clock_state::ClockState;
use crate::seq::div::{OnlyLeaves, compute_divs};
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

pub(crate) fn initialize_node_divergences(
  graph: &Graph,
  clock_state: &mut ClockState,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<(), Report> {
  let divs = compute_divs(graph, OnlyLeaves(false), branch_lengths, names)?;
  for node_ref in graph.get_nodes() {
    let node = node_ref;
    let key = node.key();
    if let Some(name) = &names[&key] {
      if let Some(&div) = divs.get(name) {
        clock_state.nodes.entry(key).or_default().div = div;
      }
    }
  }
  Ok(())
}
