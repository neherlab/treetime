use crate::commands::clock::clock_set::ClockSet;
use crate::representation::graph_ancestral::GraphAncestral;
use crate::seq::div::{OnlyLeaves, calculate_divs};
use eyre::Report;

pub fn initialize_node_divergences(graph: &GraphAncestral) {
  let divs = calculate_divs(graph, OnlyLeaves(false));
  for node_ref in graph.get_nodes() {
    let mut node = node_ref.write_arc().payload().write_arc();
    if let Some(name) = &node.name {
      if let Some(&div) = divs.get(name) {
        node.div = div;
      }
    }
  }
}

pub fn initialize_clock_totals_from_time_distributions(graph: &GraphAncestral) -> Result<(), Report> {
  for node_ref in graph.get_nodes() {
    let mut node = node_ref.write_arc().payload().write_arc();
    if let Some(dist_arc) = &node.time_distribution {
      if let Some(time) = dist_arc.likely_time() {
        node.clock_set = ClockSet::leaf_contribution(Some(time));
      }
    }
  }

  Ok(())
}
