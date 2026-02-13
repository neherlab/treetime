use crate::commands::clock::clock_set::ClockSet;
use crate::commands::clock::clock_traits::ClockNode;
use crate::commands::timetree::timetree_traits::TimetreeNode;
use treetime_graph::edge::{EdgeOptimizeOps, GraphEdge};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use crate::seq::div::{OnlyLeaves, compute_divs};
use eyre::Report;

pub fn initialize_node_divergences<N, E, D>(graph: &Graph<N, E, D>)
where
  N: GraphNode + Named + ClockNode,
  E: EdgeOptimizeOps,
  D: Send + Sync,
{
  let divs = compute_divs(graph, OnlyLeaves(false));
  for node_ref in graph.get_nodes() {
    let mut node = node_ref.write_arc().payload().write_arc();
    let name = node.name().map(|n| n.as_ref().to_owned());
    if let Some(name) = name {
      if let Some(&div) = divs.get(&name) {
        node.set_div(div);
      }
    }
  }
}

pub fn initialize_clock_totals_from_time_distributions<N, E, D>(graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode + TimetreeNode + ClockNode,
  E: GraphEdge,
  D: Send + Sync,
{
  for node_ref in graph.get_nodes() {
    let mut node = node_ref.write_arc().payload().write_arc();
    if let Some(dist_arc) = node.time_distribution() {
      if let Some(time) = dist_arc.likely_time() {
        *node.clock_set_mut() = ClockSet::leaf_contribution(Some(time));
      }
    }
  }

  Ok(())
}
