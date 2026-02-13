use eyre::Report;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey};

pub trait HasLogLh {
  /// Get the log likelihood contribution for a given node
  fn get_log_lh(&self, node_key: GraphNodeKey) -> f64;
}

/// Calculate the total log likelihood of the graph given the partitions
pub fn graph_log_lh<P, N, E, D>(graph: &Graph<N, E, D>, partitions: &[Arc<RwLock<P>>]) -> Result<f64, Report>
where
  P: HasLogLh + ?Sized,
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send + Default,
{
  let root = graph.get_exactly_one_root()?;
  let root_key = root.read_arc().key();

  let log_lh = partitions
    .iter()
    .map(|partition| partition.read_arc().get_log_lh(root_key))
    .sum();

  Ok(log_lh)
}
