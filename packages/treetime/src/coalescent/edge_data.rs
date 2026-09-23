use crate::coalescent::coalescent::CoalescentModel;
use crate::coalescent::node_time::{CoalescentNodeTime, CoalescentNodeTimes};
use crate::coalescent::time_coordinate::CalendarTime;
use eyre::Report;
use log::warn;
use treetime_graph::graph::Graph;
use treetime_primitives::LogLh;
use treetime_utils::make_error;

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub fn collect_coalescent_edges(
  graph: &Graph,
  node_times: &CoalescentNodeTimes,
) -> Result<Vec<CoalescentEdgeData>, Report> {
  let mut edges = Vec::new();

  graph.iter_breadth_first_forward(|node| {
    if node.parent_keys.is_empty() {
      return Ok(());
    }
    if node_times.get(&node.key).is_some_and(|entry| entry.bad_branch) {
      warn!(
        "Coalescent edge data: skipping node (key={:?}) with a bad branch",
        node.key
      );
      return Ok(());
    }

    let Some(child_time) = node_times.get(&node.key).and_then(node_time) else {
      warn!(
        "Coalescent edge data: skipping node (key={:?}) without an inferred date",
        node.key
      );
      return Ok(());
    };
    let parent_node_key = node.parent_keys[0].0;
    let Some(parent_time) = node_times.get(&parent_node_key).and_then(node_time) else {
      warn!(
        "Coalescent edge data: skipping node (key={:?}) whose parent has no inferred date",
        node.key
      );
      return Ok(());
    };

    if child_time < parent_time {
      return make_error!(
        "Coalescent edge has child older than parent: child key={:?}, child={child_time:.6e}, parent={parent_time:.6e}",
        node.key
      );
    }

    let n_siblings = graph
      .get_node(parent_node_key)
      .map_or(2.0, |parent| parent.outbound().len() as f64);
    edges.push(CoalescentEdgeData::new(
      CalendarTime::new(child_time),
      CalendarTime::new(parent_time),
      n_siblings,
    ));
    Ok(())
  })?;

  Ok(edges)
}

fn node_time(entry: &CoalescentNodeTime) -> Option<f64> {
  entry.time.or(entry.time_dist_likely)
}

pub fn coalescent_log_likelihood(edges: &[CoalescentEdgeData], model: &CoalescentModel) -> Result<LogLh, Report> {
  let total_contribution = edges
    .iter()
    .map(|edge| model.edge_contribution(edge))
    .sum::<Result<f64, Report>>()?;
  Ok(LogLh::new(-total_contribution))
}

#[derive(Clone, Debug)]
pub struct CoalescentEdgeData {
  child_time: CalendarTime,
  parent_time: CalendarTime,
  n_siblings: f64,
}

impl CoalescentEdgeData {
  pub fn new(child_time: CalendarTime, parent_time: CalendarTime, n_siblings: f64) -> Self {
    Self {
      child_time,
      parent_time,
      n_siblings,
    }
  }

  pub fn child_time(&self) -> CalendarTime {
    self.child_time
  }

  pub fn parent_time(&self) -> CalendarTime {
    self.parent_time
  }

  pub fn n_siblings(&self) -> f64 {
    self.n_siblings
  }
}
