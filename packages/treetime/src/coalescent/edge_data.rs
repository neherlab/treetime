use crate::coalescent::coalescent::CoalescentModel;
use crate::coalescent::time_coordinate::CalendarTime;
use crate::payload::traits::TimetreeNode;
use eyre::Report;
use log::warn;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_utils::make_error;

/// Per-edge inferred calendar dates and the parent node's child count.
#[derive(Clone, Debug)]
pub struct CoalescentEdgeData {
  child_time: CalendarTime,
  parent_time: CalendarTime,
  n_children: f64,
}

impl CoalescentEdgeData {
  pub fn new(child_time: CalendarTime, parent_time: CalendarTime, n_children: f64) -> Self {
    Self {
      child_time,
      parent_time,
      n_children,
    }
  }

  pub fn child_time(&self) -> CalendarTime {
    self.child_time
  }

  pub fn parent_time(&self) -> CalendarTime {
    self.parent_time
  }

  /// Number of children of the edge's parent node (merger events = `n_children - 1`).
  pub fn n_children(&self) -> f64 {
    self.n_children
  }
}

/// Collects inferred child and parent dates for all non-root edges.
pub fn collect_coalescent_edges<N, E, D>(graph: &Graph<N, E, D>) -> Result<Vec<CoalescentEdgeData>, Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge,
  D: Sync + Send,
{
  let mut edges = Vec::new();

  graph.iter_breadth_first_forward(|node| {
    if node.parent_keys.is_empty() {
      return Ok(());
    }

    let Some(child_time) = node
      .payload
      .time_distribution()
      .as_ref()
      .and_then(|distribution| distribution.likely_time())
    else {
      warn!(
        "Coalescent edge data: skipping node (key={:?}) without an inferred date",
        node.key
      );
      return Ok(());
    };
    let parent_node_key = node.parent_keys[0].0;
    let Some(parent_time) = graph.get_node(parent_node_key).and_then(|parent| {
      parent
        .read_arc()
        .payload()
        .read_arc()
        .time_distribution()
        .as_ref()
        .and_then(|distribution| distribution.likely_time())
    }) else {
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

    let n_children = graph
      .get_node(parent_node_key)
      .map_or(2.0, |parent| parent.read_arc().outbound().len() as f64);
    edges.push(CoalescentEdgeData::new(
      CalendarTime::new(child_time),
      CalendarTime::new(parent_time),
      n_children,
    ));
    Ok(())
  })?;

  Ok(edges)
}

/// Sums the shared model's endpoint-derived edge contributions and negates to
/// return the coalescent log-likelihood (higher is more likely).
pub fn coalescent_log_likelihood(edges: &[CoalescentEdgeData], model: &CoalescentModel) -> Result<f64, Report> {
  let total_contribution = edges
    .iter()
    .map(|edge| model.edge_contribution(edge))
    .sum::<Result<f64, Report>>()?;
  Ok(-total_contribution)
}
