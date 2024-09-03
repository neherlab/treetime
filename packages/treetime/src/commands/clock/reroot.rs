use crate::commands::clock::clock_graph::{ClockEdgePayload, ClockGraph, ClockNodePayload};
use crate::commands::clock::clock_regression::ClockOptions;
use crate::commands::clock::clock_set::ClockSet;
use crate::commands::clock::find_best_root::{find_best_root, FindRootResult};
use crate::graph::edge::{invert_edge, GraphEdgeKey, Weighted};
use crate::graph::node::GraphNodeKey;
use approx::ulps_eq;
use eyre::Report;

pub fn reroot_in_place(graph: &mut ClockGraph, options: &ClockOptions) -> Result<GraphNodeKey, Report> {
  let old_root_key = { graph.get_exactly_one_root()?.read_arc().key() };

  let FindRootResult {
    edge,
    split,
    total,
    chisq,
  } = find_best_root(graph, options)?;
  // if edge is null, we are already at the best root, return old_root_key
  if edge.is_none() {
    return Ok(old_root_key);
  }

  let edge_key = edge.expect("Edge is empty when rerooting");
  let edge = graph.get_edge(edge_key).expect("Edge not found");

  let new_root_key = if ulps_eq!(split, 0.0, max_ulps = 5) {
    edge.read_arc().target()
  } else if ulps_eq!(split, 1.0, max_ulps = 5) {
    edge.read_arc().source()
  } else {
    create_new_root_node(graph, edge_key, split, total)?
  };

  if new_root_key != old_root_key {
    apply_reroot(graph, old_root_key, new_root_key, options)?;
  }

  // TODO: remove old root node if it's trivial (i.e. having exactly 1 child and 1 parent) and merge dangling edges

  Ok(new_root_key)
}

/// Create new root node by splitting the edge into two
fn create_new_root_node(
  graph: &mut ClockGraph,
  edge_key: GraphEdgeKey,
  split: f64,
  total: ClockSet,
) -> Result<GraphNodeKey, Report> {
  let new_root_key = graph.add_node(ClockNodePayload {
    name: Some("new_root".to_owned()),
    date: None,
    bad_branch: false,
    div: 0.0,
    is_outlier: false,
    total: total.clone(),
  });

  let edge = graph.get_edge(edge_key).expect("Edge not found");
  let source_key = edge.read_arc().source();
  let target_key = edge.read_arc().target();
  let branch_length = edge.read_arc().payload().read_arc().weight().unwrap_or_default();

  graph.add_edge(
    source_key,
    new_root_key,
    ClockEdgePayload {
      branch_length: Some(split * branch_length),
      ..ClockEdgePayload::default()
    },
  )?;

  graph.add_edge(
    new_root_key,
    target_key,
    ClockEdgePayload {
      branch_length: Some((1.0 - split) * branch_length),
      ..ClockEdgePayload::default()
    },
  )?;

  graph.remove_edge(edge_key)?;

  Ok(new_root_key)
}

/// Modify graph topology to make the newly identified root the actual root.
fn apply_reroot(
  graph: &mut ClockGraph,
  old_root_key: GraphNodeKey,
  new_root_key: GraphNodeKey,
  options: &ClockOptions,
) -> Result<(), Report> {
  // Find paths from the old root to the new desired root
  let paths = graph.path_from_node_to_node(new_root_key, old_root_key)?;

  // Invert every edge on the path from old to new root.
  // This will make the desired new root into an actual root. The old root might no longer be a root.
  for (_, edge) in &paths {
    if let Some(edge) = edge {
      invert_edge(graph, edge);
      let mut edge_payload = edge.read_arc().payload().write_arc();
      let edge_len = edge_payload.weight().unwrap();
      let branch_variance = options.variance_factor * edge_len + options.variance_offset;
      let tmp_to_parent = edge_payload.to_parent.clone();
      edge_payload.to_parent = edge_payload.to_child.clone();
      edge_payload.to_child = tmp_to_parent;
      edge_payload.from_child = edge_payload.to_parent.propagate_averages(edge_len, branch_variance);
    }
  }

  // Some bookkeeping
  graph.build()?;
  Ok(())
}
