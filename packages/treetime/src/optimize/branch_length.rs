use crate::{make_internal_report, make_report};
use eyre::Report;
use treetime_graph::edge::{GraphEdge, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNodeKey, NodeOptimizeOps};

/// Whether a scalar is in the physical domain of a phylogenetic branch length.
pub fn is_valid_branch_length_value(branch_length: f64) -> bool {
  branch_length.is_finite() && branch_length >= 0.0
}

/// Whether an optional branch length is present and in the physical domain.
pub fn is_valid_branch_length(branch_length: Option<f64>) -> bool {
  branch_length.is_some_and(is_valid_branch_length_value)
}

/// Return user-facing descriptions of all invalid branch lengths in graph order.
pub fn invalid_branch_length_descriptions<N, E, D>(graph: &Graph<N, E, D>) -> Result<Vec<String>, Report>
where
  N: NodeOptimizeOps,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
{
  graph
    .get_edges()
    .iter()
    .filter_map(|edge_ref| {
      let edge = edge_ref.read_arc();
      let branch_length = edge.payload().read_arc().branch_length();
      (!is_valid_branch_length(branch_length)).then_some((edge.source(), edge.target(), branch_length))
    })
    .map(|(source, target, branch_length)| {
      let source = node_label(graph, source)?;
      let target = node_label(graph, target)?;
      let branch_length = branch_length.map_or_else(|| "missing".to_owned(), |value| value.to_string());
      Ok(format!("{source} -> {target}: {branch_length}"))
    })
    .collect()
}

/// Require a scalar to be in the physical branch-length domain.
pub fn validate_branch_length_value(branch_length: f64) -> Result<(), Report> {
  if is_valid_branch_length_value(branch_length) {
    Ok(())
  } else {
    Err(make_report!(
      "branch length must be finite and non-negative, got {branch_length}"
    ))
  }
}

fn node_label<N, E, D>(graph: &Graph<N, E, D>, key: GraphNodeKey) -> Result<String, Report>
where
  N: NodeOptimizeOps,
  E: GraphEdge,
  D: Send + Sync,
{
  let node = graph
    .get_node(key)
    .ok_or_else(|| make_internal_report!("Edge references missing node {key}"))?;
  let node = node.read_arc();
  let payload = node.payload().read_arc();
  Ok(
    payload
      .name()
      .map_or_else(|| format!("node {key}"), |name| name.as_ref().to_owned()),
  )
}
