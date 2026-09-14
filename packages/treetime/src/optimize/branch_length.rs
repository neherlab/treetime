use crate::make_report;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

/// Whether a scalar is in the physical domain of a phylogenetic branch length.
pub fn is_valid_branch_length_value(branch_length: f64) -> bool {
  branch_length.is_finite() && branch_length >= 0.0
}

/// Whether an optional branch length is present and in the physical domain.
pub fn is_valid_branch_length(branch_length: Option<f64>) -> bool {
  branch_length.is_some_and(is_valid_branch_length_value)
}

/// Return user-facing descriptions of all invalid branch lengths in graph order.
pub fn invalid_branch_length_descriptions(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<Vec<String>, Report> {
  graph
    .get_edges()
    .filter_map(|edge_ref| {
      let edge = edge_ref;
      let branch_length = branch_lengths[&edge.key()];
      (!is_valid_branch_length(branch_length)).then_some((edge.source(), edge.target(), branch_length))
    })
    .map(|(source, target, branch_length)| {
      let source = node_label(names, source);
      let target = node_label(names, target);
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

fn node_label(names: &BTreeMap<GraphNodeKey, Option<String>>, key: GraphNodeKey) -> String {
  names[&key].clone().unwrap_or_else(|| format!("node {key}"))
}
