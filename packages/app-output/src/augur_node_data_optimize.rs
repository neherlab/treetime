use crate::optimize_result::OptimizeNodeOut;
use eyre::Report;
use std::collections::BTreeMap;
use std::path::Path;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::io::json::{JsonPretty, json_write_file};
use util_augur_node_data_json::{
  AugurNodeDataJsonGeneratedBy, AugurNodeDataJsonRefine, AugurNodeDataJsonRefineMeta, AugurNodeDataJsonRefineNode,
};

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub fn build_augur_node_data_json(
  graph: &Graph,
  node_outputs: &BTreeMap<GraphNodeKey, OptimizeNodeOut>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  alignment: Option<&Path>,
  input_tree: Option<&Path>,
  mutation_counts: Option<&BTreeMap<GraphEdgeKey, usize>>,
) -> Result<AugurNodeDataJsonRefine, Report> {
  let mut nodes = BTreeMap::new();
  for node in graph.get_nodes() {
    let node_guard = node;
    let node_key = node_guard.key();
    let out = &node_outputs[&node_key];
    let node_name = out
      .name
      .as_deref()
      .map_or_else(|| format!("node_{}", node_key.0), str::to_owned);

    let branch_length = match graph.node_parent(node_key)? {
      Some((_parent_key, edge_key)) => {
        if let Some(counts) = mutation_counts {
          counts.get(&edge_key).copied().unwrap_or_default() as f64
        } else {
          branch_lengths[&edge_key].unwrap_or(0.0)
        }
      },
      None => 0.0,
    };

    let confidence = out.confidence;

    nodes.insert(
      node_name,
      AugurNodeDataJsonRefineNode {
        branch_length,
        confidence,
        numdate: None,
        clock_length: None,
        mutation_length: None,
        raw_date: None,
        date: None,
        date_inferred: None,
        num_date_confidence: None,
        other: BTreeMap::new(),
      },
    );
  }

  Ok(AugurNodeDataJsonRefine {
    generated_by: Some(AugurNodeDataJsonGeneratedBy {
      program: "treetime".to_owned(),
      version: env!("CARGO_PKG_VERSION").to_owned(),
    }),
    metadata: AugurNodeDataJsonRefineMeta {
      alignment: alignment.map(|path| path.display().to_string()),
      input_tree: input_tree.map(|path| path.display().to_string()),
      clock: None,
      other: BTreeMap::new(),
    },
    nodes,
  })
}

pub fn write_augur_node_data_json(
  graph: &Graph,
  node_outputs: &BTreeMap<GraphNodeKey, OptimizeNodeOut>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  alignment: Option<&Path>,
  input_tree: Option<&Path>,
  mutation_counts: Option<&BTreeMap<GraphEdgeKey, usize>>,
  path: &Path,
) -> Result<(), Report> {
  let data = build_augur_node_data_json(
    graph,
    node_outputs,
    branch_lengths,
    alignment,
    input_tree,
    mutation_counts,
  )?;
  json_write_file(path, &data, JsonPretty(true))?;
  Ok(())
}
