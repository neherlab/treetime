use crate::partition::traits::BranchTopology;
use crate::payload::ancestral::GraphAncestral;
use eyre::Report;
use std::collections::BTreeMap;
use std::path::Path;
use treetime_graph::edge::HasBranchLength;
use treetime_graph::node::Named;
use treetime_utils::io::json::{JsonPretty, json_write_file};
use treetime_utils::make_internal_report;
use util_augur_node_data_json::{
  AugurNodeDataJsonGeneratedBy, AugurNodeDataJsonRefine, AugurNodeDataJsonRefineMeta, AugurNodeDataJsonRefineNode,
};

/// Write augur-compatible node data JSON for the `optimize` command.
///
/// Produces the structure consumed by `augur export v2 --node-data`, equivalent
/// to `augur refine` run WITHOUT `--timetree`: top-level `alignment` and
/// `input_tree` paths and a per-node `branch_length`. No clock, date, or
/// confidence-interval fields, because `optimize` performs no temporal inference.
///
/// Augur builds this JSON itself from TreeTime tree-node attributes (it never
/// calls the TreeTime CLI), so the contract is augur's `refine.py` traced against
/// the non-timetree branch (`refine.py:233` `attributes = ['branch_length',
/// 'confidence']`; the `if args.timetree:` block at `refine.py:243` is skipped).
/// `collect_node_data` (`refine.py:94`) drops `None`-valued attributes, and under
/// the default `--divergence-units=mutations-per-site` `branch_length` stays the
/// float subs/site value.
///
/// Field semantics:
///
/// - `branch_length` = the ML-optimized branch length in substitutions/site, read
///   from the parent edge (`EdgeAncestral.branch_length` after the optimization
///   loop). `augur export v2` `node_div()` (`export_v2.py:114`) consumes
///   `branch_length` for cumulative divergence when `mutation_length` is absent.
///   The root has no incoming edge, so it carries `0.0` (the field is
///   non-optional; `export v2` sets the root divergence to 0 regardless).
/// - `confidence` is omitted: v1's Newick reader does not parse input-tree branch
///   support values (shared gap with `timetree`, tracked in
///   `kb/issues/N-timetree-node-data-confidence-not-emitted.md`).
pub fn build_augur_node_data_json(
  graph: &GraphAncestral,
  alignment: Option<&Path>,
  input_tree: Option<&Path>,
) -> Result<AugurNodeDataJsonRefine, Report> {
  let mut nodes = BTreeMap::new();
  for node in graph.get_nodes() {
    let node_guard = node.read_arc();
    let node_key = node_guard.key();
    let payload = node_guard.payload().read_arc();
    let node_name = payload
      .name()
      .map_or_else(|| format!("node_{}", node_key.0), |n| n.as_ref().to_owned());

    // branch_length comes from the parent edge. The root has no incoming branch,
    // so it carries 0.0 (matching the timetree writer and `export v2`, which sets
    // the root divergence to 0).
    let branch_length = match graph.node_parent(node_key)? {
      Some((_parent_key, edge_key)) => {
        let edge = graph
          .get_edge(edge_key)
          .ok_or_else(|| make_internal_report!("Optimize node data: missing edge {edge_key:?}"))?;
        edge.read_arc().payload().read_arc().branch_length().unwrap_or(0.0)
      },
      None => 0.0,
    };

    nodes.insert(
      node_name,
      AugurNodeDataJsonRefineNode {
        branch_length,
        // optimize emits no temporal or confidence fields.
        confidence: None,
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
      // No clock: optimize performs no molecular-clock or time-tree inference.
      clock: None,
      other: BTreeMap::new(),
    },
    nodes,
  })
}

pub fn write_augur_node_data_json(
  graph: &GraphAncestral,
  alignment: Option<&Path>,
  input_tree: Option<&Path>,
  path: &Path,
) -> Result<(), Report> {
  let data = build_augur_node_data_json(graph, alignment, input_tree)?;
  json_write_file(path, &data, JsonPretty(true))?;
  Ok(())
}
