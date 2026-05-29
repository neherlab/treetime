use crate::ancestral::mask::mask_to_string;
use crate::partition::augur::AugurNodeDataJsonAncestralPartition;
use crate::partition::traits::BranchTopology;
use crate::payload::ancestral::GraphAncestral;
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::Itertools;
use maplit::btreemap;
use std::collections::BTreeMap;
use std::path::Path;
use treetime_graph::node::Named;
use treetime_utils::io::json::{JsonPretty, json_write_file};
use util_augur_node_data_json::{
  AugurNodeDataJsonAncestral, AugurNodeDataJsonAncestralMeta, AugurNodeDataJsonAncestralNode,
  AugurNodeDataJsonAnnotationEntry, AugurNodeDataJsonAnnotations, AugurNodeDataJsonGeneratedBy,
};

/// Write augur-compatible node data JSON for the `ancestral` command.
///
/// Produces the nuc-only structure consumed by `augur export v2 --node-data`:
/// top-level `annotations.nuc`, `reference.nuc` (unmasked reconstructed root),
/// the per-position `mask` string, and per-node `muts` (mask-filtered) and
/// `sequence` (masked positions set to the ambiguous character, matching augur's
/// `collect_sequences`).
pub fn write_augur_node_data_json(
  graph: &GraphAncestral,
  partition: &dyn AugurNodeDataJsonAncestralPartition,
  mask: &[bool],
  path: &Path,
) -> Result<(), Report> {
  let alignment_length = partition.sequence_length();
  let reference_seq = partition.root_sequence(graph)?;
  let ambiguous = partition.ambiguous_char();

  let annotations = AugurNodeDataJsonAnnotations {
    nuc: Some(AugurNodeDataJsonAnnotationEntry {
      start: Some(1),
      end: Some(i64::try_from(alignment_length)?),
      strand: Some("+".to_owned()),
      entry_type: Some("source".to_owned()),
      segments: None,
      other: BTreeMap::new(),
    }),
    other: BTreeMap::new(),
  };

  let mut nodes = BTreeMap::new();
  for node in graph.get_nodes() {
    let node_guard = node.read_arc();
    let node_key = node_guard.key();
    let payload = node_guard.payload().read_arc();
    let node_name = payload
      .name()
      .map_or_else(|| format!("node_{}", node_key.0), |n| n.as_ref().to_owned());

    // Root has no parent edge, so it carries no mutations (augur emits []).
    let muts = match graph.node_parent(node_key)? {
      Some((_parent_key, edge_key)) => partition
        .edge_subs(graph, edge_key)?
        .into_iter()
        .filter(|sub| !mask[sub.pos()])
        .sorted_by_key(Sub::pos)
        .map(|sub| sub.to_string())
        .collect(),
      None => Vec::new(),
    };

    // Masked positions (every tip ambiguous) are reported as the ambiguous
    // character in the output sequence, matching augur.
    let mut sequence = partition.node_sequence(node_key);
    for (pos, &masked) in mask.iter().enumerate() {
      if masked && pos < sequence.len() {
        sequence[pos] = ambiguous;
      }
    }

    nodes.insert(
      node_name,
      AugurNodeDataJsonAncestralNode {
        muts,
        sequence: Some(sequence.as_str().to_owned()),
        aa_muts: None,
        aa_sequences: None,
        other: BTreeMap::new(),
      },
    );
  }

  let data = AugurNodeDataJsonAncestral {
    generated_by: Some(AugurNodeDataJsonGeneratedBy {
      program: "treetime".to_owned(),
      version: env!("CARGO_PKG_VERSION").to_owned(),
    }),
    metadata: AugurNodeDataJsonAncestralMeta {
      annotations: Some(annotations),
      reference: Some(btreemap! { "nuc".to_owned() => reference_seq.as_str().to_owned() }),
      mask: Some(mask_to_string(mask)),
      other: BTreeMap::new(),
    },
    nodes,
  };

  json_write_file(path, &data, JsonPretty(true))?;
  Ok(())
}
