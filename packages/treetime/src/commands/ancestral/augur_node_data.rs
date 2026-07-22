use crate::ancestral::mask::mask_to_string;
use crate::commands::ancestral::aa_node_data::AaNodeData;
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
pub fn write_augur_node_data_json<D: Send + Sync>(
  graph: &GraphAncestral<D>,
  partition: &dyn AugurNodeDataJsonAncestralPartition,
  mask: &[bool],
  path: &Path,
) -> Result<(), Report> {
  write_augur_node_data_json_with_aa(graph, partition, mask, None, path)
}

pub fn build_augur_node_data_json<D: Send + Sync>(
  graph: &GraphAncestral<D>,
  partition: &dyn AugurNodeDataJsonAncestralPartition,
  mask: &[bool],
  aa_node_data: Option<&AaNodeData>,
) -> Result<AugurNodeDataJsonAncestral, Report> {
  let alignment_length = partition.sequence_length();
  let reference_seq = partition.root_sequence(graph)?;
  let ambiguous = partition.ambiguous_char();

  let mut annotations = AugurNodeDataJsonAnnotations {
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
  if let Some(aa_node_data) = aa_node_data {
    annotations.other.extend(aa_node_data.annotations.clone());
  }

  let mut nodes = BTreeMap::new();
  let root_key = graph.root_key()?;
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

    let aa_muts = aa_node_data.and_then(|aa| aa.node_aa_muts.get(&node_key).cloned());
    nodes.insert(
      node_name,
      AugurNodeDataJsonAncestralNode {
        muts,
        sequence: Some(sequence.as_str().to_owned()),
        aa_muts,
        aa_sequences: if node_key == root_key {
          aa_node_data
            .map(|aa| aa.root_aa_sequences.clone())
            .filter(|seqs| !seqs.is_empty())
        } else {
          None
        },
        other: BTreeMap::new(),
      },
    );
  }

  let mut reference = btreemap! { "nuc".to_owned() => reference_seq.as_str().to_owned() };
  if let Some(aa_node_data) = aa_node_data {
    reference.extend(aa_node_data.reference.clone());
  }

  Ok(AugurNodeDataJsonAncestral {
    generated_by: Some(AugurNodeDataJsonGeneratedBy {
      program: "treetime".to_owned(),
      version: env!("CARGO_PKG_VERSION").to_owned(),
    }),
    metadata: AugurNodeDataJsonAncestralMeta {
      annotations: Some(annotations),
      reference: Some(reference),
      mask: Some(mask_to_string(mask)),
      other: BTreeMap::new(),
    },
    nodes,
  })
}

pub fn write_augur_node_data_json_with_aa<D: Send + Sync>(
  graph: &GraphAncestral<D>,
  partition: &dyn AugurNodeDataJsonAncestralPartition,
  mask: &[bool],
  aa_node_data: Option<&AaNodeData>,
  path: &Path,
) -> Result<(), Report> {
  let data = build_augur_node_data_json(graph, partition, mask, aa_node_data)?;
  json_write_file(path, &data, JsonPretty(true))?;
  Ok(())
}
