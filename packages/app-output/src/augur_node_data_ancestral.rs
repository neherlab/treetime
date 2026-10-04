use crate::ancestral_result::{AncestralOutputMaps, AugurOutputMaps};
use eyre::Report;
use itertools::Itertools;
use maplit::btreemap;
use std::collections::BTreeMap;
use std::path::Path;
use treetime::ancestral::aa::AaNodeData;
use treetime::ancestral::mask::mask_to_string;
use treetime::seq::mutation::{MutationEvent, mutation_event_strings};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::error::make_internal_error;
use treetime_utils::io::json::{JsonPretty, json_write_file};
use util_augur_node_data_json::{
  AugurNodeDataJsonAncestral, AugurNodeDataJsonAncestralMeta, AugurNodeDataJsonAncestralNode,
  AugurNodeDataJsonAnnotationEntry, AugurNodeDataJsonAnnotations, AugurNodeDataJsonGeneratedBy,
};

pub fn write_augur_node_data_json_with_aa(
  graph: &Graph,
  mutations: &AncestralOutputMaps,
  sequences: AugurOutputMaps,
  mask: &[bool],
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  aa_node_data: Option<&AaNodeData>,
  aa_annotations: &BTreeMap<String, AugurNodeDataJsonAnnotationEntry>,
  path: &Path,
) -> Result<(), Report> {
  let data = build_augur_node_data_json(graph, mutations, sequences, mask, names, aa_node_data, aa_annotations)?;
  json_write_file(path, &data, JsonPretty(true))?;
  Ok(())
}

pub fn build_augur_node_data_json(
  graph: &Graph,
  mutations: &AncestralOutputMaps,
  sequences: AugurOutputMaps,
  mask: &[bool],
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  aa_node_data: Option<&AaNodeData>,
  aa_annotations: &BTreeMap<String, AugurNodeDataJsonAnnotationEntry>,
) -> Result<AugurNodeDataJsonAncestral, Report> {
  let AugurOutputMaps {
    mut node_sequences,
    sequence_length: alignment_length,
    ambiguous_char: ambiguous,
  } = sequences;
  let reference_seq = &mutations.root_sequence;
  let masked_positions = mask.iter().positions(|&masked| masked).collect_vec();

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
  annotations.other.extend(aa_annotations.clone());

  let mut nodes = BTreeMap::new();
  let root_key = graph.root_key()?;
  for node in graph.get_nodes() {
    let node_guard = node;
    let node_key = node_guard.key();
    let node_name = names[&node_key]
      .as_deref()
      .map_or_else(|| format!("node_{}", node_key.0), str::to_owned);

    let muts = match graph.node_parent(node_key)? {
      Some((_parent_key, edge_key)) => mutations.edge_mutations[&edge_key]
        .iter()
        .filter_map(|mutation| match &mutation.event {
          MutationEvent::Substitution(sub) => Some(sub),
          MutationEvent::Insertion(_) | MutationEvent::Deletion(_) => None,
        })
        .filter(|sub| !mask[sub.pos()])
        .sorted_by_key(|sub| sub.pos())
        .map(|sub| sub.to_string())
        .collect(),
      None => Vec::new(),
    };

    let Some(mut sequence) = node_sequences.remove(&node_key) else {
      return make_internal_error!("Augur node data: no sequence for node {}", node_key.0);
    };
    let length = sequence.len();
    for &pos in masked_positions.iter().take_while(|&&pos| pos < length) {
      sequence[pos] = ambiguous;
    }

    let aa_muts = aa_node_data
      .and_then(|aa| aa.node_aa_mutations.get(&node_key))
      .map(aa_mutation_strings)
      .transpose()?;
    nodes.insert(
      node_name,
      AugurNodeDataJsonAncestralNode {
        muts,
        sequence: Some(sequence.into_string()),
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

fn aa_mutation_strings(tracks: &BTreeMap<String, Vec<MutationEvent>>) -> Result<BTreeMap<String, Vec<String>>, Report> {
  tracks
    .iter()
    .map(|(cds, events)| {
      let strings: Vec<String> = events.iter().map(mutation_event_strings).flatten_ok().try_collect()?;
      Ok((cds.clone(), strings))
    })
    .collect()
}
