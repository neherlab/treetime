use crate::annotated_graph::{AnnotatedTreeView, TreeSequences};
use eyre::Report;
use itertools::Itertools;
use maplit::btreemap;
use std::collections::BTreeMap;
use std::path::Path;
use treetime::ancestral::mask::mask_to_string;
use treetime::seq::mutation::{MutationEvent, mutation_event_strings};
use treetime_graph::assign_node_names::node_name_or_key;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{AsciiChar, Seq};
use treetime_utils::io::json::{JsonPretty, json_write_file};
use treetime_utils::{make_internal_error, make_internal_report};
use util_augur_node_data_json::{
  AugurNodeDataJsonAncestral, AugurNodeDataJsonAncestralMeta, AugurNodeDataJsonAncestralNode,
  AugurNodeDataJsonAnnotationEntry, AugurNodeDataJsonAnnotations, AugurNodeDataJsonGeneratedBy,
};

pub fn write_augur_node_data_ancestral(
  tree: &AnnotatedTreeView<'_>,
  sequences: AncestralNodeSequences<'_>,
  path: &Path,
) -> Result<(), Report> {
  json_write_file(
    path,
    &build_augur_node_data_ancestral(tree, sequences)?,
    JsonPretty(true),
  )
}

pub fn build_augur_node_data_ancestral(
  tree: &AnnotatedTreeView<'_>,
  sequences: AncestralNodeSequences<'_>,
) -> Result<AugurNodeDataJsonAncestral, Report> {
  let graph = tree.graph();
  let facts = graph
    .sequences
    .as_ref()
    .ok_or_else(|| make_internal_report!("Ancestral node data requires reconstructed sequences"))?;
  let AncestralNodeSequences {
    mut node_sequences,
    alignment_length,
    ambiguous_char,
    mask,
  } = sequences;
  let masked_positions = mask.iter().positions(|&masked| masked).collect_vec();
  let amino_acids = facts.amino_acids.as_ref();

  let mut nodes = BTreeMap::new();
  for node in graph.graph.get_nodes() {
    let key = node.key();
    let Some(mut sequence) = node_sequences.remove(&key) else {
      return make_internal_error!("Augur node data: no sequence for node {}", key.0);
    };
    let length = sequence.len();
    for &pos in masked_positions.iter().take_while(|&&pos| pos < length) {
      sequence[pos] = ambiguous_char;
    }
    let aa_muts = amino_acids
      .and_then(|amino_acids| amino_acids.node_data.node_aa_mutations.get(&key))
      .map(aa_mutation_strings)
      .transpose()?;
    let aa_sequences = (key == tree.tree().root())
      .then(|| amino_acids.map(|amino_acids| amino_acids.node_data.root_aa_sequences.clone()))
      .flatten()
      .filter(|sequences| !sequences.is_empty());
    nodes.insert(
      node_name_or_key(key, graph.names[&key].as_deref()),
      AugurNodeDataJsonAncestralNode {
        muts: nucleotide_substitutions(tree, facts, mask, key),
        sequence: Some(sequence.into_string()),
        aa_muts,
        aa_sequences,
        other: BTreeMap::new(),
      },
    );
  }

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
  let mut reference = btreemap! { "nuc".to_owned() => facts.root_sequence.as_str().to_owned() };
  if let Some(amino_acids) = amino_acids {
    annotations.other.extend(amino_acids.cdses.clone());
    reference.extend(amino_acids.node_data.reference.clone());
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

pub struct AncestralNodeSequences<'a> {
  pub node_sequences: BTreeMap<GraphNodeKey, Seq>,
  pub alignment_length: usize,
  pub ambiguous_char: AsciiChar,
  pub mask: &'a [bool],
}

fn nucleotide_substitutions(
  tree: &AnnotatedTreeView<'_>,
  facts: &TreeSequences<'_>,
  mask: &[bool],
  key: GraphNodeKey,
) -> Vec<String> {
  let Some((_, edge_key)) = tree.tree().parent(key) else {
    return vec![];
  };
  facts.edge_mutations[&edge_key]
    .iter()
    .filter_map(|mutation| match &mutation.event {
      MutationEvent::Substitution(sub) => Some(sub),
      MutationEvent::Insertion(_) | MutationEvent::Deletion(_) => None,
    })
    .filter(|sub| !mask[sub.pos()])
    .sorted_by_key(|sub| sub.pos())
    .map(|sub| sub.to_string())
    .collect()
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
