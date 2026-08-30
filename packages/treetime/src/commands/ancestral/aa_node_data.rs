#[cfg(test)]
mod __tests__;

use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::ancestral::attach::sanitize_to_alphabet;
use crate::make_error;
use crate::partition::io::augur::AugurNodeDataJsonAncestralPartition;
use crate::partition::traits::BranchTopology;
use crate::payload::ancestral::GraphAncestral;
use crate::seq::mutation::{Mutation, MutationEvent, MutationTrack, Sub};
use eyre::Report;
use itertools::Itertools;
use serde::Serialize;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use treetime_graph::node::{GraphNodeKey, Named};
use treetime_io::fasta::read_many_fasta;
use treetime_io::gff::{GffCdsFeature, read_gff3_cds_features_filtered};
use treetime_primitives::{AsciiChar, Seq};
use util_augur_node_data_json::{AugurNodeDataJsonAnnotationEntry, AugurNodeDataJsonAnnotationSegment};

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize)]
pub struct AaNodeData {
  pub annotations: BTreeMap<String, AugurNodeDataJsonAnnotationEntry>,
  pub reference: BTreeMap<String, String>,
  // Keyed by graph node key, not node name: every amino-acid partition is reconstructed on the same
  // shared tree as the nucleotide partition, so per-node results join by node identity (key) rather
  // than by reconstructing identity from a synthesized node name across independent graphs.
  pub node_aa_muts: BTreeMap<GraphNodeKey, BTreeMap<String, Vec<String>>>,
  pub node_aa_mutations: BTreeMap<GraphNodeKey, BTreeMap<String, Vec<MutationEvent>>>,
  pub root_aa_sequences: BTreeMap<String, String>,
}

impl AaNodeData {
  pub fn add_cds(&mut self, cds: &str, cds_data: AaCdsNodeData, annotation: Option<AugurNodeDataJsonAnnotationEntry>) {
    if let Some(annotation) = annotation {
      self.annotations.insert(cds.to_owned(), annotation);
    }
    self.reference.insert(cds.to_owned(), cds_data.reference);
    self.root_aa_sequences.insert(cds.to_owned(), cds_data.root_sequence);
    for (node_key, muts) in cds_data.node_muts {
      self
        .node_aa_muts
        .entry(node_key)
        .or_default()
        .insert(cds.to_owned(), muts);
    }
    for (node_key, mutations) in cds_data.node_mutations {
      self
        .node_aa_mutations
        .entry(node_key)
        .or_default()
        .insert(cds.to_owned(), mutations);
    }
  }
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize)]
pub struct AaCdsNodeData {
  pub reference: String,
  pub root_sequence: String,
  pub node_muts: BTreeMap<GraphNodeKey, Vec<String>>,
  pub node_mutations: BTreeMap<GraphNodeKey, Vec<MutationEvent>>,
}

pub fn validate_aa_args(
  translations: &Option<String>,
  cdses: &[String],
  annotation: &Option<PathBuf>,
  aa_root_sequence: &Option<PathBuf>,
) -> Result<(), Report> {
  if translations.is_none() && cdses.is_empty() && annotation.is_none() && aa_root_sequence.is_none() {
    return Ok(());
  }

  let Some(template) = translations else {
    return make_error!("--translations is required when using --cdses, --annotation, or --aa-root-sequence");
  };

  #[allow(clippy::literal_string_with_formatting_args)]
  if !template_has_cds_placeholder(template) {
    return make_error!("--translations must contain a CDS placeholder ('{{cds}}' or '%GENE')");
  }

  // The CDS set is either listed explicitly or derived from the annotation; require at least one source.
  if cdses.is_empty() && annotation.is_none() {
    return make_error!("--cdses must list at least one CDS, or pass --annotation to derive the CDS set");
  }

  validate_file_arg("--annotation", annotation.as_deref())?;
  validate_file_arg("--aa-root-sequence", aa_root_sequence.as_deref())?;

  let mut seen = BTreeSet::new();
  for cds in cdses {
    if cds.is_empty() {
      return make_error!("--cdses must not contain an empty CDS name");
    }
    if !seen.insert(cds) {
      return make_error!("--cdses contains duplicate CDS '{cds}'");
    }
  }

  Ok(())
}

fn validate_file_arg(arg_name: &str, path: Option<&Path>) -> Result<(), Report> {
  if let Some(path) = path
    && !path.is_file()
  {
    return make_error!("{arg_name} '{}' does not exist or is not a file", path.display());
  }
  Ok(())
}

/// CDS placeholders accepted in path templates: `{cds}` (Nextclade) and `%GENE` (augur).
const CDS_PLACEHOLDERS: &[&str] = &["{cds}", "%GENE"];

pub fn template_has_cds_placeholder(template: &str) -> bool {
  CDS_PLACEHOLDERS
    .iter()
    .any(|placeholder| template.contains(placeholder))
}

pub fn translation_path(template: &str, cds: &str) -> PathBuf {
  let mut path = template.to_owned();
  for placeholder in CDS_PLACEHOLDERS {
    path = path.replace(placeholder, cds);
  }
  PathBuf::from(path)
}

pub fn read_aa_root_sequences(
  path: Option<&Path>,
  cdses: &[String],
  recon_alphabet: &Alphabet,
) -> Result<BTreeMap<String, Seq>, Report> {
  let Some(path) = path else {
    return Ok(BTreeMap::new());
  };

  // Read with the stop-inclusive alphabet, then fold out-of-alphabet characters into the unknown
  // state of the reconstruction alphabet so the root sequence shares its alphabet with the partition.
  let read_alphabet = Alphabet::new(AlphabetName::Aa)?;
  let records = read_many_fasta(&[path], &read_alphabet)?;
  let mut by_cds = BTreeMap::new();
  for record in records {
    let (seq, _changed) = sanitize_to_alphabet(&record.seq, recon_alphabet);
    by_cds.insert(record.seq_name, seq);
  }

  validate_aa_root_sequence_cdses(path, &by_cds, cdses)?;

  Ok(by_cds)
}

fn validate_aa_root_sequence_cdses(
  path: &Path,
  by_cds: &BTreeMap<String, Seq>,
  cdses: &[String],
) -> Result<(), Report> {
  for cds in cdses {
    if !by_cds.contains_key(cds) {
      return make_error!(
        "--aa-root-sequence '{}' does not contain a FASTA record for CDS '{}'",
        path.display(),
        cds
      );
    }
  }
  Ok(())
}

pub fn read_gff3_annotations(
  path: Option<&Path>,
  cdses: &[String],
) -> Result<BTreeMap<String, AugurNodeDataJsonAnnotationEntry>, Report> {
  let Some(path) = path else {
    return Ok(BTreeMap::new());
  };

  let features = read_gff3_cds_features_filtered(path, cdses)?;
  Ok(
    features
      .into_iter()
      .map(|feature| {
        let annotation = gff_cds_to_annotation(&feature);
        (feature.name, annotation)
      })
      .collect(),
  )
}

fn gff_cds_to_annotation(feature: &GffCdsFeature) -> AugurNodeDataJsonAnnotationEntry {
  let mut other = BTreeMap::new();
  other.insert("seqid".to_owned(), json!(feature.seqid));

  if feature.segments.len() == 1 {
    let seg = &feature.segments[0];
    AugurNodeDataJsonAnnotationEntry {
      start: Some(seg.start),
      end: Some(seg.end),
      strand: Some(feature.strand.clone()),
      entry_type: Some("CDS".to_owned()),
      segments: None,
      other,
    }
  } else {
    AugurNodeDataJsonAnnotationEntry {
      start: None,
      end: None,
      strand: Some(feature.strand.clone()),
      entry_type: Some("CDS".to_owned()),
      segments: Some(
        feature
          .segments
          .iter()
          .map(|seg| AugurNodeDataJsonAnnotationSegment {
            start: seg.start,
            end: seg.end,
            other: BTreeMap::new(),
          })
          .collect(),
      ),
      other,
    }
  }
}

pub fn collect_aa_cds_node_data(
  graph: &GraphAncestral,
  partition: &dyn AugurNodeDataJsonAncestralPartition,
  cds: &str,
  reference_override: Option<&Seq>,
) -> Result<AaCdsNodeData, Report> {
  let root_key = graph.root_key()?;
  let inferred_root = partition.root_sequence(graph)?;
  let reference = reference_override.cloned().unwrap_or_else(|| inferred_root.clone());

  if reference.len() != inferred_root.len() {
    return make_error!(
      "AA root/reference sequence for CDS '{cds}' has length {}, but inferred root has length {}",
      reference.len(),
      inferred_root.len()
    );
  }

  let mut node_muts = BTreeMap::new();
  let mut node_mutations = BTreeMap::new();
  for node in graph.get_nodes() {
    let node_guard = node.read_arc();
    let node_key = node_guard.key();
    let payload = node_guard.payload().read_arc();
    let node_name = payload
      .name()
      .map_or_else(|| format!("node_{}", node_key.0), |n| n.as_ref().to_owned());

    let mutations = if node_key == root_key {
      diff_sequences(&reference, &inferred_root, partition.ambiguous_char())?
        .into_iter()
        .map(MutationEvent::Substitution)
        .collect()
    } else {
      let (_parent_key, edge_key) = graph
        .node_parent(node_key)?
        .ok_or_else(|| eyre::eyre!("Non-root node '{node_name}' has no parent while collecting AA node data"))?;
      let substitutions = partition
        .edge_subs(graph, edge_key)?
        .into_iter()
        .sorted_by_key(Sub::pos)
        .map(MutationEvent::Substitution)
        .map(Ok);
      let indels = partition
        .edge_indels(edge_key)
        .into_iter()
        .map(|indel| Mutation::indel(MutationTrack::AminoAcid(cds.to_owned()), &indel).map(|mutation| mutation.event));
      substitutions.chain(indels).collect::<Result<Vec<_>, Report>>()?
    };
    let muts = mutations.iter().flat_map(mutation_event_strings).collect();

    node_muts.insert(node_key, muts);
    node_mutations.insert(node_key, mutations);
  }

  Ok(AaCdsNodeData {
    reference: reference.as_str().to_owned(),
    root_sequence: inferred_root.as_str().to_owned(),
    node_muts,
    node_mutations,
  })
}

fn diff_sequences(reference: &Seq, query: &Seq, unknown: AsciiChar) -> Result<Vec<Sub>, Report> {
  if reference.len() != query.len() {
    return make_error!(
      "Cannot diff sequences with lengths {} and {}",
      reference.len(),
      query.len()
    );
  }

  reference
    .iter()
    .zip(query.iter())
    .enumerate()
    .filter(|(_pos, (reff, qry))| reff != qry && is_reportable_sub(**reff, **qry, unknown))
    .map(|(pos, (reff, qry))| Sub::new(*reff, pos, *qry))
    .collect()
}

fn mutation_event_strings(event: &MutationEvent) -> Vec<String> {
  match event {
    MutationEvent::Substitution(substitution) => vec![substitution.to_string()],
    MutationEvent::Insertion(segment) => segment
      .sequence
      .iter()
      .enumerate()
      .map(|(offset, state)| format!("-{}{state}", segment.range.0 + offset + 1))
      .collect(),
    MutationEvent::Deletion(segment) => segment
      .sequence
      .iter()
      .enumerate()
      .map(|(offset, state)| format!("{state}{}-", segment.range.0 + offset + 1))
      .collect(),
  }
}

fn is_reportable_sub(reff: AsciiChar, qry: AsciiChar, unknown: AsciiChar) -> bool {
  let gap = AsciiChar::from_byte_unchecked(b'-');
  reff != gap && qry != gap && reff != unknown && qry != unknown
}

/// Total nucleotide length of a CDS annotation: the sum of its segment lengths, or the single
/// `start..=end` span, in 1-based inclusive coordinates. `None` when the entry carries neither.
pub fn annotation_cds_nuc_length(entry: &AugurNodeDataJsonAnnotationEntry) -> Option<i64> {
  if let Some(segments) = &entry.segments {
    Some(segments.iter().map(|segment| segment.end - segment.start + 1).sum())
  } else if let (Some(start), Some(end)) = (entry.start, entry.end) {
    Some(end - start + 1)
  } else {
    None
  }
}
