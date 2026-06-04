use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::ancestral::attach::sanitize_to_alphabet;
use crate::make_error;
use crate::partition::augur::AugurNodeDataJsonAncestralPartition;
use crate::partition::traits::BranchTopology;
use crate::payload::ancestral::GraphAncestral;
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::Itertools;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use treetime_graph::node::{GraphNodeKey, Named};
use treetime_io::fasta::read_many_fasta;
use treetime_io::gff::{GffCdsFeature, read_gff3_cds_features_filtered};
use treetime_primitives::{AsciiChar, Seq};
use util_augur_node_data_json::{AugurNodeDataJsonAnnotationEntry, AugurNodeDataJsonAnnotationSegment};

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AaNodeData {
  pub annotations: BTreeMap<String, AugurNodeDataJsonAnnotationEntry>,
  pub reference: BTreeMap<String, String>,
  // Keyed by graph node key, not node name: every amino-acid partition is reconstructed on the same
  // shared tree as the nucleotide partition, so per-node results join by node identity (key) rather
  // than by reconstructing identity from a synthesized node name across independent graphs.
  pub node_aa_muts: BTreeMap<GraphNodeKey, BTreeMap<String, Vec<String>>>,
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
  }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AaCdsNodeData {
  pub reference: String,
  pub root_sequence: String,
  pub node_muts: BTreeMap<GraphNodeKey, Vec<String>>,
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
  for node in graph.get_nodes() {
    let node_guard = node.read_arc();
    let node_key = node_guard.key();
    let payload = node_guard.payload().read_arc();
    let node_name = payload
      .name()
      .map_or_else(|| format!("node_{}", node_key.0), |n| n.as_ref().to_owned());

    let muts = if node_key == root_key {
      diff_sequences(&reference, &inferred_root, partition.ambiguous_char())?
    } else {
      let (_parent_key, edge_key) = graph
        .node_parent(node_key)?
        .ok_or_else(|| eyre::eyre!("Non-root node '{node_name}' has no parent while collecting AA node data"))?;
      partition
        .edge_subs(graph, edge_key)?
        .into_iter()
        .sorted_by_key(Sub::pos)
        .map(|sub| sub.to_string())
        .collect()
    };

    node_muts.insert(node_key, muts);
  }

  Ok(AaCdsNodeData {
    reference: reference.as_str().to_owned(),
    root_sequence: inferred_root.as_str().to_owned(),
    node_muts,
  })
}

fn diff_sequences(reference: &Seq, query: &Seq, unknown: AsciiChar) -> Result<Vec<String>, Report> {
  if reference.len() != query.len() {
    return make_error!(
      "Cannot diff sequences with lengths {} and {}",
      reference.len(),
      query.len()
    );
  }

  Ok(
    reference
      .iter()
      .zip(query.iter())
      .enumerate()
      .filter(|(_pos, (reff, qry))| reff != qry && is_reportable_sub(**reff, **qry, unknown))
      .map(|(pos, (reff, qry))| format!("{}{}{}", char::from(*reff), pos + 1, char::from(*qry)))
      .collect(),
  )
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

#[cfg(test)]
mod tests {
  use super::*;
  use crate::seq::mutation::Sub;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_utils::o;

  #[test]
  fn test_validate_aa_args_requires_cds_placeholder() {
    let err = validate_aa_args(&Some("translations.fasta".to_owned()), &["S".to_owned()], &None, &None).unwrap_err();
    assert!(err.to_string().contains("CDS placeholder"));
  }

  #[test]
  fn test_validate_aa_args_accepts_percent_gene_placeholder() {
    let result = validate_aa_args(&Some("out/%GENE.fasta".to_owned()), &["S".to_owned()], &None, &None);
    result.unwrap();
  }

  #[test]
  fn test_validate_aa_args_accepts_cds_placeholder() {
    let template = ["{", "cds", "}"].concat();
    let result = validate_aa_args(&Some(format!("out/{template}.fasta")), &["S".to_owned()], &None, &None);
    result.unwrap();
  }

  #[test]
  fn test_validate_aa_args_empty_cdses_with_annotation_ok() {
    let template = ["{", "cds", "}"].concat();
    let result = validate_aa_args(
      &Some(format!("out/{template}.fasta")),
      &[],
      &Some(PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/Cargo.toml"))),
      &None,
    );
    result.unwrap();
  }

  #[test]
  fn test_validate_aa_args_empty_cdses_no_annotation_errors() {
    let template = ["{", "cds", "}"].concat();
    let err = validate_aa_args(&Some(format!("out/{template}.fasta")), &[], &None, &None).unwrap_err();
    assert!(err.to_string().contains("--cdses"));
  }

  #[allow(clippy::literal_string_with_formatting_args)]
  #[rustfmt::skip]
  #[rstest]
  #[case::cds_placeholder( "out/{cds}.fasta",  "S",  "out/S.fasta")]
  #[case::gene_placeholder("out/%GENE.fasta",   "S",  "out/S.fasta")]
  #[allow(clippy::literal_string_with_formatting_args)]
  #[case::both_placeholders("out/{cds}/%GENE.fasta", "ORF1a", "out/ORF1a/ORF1a.fasta")]
  fn test_translation_path_expands_placeholders(
    #[case] template: &str,
    #[case] cds: &str,
    #[case] expected: &str,
  ) {
    assert_eq!(PathBuf::from(expected), translation_path(template, cds));
  }

  #[test]
  fn test_validate_aa_root_sequence_cdses_requires_every_cds() {
    let by_cds = btreemap! {
      o!("S") => Seq::try_from_str("AC").unwrap(),
    };

    let err = validate_aa_root_sequence_cdses(Path::new("roots.fasta"), &by_cds, &[o!("S"), o!("M")]).unwrap_err();

    assert!(err.to_string().contains("CDS 'M'"));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::single_span(
    AugurNodeDataJsonAnnotationEntry { start: Some(100), end: Some(400), segments: None, ..Default::default() },
    Some(301)
  )]
  #[case::segments(
    AugurNodeDataJsonAnnotationEntry {
      segments: Some(vec![
        AugurNodeDataJsonAnnotationSegment { start: 1, end: 100, other: btreemap! {} },
        AugurNodeDataJsonAnnotationSegment { start: 200, end: 300, other: btreemap! {} },
      ]),
      ..Default::default()
    },
    Some(201)
  )]
  #[case::neither(
    AugurNodeDataJsonAnnotationEntry::default(),
    None
  )]
  fn test_annotation_cds_nuc_length(
    #[case] entry: AugurNodeDataJsonAnnotationEntry,
    #[case] expected: Option<i64>,
  ) {
    assert_eq!(expected, annotation_cds_nuc_length(&entry));
  }

  #[test]
  fn test_diff_sequences_skips_gap_and_unknown_states() {
    let reference = Seq::try_from_str("ACDX-").unwrap();
    let query = Seq::try_from_str("ADQXF").unwrap();

    let actual = diff_sequences(&reference, &query, AsciiChar::from_byte_unchecked(b'X')).unwrap();

    let expected = vec![o!("C2D"), o!("D3Q")];
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_collect_aa_cds_node_data_keeps_inferred_root_sequence() {
    let graph = helpers::named_tree();
    let name_to_key = helpers::node_name_to_key(&graph);
    let partition = helpers::StubAugurPartition::new(
      &graph,
      &btreemap! {
        o!("A") => o!("AD"),
        o!("B") => o!("AC"),
        o!("root") => o!("AC"),
      },
    );
    let reference = Seq::try_from_str("AA").unwrap();

    let actual = collect_aa_cds_node_data(&graph, &partition, "S", Some(&reference)).unwrap();

    let expected = AaCdsNodeData {
      reference: o!("AA"),
      root_sequence: o!("AC"),
      node_muts: btreemap! {
        name_to_key["A"] => vec![],
        name_to_key["B"] => vec![],
        name_to_key["root"] => vec![o!("A2C")],
      },
    };
    assert_eq!(expected, actual);
  }

  mod helpers {
    use super::*;
    use treetime_graph::node::{GraphNodeKey, Named};
    use treetime_io::nwk::nwk_read_str;

    pub fn node_name_to_key(graph: &GraphAncestral) -> BTreeMap<String, GraphNodeKey> {
      graph
        .get_nodes()
        .into_iter()
        .map(|node| {
          let node = node.read_arc();
          let key = node.key();
          let name = node.payload().read_arc().name().unwrap().as_ref().to_owned();
          (name, key)
        })
        .collect()
    }

    pub fn named_tree() -> GraphAncestral {
      nwk_read_str("(A:0.1,B:0.1)root;").unwrap()
    }

    pub struct StubAugurPartition {
      sequences: BTreeMap<GraphNodeKey, Seq>,
      unknown: AsciiChar,
    }

    impl StubAugurPartition {
      pub fn new(graph: &GraphAncestral, sequences_by_name: &BTreeMap<String, String>) -> Self {
        let sequences = graph
          .get_nodes()
          .into_iter()
          .map(|node| {
            let node = node.read_arc();
            let payload = node.payload().read_arc();
            let name = payload.name().unwrap();
            (
              node.key(),
              Seq::try_from_str(&sequences_by_name[name.as_ref()]).unwrap(),
            )
          })
          .collect();
        Self {
          sequences,
          unknown: AsciiChar::from_byte_unchecked(b'X'),
        }
      }
    }

    impl AugurNodeDataJsonAncestralPartition for StubAugurPartition {
      fn sequence_length(&self) -> usize {
        self.sequences.values().next().unwrap().len()
      }

      fn node_sequence(&self, node_key: GraphNodeKey) -> Seq {
        self.sequences[&node_key].clone()
      }

      fn edge_subs(&self, _graph: &dyn BranchTopology, _edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
        Ok(vec![])
      }

      fn ambiguous_char(&self) -> AsciiChar {
        self.unknown
      }
    }
  }
}
