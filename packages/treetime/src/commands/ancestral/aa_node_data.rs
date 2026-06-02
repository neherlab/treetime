use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::ancestral::attach::sanitize_to_alphabet;
use crate::make_error;
use crate::partition::augur::AugurNodeDataJsonAncestralPartition;
use crate::partition::traits::BranchTopology;
use crate::payload::ancestral::GraphAncestral;
use crate::seq::mutation::Sub;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use percent_encoding::percent_decode_str;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use treetime_graph::node::{GraphNodeKey, Named};
use treetime_io::fasta::read_many_fasta;
use treetime_primitives::{AsciiChar, Seq};
use treetime_utils::io::fs::read_file_to_string;
use util_augur_node_data_json::{AugurNodeDataJsonAnnotationEntry, AugurNodeDataJsonAnnotationSegment};

/// CDS-name attribute priority matching nextclade's `NAME_ATTRS_CDS` in `gff3_reader.rs`.
const NAME_ATTRS_CDS: &[&str] = &[
  "Name",
  "name",
  "Alias",
  "alias",
  "standard_name",
  "old-name",
  "Gene",
  "gene",
  "gene_name",
  "locus_tag",
  "product",
  "gene_synonym",
  "gb-synonym",
  "acronym",
  "gb-acronym",
  "protein_id",
  "ID",
];

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

  let contents = read_file_to_string(path)?;
  parse_gff3_annotations(&contents, cdses, path)
}

fn parse_gff3_annotations(
  contents: &str,
  cdses: &[String],
  path: &Path,
) -> Result<BTreeMap<String, AugurNodeDataJsonAnnotationEntry>, Report> {
  let wanted: BTreeSet<&str> = cdses.iter().map(String::as_str).collect();
  let mut features: BTreeMap<String, Vec<GffCdsRow>> = BTreeMap::new();

  for (line_no, line) in contents.lines().enumerate() {
    let line = line.trim();
    if line == "##FASTA" {
      break;
    }
    if line.is_empty() || line.starts_with('#') {
      continue;
    }

    let cols = line.split('\t').collect_vec();
    if cols.len() != 9 {
      return make_error!(
        "Invalid GFF3 line {} in '{}': expected 9 tab-separated columns",
        line_no + 1,
        path.display()
      );
    }
    if cols[2] != "CDS" {
      continue;
    }

    let attrs = parse_gff_attributes(cols[8]).wrap_err_with(|| {
      format!(
        "Invalid GFF3 attributes on line {} in '{}'",
        line_no + 1,
        path.display()
      )
    })?;
    let Some(cds) = resolve_cds_name(&attrs) else {
      continue;
    };

    // Empty `wanted` means derive the CDS set from the annotation itself (every CDS feature).
    if !wanted.is_empty() && !wanted.contains(cds.as_str()) {
      continue;
    }

    let start = cols[3].parse::<i64>()?;
    let end = cols[4].parse::<i64>()?;
    if start < 1 || end < start {
      return make_error!("Invalid CDS coordinates for CDS '{cds}' on GFF3 line {}", line_no + 1);
    }

    let strand = match cols[6] {
      "+" | "-" => cols[6].to_owned(),
      strand => {
        return make_error!(
          "Invalid CDS strand '{strand}' for CDS '{cds}' on GFF3 line {}",
          line_no + 1
        );
      },
    };

    features.entry(cds).or_default().push(GffCdsRow {
      seqid: cols[0].to_owned(),
      start,
      end,
      strand,
    });
  }

  for cds in cdses {
    if !features.contains_key(cds) {
      return make_error!(
        "--annotation '{}' does not contain a CDS feature for CDS '{}'",
        path.display(),
        cds
      );
    }
  }

  features
    .into_iter()
    .map(|(cds, mut rows)| {
      rows.sort_by_key(|row| (row.start, row.end));
      let strand = rows.first().expect("features are non-empty").strand.clone();
      let seqid = rows.first().expect("features are non-empty").seqid.clone();

      if rows.iter().any(|row| row.strand != strand) {
        return make_error!("CDS feature rows for CDS '{cds}' use multiple strands");
      }

      if rows.iter().any(|row| row.seqid != seqid) {
        return make_error!("CDS feature rows for CDS '{cds}' span multiple sequence IDs");
      }

      // Auspice requires segments in 5'-to-3' order. GFF coordinates are always ascending, so the
      // 5' end is the smallest coordinate on the plus strand and the largest on the minus strand.
      // After the ascending sort above, reverse for the minus strand to get 5'-to-3' order. This
      // matches BioPython's `CompoundLocation.parts` ordering that augur relies on.
      if strand == "-" {
        rows.reverse();
      }

      // The CDS nucleotide length must be a multiple of 3 to translate to whole codons. Auspice
      // silently drops a non-divisible CDS from the genome map; augur errors. We error.
      let cds_length: i64 = rows.iter().map(|row| row.end - row.start + 1).sum();
      if cds_length % 3 != 0 {
        return make_error!(
          "CDS feature for CDS '{cds}' has nucleotide length {cds_length}, which is not a multiple of 3. \
           Check that the annotation matches the translations."
        );
      }

      let mut other = BTreeMap::new();
      other.insert("seqid".to_owned(), json!(seqid));

      let annotation = if rows.len() == 1 {
        let row = rows.first().expect("features are non-empty");
        AugurNodeDataJsonAnnotationEntry {
          start: Some(row.start),
          end: Some(row.end),
          strand: Some(strand),
          entry_type: Some("CDS".to_owned()),
          segments: None,
          other,
        }
      } else {
        AugurNodeDataJsonAnnotationEntry {
          start: None,
          end: None,
          strand: Some(strand),
          entry_type: Some("CDS".to_owned()),
          segments: Some(
            rows
              .iter()
              .map(|row| AugurNodeDataJsonAnnotationSegment {
                start: row.start,
                end: row.end,
                other: BTreeMap::new(),
              })
              .collect(),
          ),
          other,
        }
      };

      Ok((cds, annotation))
    })
    .collect()
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

fn resolve_cds_name(attrs: &BTreeMap<String, String>) -> Option<String> {
  NAME_ATTRS_CDS.iter().find_map(|key| attrs.get(*key).cloned())
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

fn parse_gff_attributes(raw: &str) -> Result<BTreeMap<String, String>, Report> {
  raw
    .split(';')
    .filter(|attr| !attr.is_empty())
    .map(|attr| {
      let Some((key, value)) = attr.split_once('=') else {
        return make_error!("Invalid GFF3 attribute '{attr}': expected key=value");
      };
      Ok((decode_gff3_attribute(key)?, decode_gff3_attribute(value)?))
    })
    .collect()
}

fn decode_gff3_attribute(raw: &str) -> Result<String, Report> {
  Ok(percent_decode_str(raw).decode_utf8()?.into_owned())
}

#[derive(Clone, Debug)]
struct GffCdsRow {
  seqid: String,
  start: i64,
  end: i64,
  strand: String,
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

  #[test]
  fn test_parse_gff3_annotations_returns_exact_annotations() {
    let input = "\
##gff-version 3
MN908947.3\tNextclade\tCDS\t21563\t25384\t.\t+\t0\tID=cds-S;Name=S
MN908947.3\tNextclade\tCDS\t266\t13468\t.\t+\t0\tID=cds-ORF1a;Name=ORF%31a
MN908947.3\tNextclade\tCDS\t13468\t21555\t.\t+\t0\tID=cds-ORF1a-2;Name=ORF%31a
MN908947.3\tNextclade\tgene\t100\t200\t.\t+\t.\tName=ignored
";
    let cdses = vec![o!("S"), o!("ORF1a")];

    let actual = parse_gff3_annotations(input, &cdses, Path::new("annotations.gff3")).unwrap();

    let expected = btreemap! {
      o!("ORF1a") => AugurNodeDataJsonAnnotationEntry {
        start: None,
        end: None,
        strand: Some(o!("+")),
        entry_type: Some(o!("CDS")),
        segments: Some(vec![
          AugurNodeDataJsonAnnotationSegment {
            start: 266,
            end: 13468,
            other: btreemap! {},
          },
          AugurNodeDataJsonAnnotationSegment {
            start: 13468,
            end: 21555,
            other: btreemap! {},
          },
        ]),
        other: btreemap! {
          o!("seqid") => json!("MN908947.3"),
        },
      },
      o!("S") => AugurNodeDataJsonAnnotationEntry {
        start: Some(21563),
        end: Some(25384),
        strand: Some(o!("+")),
        entry_type: Some(o!("CDS")),
        segments: None,
        other: btreemap! {
          o!("seqid") => json!("MN908947.3"),
        },
      },
    };
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_parse_gff3_annotations_orders_minus_strand_segments_5_to_3() {
    let input = "\
##gff-version 3
seq1\tNextclade\tCDS\t100\t150\t.\t-\t0\tName=ORF
seq1\tNextclade\tCDS\t200\t250\t.\t-\t0\tName=ORF
";
    let cdses = vec![o!("ORF")];

    let actual = parse_gff3_annotations(input, &cdses, Path::new("annotations.gff3")).unwrap();

    let expected = btreemap! {
      o!("ORF") => AugurNodeDataJsonAnnotationEntry {
        start: None,
        end: None,
        strand: Some(o!("-")),
        entry_type: Some(o!("CDS")),
        segments: Some(vec![
          AugurNodeDataJsonAnnotationSegment { start: 200, end: 250, other: btreemap! {} },
          AugurNodeDataJsonAnnotationSegment { start: 100, end: 150, other: btreemap! {} },
        ]),
        other: btreemap! { o!("seqid") => json!("seq1") },
      },
    };
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_parse_gff3_annotations_rejects_non_multiple_of_three() {
    let input = "\
##gff-version 3
seq1\tNextclade\tCDS\t1\t5\t.\t+\t0\tName=BAD
";
    let cdses = vec![o!("BAD")];

    let err = parse_gff3_annotations(input, &cdses, Path::new("annotations.gff3")).unwrap_err();

    assert!(err.to_string().contains("not a multiple of 3"));
  }

  #[test]
  fn test_parse_gff3_annotations_derives_all_cds_when_cdses_empty() {
    let input = "\
##gff-version 3
seq1\tNextclade\tCDS\t1\t6\t.\t+\t0\tName=A
seq1\tNextclade\tCDS\t10\t18\t.\t+\t0\tName=B
";

    let actual = parse_gff3_annotations(input, &[], Path::new("annotations.gff3")).unwrap();

    assert_eq!(2, actual.len());
    assert_eq!(Some(1), actual["A"].start);
    assert_eq!(Some(6), actual["A"].end);
    assert_eq!(Some(10), actual["B"].start);
    assert_eq!(Some(18), actual["B"].end);
  }

  #[test]
  fn test_parse_gff3_annotations_stops_at_fasta_directive() {
    let input = "\
##gff-version 3
seq1\tNextclade\tCDS\t1\t6\t.\t+\t0\tName=A
##FASTA
>seq1
ACGTAC
";

    let actual = parse_gff3_annotations(input, &[], Path::new("annotations.gff3")).unwrap();
    assert_eq!(1, actual.len());
    assert!(actual.contains_key("A"));
  }

  #[test]
  fn test_parse_gff3_annotations_rejects_mixed_seqids() {
    let input = "\
##gff-version 3
chr1\tNextclade\tCDS\t1\t3\t.\t+\t0\tName=X
chr2\tNextclade\tCDS\t4\t6\t.\t+\t0\tName=X
";
    let err = parse_gff3_annotations(input, &[o!("X")], Path::new("annotations.gff3")).unwrap_err();
    assert!(err.to_string().contains("multiple sequence IDs"));
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

  #[rustfmt::skip]
  #[rstest]
  #[case::name(        btreemap! { o!("Name") => o!("S") },                          Some(o!("S")))]
  #[case::locus_tag(   btreemap! { o!("locus_tag") => o!("LT1"), o!("Parent") => o!("p") }, Some(o!("LT1")))]
  #[case::gene_over_id(btreemap! { o!("gene") => o!("G"), o!("ID") => o!("cds-G") },  Some(o!("G")))]
  #[case::parent_only( btreemap! { o!("Parent") => o!("gene-Y") },                    None)]
  #[case::id_fallback( btreemap! { o!("ID") => o!("cds-X") },                         Some(o!("cds-X")))]
  fn test_resolve_cds_name_priority(#[case] attrs: BTreeMap<String, String>, #[case] expected: Option<String>) {
    assert_eq!(expected, resolve_cds_name(&attrs));
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
