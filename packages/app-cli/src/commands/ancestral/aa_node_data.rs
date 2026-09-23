#[cfg(test)]
mod __tests__;

use eyre::Report;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use treetime::alphabet::alphabet::{Alphabet, AlphabetName};
use treetime::ancestral::attach::sanitize_to_alphabet;
use treetime::make_error;
use treetime_io::fasta::read_many_fasta_path;
use treetime_io::gff::{GffCdsFeature, read_gff3_cds_features_filtered};
use treetime_primitives::Seq;
use util_augur_node_data_json::{AugurNodeDataJsonAnnotationEntry, AugurNodeDataJsonAnnotationSegment};

pub fn validate_aa_args(
  translations: Option<&str>,
  cdses: &[String],
  annotation: Option<&Path>,
  aa_root_sequence: Option<&Path>,
) -> Result<(), Report> {
  if translations.is_none() && cdses.is_empty() && annotation.is_none() && aa_root_sequence.is_none() {
    return Ok(());
  }

  let Some(template) = translations else {
    return make_error!("--translations is required when using --cdses, --annotation, or --aa-root-sequence");
  };

  #[expect(
    clippy::literal_string_with_formatting_args,
    reason = "the braces are a path template placeholder, not a format argument"
  )]
  if !template_has_cds_placeholder(template) {
    return make_error!("--translations must contain a CDS placeholder ('{{cds}}' or '%GENE')");
  }

  if cdses.is_empty() && annotation.is_none() {
    return make_error!("--cdses must list at least one CDS, or pass --annotation to derive the CDS set");
  }

  validate_file_arg("--annotation", annotation)?;
  validate_file_arg("--aa-root-sequence", aa_root_sequence)?;

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

  let read_alphabet = Alphabet::new(AlphabetName::Aa)?;
  let records = read_many_fasta_path(&[path], &read_alphabet)?;
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
