#[cfg(test)]
mod __tests__;

use eyre::{Report, WrapErr};
use itertools::Itertools;
use percent_encoding::percent_decode_str;
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;
use treetime_utils::io::fs::read_file_to_string;
use treetime_utils::make_error;

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

#[derive(Clone, Debug)]
pub struct GffCdsFeature {
  pub name: String,
  pub seqid: String,
  pub segments: Vec<GffCdsSegment>,
  pub strand: String,
}

#[derive(Clone, Debug)]
pub struct GffCdsSegment {
  pub start: i64,
  pub end: i64,
}

pub fn read_gff3_cds_features(path: &Path) -> Result<Vec<GffCdsFeature>, Report> {
  let contents = read_file_to_string(path)?;
  parse_gff3_cds_features(&contents, path)
}

pub fn read_gff3_cds_features_filtered(path: &Path, cdses: &[String]) -> Result<Vec<GffCdsFeature>, Report> {
  let contents = read_file_to_string(path)?;
  let features = parse_gff3_cds_features(&contents, path)?;

  if cdses.is_empty() {
    return Ok(features);
  }

  let wanted: BTreeSet<&str> = cdses.iter().map(String::as_str).collect();
  for cds in &wanted {
    if !features.iter().any(|f| f.name == *cds) {
      return make_error!(
        "--annotation '{}' does not contain a CDS feature for CDS '{cds}'",
        path.display()
      );
    }
  }

  Ok(
    features
      .into_iter()
      .filter(|f| wanted.contains(f.name.as_str()))
      .collect(),
  )
}

fn parse_gff3_cds_features(contents: &str, path: &Path) -> Result<Vec<GffCdsFeature>, Report> {
  let mut raw_features: BTreeMap<String, Vec<RawCdsRow>> = BTreeMap::new();

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

    raw_features.entry(cds).or_default().push(RawCdsRow {
      seqid: cols[0].to_owned(),
      start,
      end,
      strand,
    });
  }

  raw_features
    .into_iter()
    .map(|(name, mut rows)| {
      rows.sort_by_key(|row| (row.start, row.end));
      let strand = rows.first().expect("features are non-empty").strand.clone();
      let seqid = rows.first().expect("features are non-empty").seqid.clone();

      if rows.iter().any(|row| row.strand != strand) {
        return make_error!("CDS feature rows for CDS '{name}' use multiple strands");
      }

      if rows.iter().any(|row| row.seqid != seqid) {
        return make_error!("CDS feature rows for CDS '{name}' span multiple sequence IDs");
      }

      // 5'-to-3' order: ascending for plus strand (already sorted), reversed for minus strand.
      if strand == "-" {
        rows.reverse();
      }

      let cds_length: i64 = rows.iter().map(|row| row.end - row.start + 1).sum();
      if cds_length % 3 != 0 {
        return make_error!(
          "CDS feature for CDS '{name}' has nucleotide length {cds_length}, which is not a multiple of 3. \
           Check that the annotation matches the translations."
        );
      }

      let segments = rows
        .iter()
        .map(|row| GffCdsSegment {
          start: row.start,
          end: row.end,
        })
        .collect();

      Ok(GffCdsFeature {
        name,
        seqid,
        segments,
        strand,
      })
    })
    .collect()
}

pub fn resolve_cds_name(attrs: &BTreeMap<String, String>) -> Option<String> {
  NAME_ATTRS_CDS.iter().find_map(|key| attrs.get(*key).cloned())
}

pub fn parse_gff_attributes(raw: &str) -> Result<BTreeMap<String, String>, Report> {
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
struct RawCdsRow {
  seqid: String,
  start: i64,
  end: i64,
  strand: String,
}
