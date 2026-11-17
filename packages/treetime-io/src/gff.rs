use eyre::{Report, WrapErr};
use itertools::Itertools;
use noodles_gff::directive_buf::key::{FASTA, SEQUENCE_REGION};
use noodles_gff::directive_buf::value::SequenceRegion;
use noodles_gff::feature::record::{Phase, Strand};
use noodles_gff::{Line, Record};
use std::collections::BTreeMap;
use std::io::Read;
use std::path::Path;
use std::str::from_utf8;
use treetime_utils::io::file::read_file_with;
use treetime_utils::make_error;

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

pub fn gff3_cds_read_file(filepath: impl AsRef<Path>) -> Result<GffCdsAnnotation, Report> {
  read_file_with(filepath, gff3_cds_read)
}

pub fn gff3_cds_read(mut reader: impl Read) -> Result<GffCdsAnnotation, Report> {
  let mut contents = String::new();
  reader.read_to_string(&mut contents).wrap_err("When reading GFF3")?;
  parse_gff3_cds_annotation(&contents).wrap_err("When reading GFF3")
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GffCdsAnnotation {
  pub features: Vec<GffCdsFeature>,
  pub sequence_regions: Vec<GffSequenceRegion>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GffSequenceRegion {
  pub seqid: String,
  pub start: usize,
  pub end: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GffCdsFeature {
  pub label: String,
  pub id: Option<String>,
  pub name: String,
  pub seqid: String,
  pub strand: GffStrand,
  pub segments: Vec<GffCdsSegment>,
}

impl GffCdsFeature {
  pub fn coding_length(&self) -> usize {
    self.segments.iter().map(GffCdsSegment::length).sum()
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GffStrand {
  Forward,
  Reverse,
}

impl GffStrand {
  pub const fn symbol(self) -> &'static str {
    match self {
      Self::Forward => "+",
      Self::Reverse => "-",
    }
  }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GffCdsSegment {
  pub start: usize,
  pub end: usize,
  pub phase: GffPhase,
}

impl GffCdsSegment {
  pub const fn length(&self) -> usize {
    self.end - self.start + 1
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GffPhase {
  Zero,
  One,
  Two,
}

impl GffPhase {
  pub const fn offset(self) -> usize {
    match self {
      Self::Zero => 0,
      Self::One => 1,
      Self::Two => 2,
    }
  }
}

fn parse_gff3_cds_annotation(contents: &str) -> Result<GffCdsAnnotation, Report> {
  let mut rows = vec![];
  let mut sequence_regions = vec![];

  for (index, text) in contents.lines().enumerate() {
    let line_no = index + 1;
    if text.trim().is_empty() {
      continue;
    }
    let line = read_gff_line(text).wrap_err_with(|| format!("Invalid GFF3 line {line_no}"))?;
    if let Some(directive) = line.as_directive() {
      let key: &[u8] = directive.key();
      if key == FASTA {
        break;
      }
      if key == SEQUENCE_REGION {
        sequence_regions.push(parse_sequence_region(directive.value().map(|value| &**value), line_no)?);
      }
      continue;
    }
    if let Some(record) = line.as_record() {
      let record = record.wrap_err_with(|| format!("Invalid GFF3 line {line_no}"))?;
      if let Some(row) = parse_cds_row(&record, text, line_no)? {
        rows.push(row);
      }
    }
  }

  if rows.is_empty() {
    return make_error!(
      "GFF3 file has no CDS rows: TreeTime reads the coding regions from CDS rows, and ignores gene and other feature rows"
    );
  }

  Ok(GffCdsAnnotation {
    features: assemble_features(rows)?,
    sequence_regions,
  })
}

fn read_gff_line(text: &str) -> Result<Line, Report> {
  let mut reader = noodles_gff::io::Reader::new(text.as_bytes());
  let mut line = Line::default();
  reader.read_line(&mut line)?;
  Ok(line)
}

fn parse_sequence_region(value: Option<&[u8]>, line_no: usize) -> Result<GffSequenceRegion, Report> {
  let value = from_utf8(value.unwrap_or_default())?;
  let region: SequenceRegion = value
    .parse()
    .wrap_err_with(|| format!("Invalid ##sequence-region directive '{value}' on GFF3 line {line_no}"))?;
  Ok(GffSequenceRegion {
    seqid: from_utf8(region.reference_sequence_name())?.to_owned(),
    start: usize::from(region.start()),
    end: usize::from(region.end()),
  })
}

fn parse_cds_row(record: &Record<'_>, text: &str, line_no: usize) -> Result<Option<CdsRow>, Report> {
  if &**record.ty() != b"CDS" {
    return Ok(None);
  }

  let attrs = record_attributes(record).wrap_err_with(|| format!("Invalid GFF3 attributes on line {line_no}"))?;
  let Some(name) = resolve_cds_name(&attrs) else {
    return make_error!(
      "CDS on GFF3 line {line_no} has no name attribute and no ID: TreeTime needs one of them to name the CDS"
    );
  };

  let start = usize::from(
    record
      .start()
      .wrap_err_with(|| format!("Invalid CDS start on GFF3 line {line_no}"))?,
  );
  let end = usize::from(
    record
      .end()
      .wrap_err_with(|| format!("Invalid CDS end on GFF3 line {line_no}"))?,
  );
  if end < start {
    return make_error!("Invalid CDS coordinates {start}-{end} for CDS '{name}' on GFF3 line {line_no}");
  }

  let strand = match record.strand() {
    Ok(Strand::Forward) => GffStrand::Forward,
    Ok(Strand::Reverse) => GffStrand::Reverse,
    Ok(Strand::None | Strand::Unknown) | Err(_) => {
      return make_error!(
        "Invalid CDS strand '{}' for CDS '{name}' on GFF3 line {line_no}: a CDS needs strand '+' or '-'",
        raw_column(text, 6)
      );
    },
  };

  let phase = match record.phase() {
    Some(Ok(Phase::Zero)) => GffPhase::Zero,
    Some(Ok(Phase::One)) => GffPhase::One,
    Some(Ok(Phase::Two)) => GffPhase::Two,
    Some(Err(_)) | None => {
      return make_error!(
        "Invalid CDS phase '{}' for CDS '{name}' on GFF3 line {line_no}: a CDS needs phase 0, 1 or 2",
        raw_column(text, 7)
      );
    },
  };

  Ok(Some(CdsRow {
    line_no,
    id: attrs.get("ID").cloned(),
    name,
    seqid: from_utf8(&record.reference_sequence_name())?.to_owned(),
    strand,
    segment: GffCdsSegment { start, end, phase },
  }))
}

fn raw_column(text: &str, index: usize) -> &str {
  text.split('\t').nth(index).unwrap_or_default()
}

fn record_attributes(record: &Record<'_>) -> Result<BTreeMap<String, String>, Report> {
  record
    .attributes()
    .iter()
    .map(|field| {
      let (tag, value) = field?;
      let value = match value {
        noodles_gff::record::attributes::field::Value::String(value) => from_utf8(&value)?.to_owned(),
        noodles_gff::record::attributes::field::Value::Array(values) => values
          .iter()
          .map(|value| Ok(from_utf8(&value)?.to_owned()))
          .collect::<Result<Vec<String>, Report>>()?
          .join(","),
      };
      Ok((from_utf8(&tag)?.to_owned(), value))
    })
    .collect()
}

pub(super) fn resolve_cds_name(attrs: &BTreeMap<String, String>) -> Option<String> {
  NAME_ATTRS_CDS.iter().find_map(|key| attrs.get(*key).cloned())
}

fn assemble_features(rows: Vec<CdsRow>) -> Result<Vec<GffCdsFeature>, Report> {
  let groups = rows
    .into_iter()
    .fold(BTreeMap::<FeatureKey, Vec<CdsRow>>::new(), |mut groups, row| {
      let key = row.id.clone().map_or(FeatureKey::Row(row.line_no), FeatureKey::Id);
      groups.entry(key).or_default().push(row);
      groups
    });

  let features = groups
    .into_values()
    .map(assemble_feature)
    .collect::<Result<Vec<_>, Report>>()?;

  label_features(features)
}

fn assemble_feature(mut rows: Vec<CdsRow>) -> Result<RawFeature, Report> {
  rows.sort_by_key(|row| (row.segment.start, row.segment.end));
  let Some(first) = rows.iter().min_by_key(|row| row.line_no) else {
    return make_error!("CDS feature has no rows");
  };
  let name = first.name.clone();
  let id = first.id.clone();
  let strand = first.strand;
  let seqid = first.seqid.clone();
  let lines = rows.iter().map(|row| row.line_no).sorted().collect_vec();

  if rows.iter().any(|row| row.strand != strand) {
    return make_error!("CDS feature rows for CDS '{name}' use multiple strands");
  }
  if rows.iter().any(|row| row.seqid != seqid) {
    return make_error!("CDS feature rows for CDS '{name}' span multiple sequence IDs");
  }
  if strand == GffStrand::Reverse {
    rows.reverse();
  }

  let segments = rows.into_iter().map(|row| row.segment).collect_vec();
  let cds_length: usize = segments.iter().map(GffCdsSegment::length).sum();
  if !cds_length.is_multiple_of(3) {
    return make_error!(
      "CDS feature for CDS '{name}' has nucleotide length {cds_length}, which is not a multiple of 3. \
       Check that the annotation matches the translations."
    );
  }

  Ok(RawFeature {
    feature: GffCdsFeature {
      label: name.clone(),
      id,
      name,
      seqid,
      strand,
      segments,
    },
    lines,
  })
}

fn label_features(features: Vec<RawFeature>) -> Result<Vec<GffCdsFeature>, Report> {
  let name_counts = features.iter().map(|raw| raw.feature.name.as_str()).counts();

  let labeled = features
    .iter()
    .map(|raw| {
      if name_counts[raw.feature.name.as_str()] == 1 {
        return Ok(raw.feature.name.clone());
      }
      raw.feature.id.clone().map_or_else(
        || {
          let lines = features
            .iter()
            .filter(|other| other.feature.name == raw.feature.name)
            .flat_map(|other| other.lines.iter())
            .sorted()
            .join(", ");
          make_error!(
            "CDS name '{}' is used by several CDS features (GFF3 lines {lines}), and the CDS on GFF3 line {} has no ID \
             to tell them apart",
            raw.feature.name,
            raw.lines[0]
          )
        },
        Ok,
      )
    })
    .collect::<Result<Vec<String>, Report>>()?;

  if let Some(label) = labeled.iter().duplicates().next() {
    return make_error!("CDS label '{label}' names more than one CDS feature: each CDS needs a unique name or ID");
  }

  Ok(
    features
      .into_iter()
      .zip(labeled)
      .map(|(raw, label)| GffCdsFeature { label, ..raw.feature })
      .sorted_by(|a, b| a.label.cmp(&b.label))
      .collect(),
  )
}

#[derive(Clone, Debug)]
struct CdsRow {
  line_no: usize,
  id: Option<String>,
  name: String,
  seqid: String,
  strand: GffStrand,
  segment: GffCdsSegment,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum FeatureKey {
  Id(String),
  Row(usize),
}

#[derive(Clone, Debug)]
struct RawFeature {
  feature: GffCdsFeature,
  lines: Vec<usize>,
}
