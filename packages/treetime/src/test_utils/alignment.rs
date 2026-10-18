use crate::alphabet::alphabet::Alphabet;
use crate::seq::alignment::common_length;
use eyre::Report;
use std::collections::BTreeMap;
use std::io::BufReader;
use std::ops::Range;
use treetime_io::csv::default_name_candidates;
use treetime_io::dates_csv::{DateConstraint, metadata_read};
use treetime_io::fasta::fasta_read;
use treetime_primitives::{AlignmentRecord, Seq};
use treetime_utils::io::compression::{CompressionType, Decompressor};
use treetime_utils::make_internal_report;

#[derive(Clone, Copy, Debug)]
pub(crate) enum TestDataset {
  Zika20,
  Zika86,
  Ebola20,
  Flu20,
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum Columns {
  All,
  Leading(usize),
}

struct DatasetFiles {
  tree: &'static str,
  alignment_xz: &'static [u8],
  metadata: &'static str,
}

macro_rules! dataset_files {
  ($dir:literal) => {
    DatasetFiles {
      tree: include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/", $dir, "/tree.nwk")),
      alignment_xz: include_bytes!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../data/",
        $dir,
        "/aln.fasta.xz"
      )),
      metadata: include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../data/",
        $dir,
        "/metadata.tsv"
      )),
    }
  };
}

const ZIKA_20: DatasetFiles = dataset_files!("zika/20");
const ZIKA_86: DatasetFiles = dataset_files!("zika/86");
const EBOLA_20: DatasetFiles = dataset_files!("ebola/20");
const FLU_20: DatasetFiles = dataset_files!("flu/h3n2/20");

impl TestDataset {
  fn files(self) -> &'static DatasetFiles {
    match self {
      Self::Zika20 => &ZIKA_20,
      Self::Zika86 => &ZIKA_86,
      Self::Ebola20 => &EBOLA_20,
      Self::Flu20 => &FLU_20,
    }
  }

  pub(crate) fn tree(self) -> &'static str {
    self.files().tree
  }

  pub(crate) fn alignment(self) -> Result<Vec<AlignmentRecord>, Report> {
    fasta_xz_read(self.files().alignment_xz)
  }

  pub(crate) fn alignment_columns(self, selection: Columns) -> Result<Vec<AlignmentRecord>, Report> {
    let records = self.alignment()?;
    match selection {
      Columns::All => Ok(records),
      Columns::Leading(count) => {
        let length = get_common_length(&records)?;
        Ok(columns(&records, 0..length.min(count)))
      },
    }
  }

  pub(crate) fn dates(self) -> Result<BTreeMap<String, Option<DateConstraint>>, Report> {
    let table = metadata_read(
      self.files().metadata.as_bytes(),
      None,
      &['\t'],
      &default_name_candidates(),
      None,
      None,
    )?;
    Ok(table.dates()?.into_iter().collect())
  }
}

pub(crate) fn fasta_xz_read(bytes: &[u8]) -> Result<Vec<AlignmentRecord>, Report> {
  let decompressor = Decompressor::new(bytes, &CompressionType::Xz)?;
  Ok(
    fasta_read(BufReader::new(decompressor), &Alphabet::default())?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect(),
  )
}

pub(crate) fn get_common_length(aln: &[AlignmentRecord]) -> Result<usize, Report> {
  common_length(aln.iter().map(|record| (record.name.as_str(), &record.seq)))
}

pub(crate) fn split_at_gap_free_column(records: &[AlignmentRecord]) -> Result<[Vec<AlignmentRecord>; 2], Report> {
  let length = get_common_length(records)?;
  let gap = Alphabet::default().gap();
  let gap_free = |column: usize| records.iter().all(|record| record.seq[column] != gap);
  let cut = (length.div_ceil(2)..length)
    .find(|&cut| gap_free(cut - 1) && gap_free(cut))
    .ok_or_else(|| make_internal_report!("The alignment has no gap-free cut"))?;
  Ok([columns(records, 0..cut), columns(records, cut..length)])
}

pub(crate) fn columns(records: &[AlignmentRecord], range: Range<usize>) -> Vec<AlignmentRecord> {
  records
    .iter()
    .map(|record| AlignmentRecord {
      name: record.name.clone(),
      seq: Seq::from(&record.seq.as_slice()[range.clone()]),
    })
    .collect()
}
