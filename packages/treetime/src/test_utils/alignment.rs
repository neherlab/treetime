use crate::alphabet::alphabet::Alphabet;
use crate::seq::alignment::common_length;
use eyre::Report;
use std::ops::Range;
use std::path::{Path, PathBuf};
use treetime_io::fasta::fasta_read_file;
use treetime_primitives::{AlignmentRecord, Seq};
use treetime_utils::make_internal_report;

pub(crate) fn get_common_length(aln: &[AlignmentRecord]) -> Result<usize, Report> {
  common_length(aln.iter().map(|record| (record.name.as_str(), &record.seq)))
}

pub(crate) fn dataset_dir(name: &str) -> PathBuf {
  Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data").join(name)
}

pub(crate) fn read_dataset_alignment(name: &str) -> Result<Vec<AlignmentRecord>, Report> {
  Ok(
    fasta_read_file(dataset_dir(name).join("aln.fasta.xz"), &Alphabet::default())?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect(),
  )
}

pub(crate) fn split_at_gap_free_column(records: &[AlignmentRecord]) -> Result<[Vec<AlignmentRecord>; 2], Report> {
  let length = get_common_length(records)?;
  let gap = Alphabet::default().gap();
  let gap_free = |column: usize| records.iter().all(|record| record.seq[column] != gap);
  let cut = (length.div_ceil(2)..length)
    .find(|&cut| gap_free(cut - 1) && gap_free(cut))
    .ok_or_else(|| make_internal_report!("The alignment has no gap-free cut"))?;
  let columns = |range: Range<usize>| {
    records
      .iter()
      .map(|record| AlignmentRecord {
        name: record.name.clone(),
        seq: Seq::from(&record.seq.as_slice()[range.clone()]),
      })
      .collect()
  };
  Ok([columns(0..cut), columns(cut..length)])
}
