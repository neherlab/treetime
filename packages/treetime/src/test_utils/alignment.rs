use crate::seq::alignment::common_length;
use eyre::Report;
use treetime_primitives::AlignmentRecord;

pub(crate) fn get_common_length(aln: &[AlignmentRecord]) -> Result<usize, Report> {
  common_length(aln.iter().map(|record| (record.name.as_str(), &record.seq)))
}
