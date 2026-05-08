use crate::make_error;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use treetime_io::fasta::FastaRecord;

pub fn get_common_length(aln: &[FastaRecord]) -> Result<usize, Report> {
  let lengths = aln
    .iter()
    .into_group_map_by(|aln| aln.seq.len())
    .into_iter()
    .collect_vec();

  match lengths[..] {
    [] => Ok(0),
    [(length, _)] => Ok(length),
    _ => {
      let message = lengths
        .into_iter()
        .sorted_by_key(|(length, _)| *length)
        .map(|(length, entries)| {
          let names = entries.iter().map(|aln| format!("    \"{}\"", aln.seq_name)).join("\n");
          format!("Length {length}:\n{names}")
        })
        .join("\n\n");

      make_error!("Sequences are expected to all have the same length, but found the following lengths:\n\n{message}")
    },
  }
  .wrap_err("When calculating length of sequences")
}
