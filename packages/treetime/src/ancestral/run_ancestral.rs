use crate::cli::treetime_cli::TreetimeAncestralArgs;
use crate::io::fasta::read_many_fasta;
use eyre::Report;

pub fn run_ancestral(ancestral_args: &TreetimeAncestralArgs) -> Result<(), Report> {
  let fasta_records = read_many_fasta(&ancestral_args.input_fastas)?;

  dbg!(fasta_records);

  Ok(())
}
