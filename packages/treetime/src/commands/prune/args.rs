use crate::alphabet::alphabet::AlphabetName;
use clap::{Parser, ValueHint};
use std::fmt::Debug;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct TreetimePruneArgs {
  /// Path to one or multiple FASTA files with aligned input sequences
  ///
  /// Accepts plain or compressed FASTA files. If a compressed fasta file is provided, it will be transparently
  /// decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file
  /// extension. If there's multiple input files, then different files can have different compression formats.
  ///
  /// Use '-' to read uncompressed FASTA from standard input (stdin).
  ///
  /// See: https://en.wikipedia.org/wiki/FASTA_format
  ///
  /// Required if --prune-empty is set.
  #[clap(long = "aln", value_hint = ValueHint::FilePath, value_name = "FILEPATH")]
  pub input_fastas: Vec<PathBuf>,

  /// Name of file containing the tree in newick, nexus, or phylip format.
  #[clap(long, short = 't')]
  #[clap(value_hint = ValueHint::FilePath)]
  pub tree: PathBuf,

  /// Alphabet
  ///
  #[clap(long, short = 'a', value_enum)]
  pub alphabet: Option<AlphabetName>,

  /// Directory to write the output to
  #[clap(long, short = 'O')]
  pub outdir: PathBuf,

  /// Threshold value for pruning of branches
  ///
  /// If set, prune branches with a length below this value
  #[clap(long, short = 's', value_name = "THRESHOLD")]
  pub prune_short: Option<f64>,

  /// Prune empty branches
  ///
  /// If set, prune any branch that does not have a mutation or other state transition mapped to it.
  ///
  /// Requires --aln
  #[clap(long, short = 'e')]
  pub prune_empty: bool,
}
