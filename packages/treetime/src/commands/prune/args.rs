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

  /// List of node names to prune
  ///
  /// List of node names to remove from the tree, comma-separated (,)
  ///
  /// Use --prune-nodes-list-delimiter to specify a different delimiter.
  #[clap(long, short = 'n', value_name = "NODE_NAMES")]
  pub prune_nodes_list: Option<String>,

  /// Name separator for `--prune-nodes-list`
  ///
  /// String used to separate node names in the list given to (--prune-nodes-list). Make sure to correctly quote and escape the delimiter according to your shell.
  #[clap(long, default_value = ",", value_name = "DELIMITER")]
  pub prune_nodes_list_delimiter: char,

  /// File containing list of node names to prune
  ///
  /// Path to a file containing node names to remove from the tree, newline-delimited (\n).
  ///
  /// Use '-' to read from standard input (stdin).
  ///
  /// Use --prune-nodes-list-file-delimiter to specify a different delimiter.
  #[clap(long, short = 'N', value_hint = ValueHint::FilePath, value_name = "FILEPATH")]
  pub prune_nodes_list_file: Option<PathBuf>,

  /// Separator for node names in the list file
  ///
  /// Character or string used to separate node names in the list file.
  #[clap(long, default_value = "\n", value_name = "DELIMITER")]
  pub prune_nodes_list_file_delimiter: char,
}
