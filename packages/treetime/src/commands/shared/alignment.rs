#[cfg(feature = "clap")]
use clap::ValueHint;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::PathBuf;

/// Sequence alignment input shared by all commands that read sequences.
///
/// One flag name (`--alignment`, short `-a`, alias `--aln`) is used across every command, replacing
/// the earlier mix of positional arguments and `--aln`. Multiple files are accepted; each is one
/// input alignment. When the list is empty, callers read uncompressed FASTA from standard input.
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct AlignmentArgs {
  /// Path to one or multiple FASTA files with aligned input sequences
  ///
  /// Accepts plain or compressed FASTA files. If a compressed fasta file is provided, it will be transparently
  /// decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file
  /// extension. If there's multiple input files, then different files can have different compression formats.
  ///
  /// If no input files provided, the plain fasta input is read from standard input (stdin).
  ///
  /// See: https://en.wikipedia.org/wiki/FASTA_format
  #[cfg_attr(
    feature = "clap",
    clap(
      long = "alignment",
      short = 'a',
      visible_alias = "aln",
      value_hint = ValueHint::FilePath,
      value_name = "FILEPATH",
      display_order = 1,
    )
  )]
  pub alignment: Vec<PathBuf>,
}
