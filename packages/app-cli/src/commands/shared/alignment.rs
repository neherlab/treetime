#[cfg(feature = "clap")]
use clap::ValueHint;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::PathBuf;

/// Sequence alignment input shared by all commands that read sequences.
///
/// One flag name (`--alignment`, short `-a`, alias `--aln`) is used across every command, replacing
/// the earlier mix of positional arguments and `--aln`. Multiple files are accepted; each is one
/// input alignment. When the list is empty, callers read uncompressed FASTA from standard input.
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub(crate) struct AlignmentArgs {
  /// Aligned FASTA input. Accepts multiple plain or compressed (`gz`, `bz2`,
  /// `xz`, `zstd`) files and detects compression by extension. With no files,
  /// reads uncompressed FASTA from standard input.
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
