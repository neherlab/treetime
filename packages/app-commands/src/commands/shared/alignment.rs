#[cfg(feature = "clap")]
use clap::ValueHint;
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime_io::fasta::{FastaRecord, fasta_read_file};
use treetime_primitives::AlphabetLike;

/// Sequence alignment input shared by all commands that read sequences.
///
/// One flag name (`--alignment`, short `-a`, alias `--aln`) serves every command. Multiple files are
/// accepted; their records form one alignment. The path `-` reads uncompressed FASTA from standard
/// input.
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct AlignmentArgs {
  /// Aligned FASTA input. Accepts multiple plain or compressed (`gz`, `bz2`,
  /// `xz`, `zstd`) files and detects compression by extension. The records of all
  /// files form one alignment. Use `-` to read uncompressed FASTA from standard input.
  #[cfg_attr(
    feature = "clap",
    clap(
      long = "alignment",
      short = 'a',
      visible_alias = "aln",
      value_hint = ValueHint::FilePath,
      value_name = "FILEPATH",
      display_order = 1,
      help_heading = "Input data",
    )
  )]
  #[schemars(extend("x-path" = "input"))]
  pub alignment: Vec<PathBuf>,
}

pub(crate) fn sequence_descriptions<'a>(
  records: impl IntoIterator<Item = &'a FastaRecord>,
) -> BTreeMap<String, Option<String>> {
  records.into_iter().fold(BTreeMap::new(), |mut descs, record| {
    descs
      .entry(record.seq_name.clone())
      .or_insert_with(|| record.desc.clone());
    descs
  })
}

pub(crate) fn read_alignment<A: AlphabetLike>(paths: &[PathBuf], alphabet: &A) -> Result<Vec<FastaRecord>, Report> {
  paths
    .iter()
    .map(|path| fasta_read_file(path, alphabet))
    .flatten_ok()
    .collect()
}
