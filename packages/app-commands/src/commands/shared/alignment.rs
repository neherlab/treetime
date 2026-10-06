use crate::commands::shared::input_warnings::warn_duplicate_names;
use clap::ValueHint;
use deser::{Deserialize, Serialize};
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime::progress::{LogSink, RunWarningKind};
use treetime::seq::alignment::LeafSequences;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pair_by_name::pair_by_name;
use treetime_io::fasta::{FastaRecord, fasta_read_file};
use treetime_primitives::{AlphabetLike, Seq};
#[cfg(feature = "clap")]
use treetime_schema::schema_defaults;

/// Sequence alignment input shared by all commands that read sequences.
///
/// One flag name (`--alignment`, short `-a`, alias `--aln`) serves every command. Multiple files are
/// accepted; their records form one alignment. The path `-` reads uncompressed FASTA from standard
/// input.
#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
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

pub(crate) fn pair_alignment(
  records: Vec<FastaRecord>,
  paths: &[PathBuf],
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> PairedAlignment {
  let pairing = pair_by_name(
    graph.get_leaves().map(|leaf| leaf.key()),
    names,
    records
      .into_iter()
      .map(|record| (record.seq_name, (record.seq, record.desc))),
  );
  let files = paths.iter().map(|path| format!("'{}'", path.display())).join(", ");
  warn_duplicate_names(
    log,
    RunWarningKind::DuplicateSequenceNames,
    &format!("The alignment {files} has more than one sequence named"),
    "TreeTime uses the first sequence of each name.",
    &pairing.duplicate_entry_names,
  );
  let (seqs, descs): (BTreeMap<GraphNodeKey, Seq>, BTreeMap<GraphNodeKey, Option<String>>) = pairing
    .by_node
    .into_iter()
    .map(|(key, (seq, desc))| ((key, seq), (key, desc)))
    .unzip();
  let unmatched = pairing
    .unmatched
    .into_iter()
    .map(|(name, (seq, _))| (name, seq))
    .collect();
  PairedAlignment {
    sequences: LeafSequences::new(names, seqs, unmatched),
    descs,
  }
}

pub(crate) struct PairedAlignment {
  pub(crate) sequences: LeafSequences,
  pub(crate) descs: BTreeMap<GraphNodeKey, Option<String>>,
}

pub(crate) fn read_alignment<A: AlphabetLike>(paths: &[PathBuf], alphabet: &A) -> Result<Vec<FastaRecord>, Report> {
  paths
    .iter()
    .map(|path| fasta_read_file(path, alphabet))
    .flatten_ok()
    .collect()
}
