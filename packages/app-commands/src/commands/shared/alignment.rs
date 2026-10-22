use crate::commands::shared::argument_error::argument_error;
use crate::commands::shared::input_warnings::warn_duplicate_names;
use clap::ValueHint;
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use treetime::progress::{LogSink, RunWarningKind};
use treetime::seq::alignment::LeafSequences;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_graph::pair_by_name::pair_by_name;
use treetime_io::fasta::FastaRecord;
use treetime_primitives::Seq;
#[cfg(feature = "clap")]
use treetime_schema::schema_defaults;
use treetime_utils::io::file::is_path_stdin;

/// Sequence alignment input shared by all commands that accept an alignment.
///
/// Each `--alignment` file is one partition, numbered by its position in the list (1, 2, ...). The path `-` reads
/// uncompressed FASTA from standard input and may appear at most once.
#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct AlignmentArgs {
  /// Aligned FASTA input, plain or compressed (`gz`, `bz2`, `xz`, `zst`); compression is detected by extension.
  ///
  /// Each `--alignment` file is one partition, numbered by its position in the list (1, 2, ...). Use `-` to read
  /// uncompressed FASTA from standard input; `-` may appear at most once.
  #[cfg_attr(
    feature = "clap",
    clap(
      long,
      short = 'a',
      value_hint = ValueHint::FilePath,
      value_name = "FILEPATH",
      display_order = 1,
      help_heading = "Input data",
    )
  )]
  #[schemars(extend("x-path" = "input"))]
  pub alignment: Vec<PathBuf>,
}

impl AlignmentArgs {
  pub fn paths(&self) -> Result<&[PathBuf], Report> {
    let n_stdin = self.alignment.iter().filter(|path| is_path_stdin(path)).count();
    if n_stdin > 1 {
      return Err(argument_error(
        &["alignment"],
        format!("Standard input ('-') can be read only once, but --alignment lists it {n_stdin} times"),
      ));
    }
    Ok(&self.alignment)
  }

  pub fn required_paths(&self) -> Result<&[PathBuf], Report> {
    let paths = self.paths()?;
    if paths.is_empty() {
      return Err(argument_error(
        &["alignment"],
        "--alignment is required: pass one or more FASTA files, or '-' to read standard input",
      ));
    }
    Ok(paths)
  }
}

pub(crate) fn pair_alignment(
  records: Vec<FastaRecord>,
  path: &Path,
  graph: &Graph,
  names: &NodeNames,
  log: &dyn LogSink,
) -> PairedAlignment {
  let pairing = pair_by_name(
    graph.get_leaves().map(|leaf| leaf.key()),
    names,
    records
      .into_iter()
      .map(|record| (record.seq_name, (record.seq, record.desc))),
  );
  warn_duplicate_names(
    log,
    RunWarningKind::DuplicateSequenceNames,
    &format!("The alignment '{}' has more than one sequence named", path.display()),
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
    sequences: LeafSequences::new(graph, names, seqs, unmatched),
    descs,
  }
}

pub(crate) struct PairedAlignment {
  pub(crate) sequences: LeafSequences,
  pub(crate) descs: BTreeMap<GraphNodeKey, Option<String>>,
}
