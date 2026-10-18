use crate::commands::shared::alignment::{AlignmentArgs, PairedAlignment, pair_alignment};
use crate::commands::shared::alphabet::AlphabetArgs;
use crate::commands::shared::gap_fill::GapFillArgs;
use crate::commands::shared::tree_input::read_input_tree;
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use std::path::Path;
use treetime::alphabet::alphabet::Alphabet;
use treetime::ancestral::attach::{MissingLeafLimit, complete_alignment_for_leaves};
use treetime::cancel::Cancel;
use treetime::progress::{LogSink, StageSink};
use treetime::seq::alignment::{AncestralInput, EdgeSeqInput, SequenceSource};
use treetime::seq::gap_fill::{GapFill, apply_gap_fill};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::fasta_read_file;
use treetime_io::nwk::{NewickDialect, NwkParse};

pub(crate) fn read_sequence_sources(
  alignment: &AlignmentArgs,
  settings: &SourceSettings<'_>,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> Result<Vec<AlignmentSource>, Report> {
  alignment
    .paths()?
    .iter()
    .map(|path| read_sequence_source(path, settings, graph, names, log))
    .collect()
}

pub(crate) fn read_tree_and_sources(
  args: &SequenceInputArgs<'_>,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<TreeAndSources, Report> {
  args.alignment.required_paths()?;
  let alphabet = Alphabet::new(args.alphabet_args.alphabet_name().unwrap_or_default())?;

  cancel.check()?;
  stages.report("Reading input", 0.0, "");
  let tree = read_input_tree(args.tree, args.tree_dialect, log)?;
  let names = tree.names();
  let settings = SourceSettings {
    alphabet: &alphabet,
    gap_fill: Some(args.gap_fill_args.effective_gap_fill()),
    missing_leaves: MissingLeafLimit::Overridable {
      ignored: args.ignore_missing_alns,
    },
  };
  let sources = read_sequence_sources(args.alignment, &settings, &tree.graph, &names, log)?;
  Ok(TreeAndSources { tree, sources })
}

pub(crate) fn ancestral_input(tree: NwkParse, sources: Vec<SequenceSource>) -> AncestralInput {
  let edges = tree
    .branch_lengths
    .into_iter()
    .map(|(key, branch_length)| (key, EdgeSeqInput { branch_length }))
    .collect();
  AncestralInput {
    graph: tree.graph,
    edges,
    sources,
  }
}

pub(crate) struct SourceSettings<'a> {
  pub(crate) alphabet: &'a Alphabet,
  pub(crate) gap_fill: Option<GapFill>,
  pub(crate) missing_leaves: MissingLeafLimit,
}

pub(crate) struct AlignmentSource {
  pub(crate) source: SequenceSource,
  pub(crate) descs: BTreeMap<GraphNodeKey, Option<String>>,
  pub(crate) n_records: usize,
}

pub(crate) struct SequenceInputArgs<'a> {
  pub(crate) alignment: &'a AlignmentArgs,
  pub(crate) tree: &'a Path,
  pub(crate) tree_dialect: NewickDialect,
  pub(crate) alphabet_args: &'a AlphabetArgs,
  pub(crate) gap_fill_args: &'a GapFillArgs,
  pub(crate) ignore_missing_alns: bool,
}

pub(crate) struct TreeAndSources {
  pub(crate) tree: NwkParse,
  pub(crate) sources: Vec<AlignmentSource>,
}

fn read_sequence_source(
  path: &Path,
  settings: &SourceSettings<'_>,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> Result<AlignmentSource, Report> {
  let alphabet = settings.alphabet;
  let mut records = fasta_read_file(path, alphabet)?;
  if let Some(gap_fill) = settings.gap_fill {
    for record in &mut records {
      apply_gap_fill(&mut record.seq, gap_fill, alphabet.gap(), alphabet.unknown());
    }
  }
  let n_records = records.len();
  let PairedAlignment { mut sequences, descs } = pair_alignment(records, path, graph, names, log);
  let in_file = || format!("In alignment '{}'", path.display());
  let length = sequences.common_length().wrap_err_with(in_file)?;
  complete_alignment_for_leaves(
    graph,
    &mut sequences.nodes,
    length,
    alphabet,
    settings.missing_leaves,
    &format!("alignment '{}'", path.display()),
    log,
  )
  .wrap_err_with(in_file)?;
  let mask = sequences.mask(length, alphabet);
  Ok(AlignmentSource {
    source: SequenceSource {
      path: path.to_path_buf(),
      alphabet: alphabet.clone(),
      length,
      nodes: sequences.nodes,
      mask,
    },
    descs,
    n_records,
  })
}
