use crate::commands::shared::alignment::{AlignmentArgs, PairedAlignment, pair_alignment, read_alignment};
use crate::commands::shared::alphabet::AlphabetArgs;
use crate::commands::shared::gap_fill::GapFillArgs;
use crate::commands::shared::tree_input::read_input_tree;
use eyre::Report;
use std::collections::BTreeMap;
use std::path::Path;
use treetime::alphabet::alphabet::Alphabet;
use treetime::ancestral::attach::complete_alignment_for_leaves;
use treetime::cancel::Cancel;
use treetime::make_error;
use treetime::progress::{LogSink, StageSink};
use treetime::seq::alignment::{AncestralInput, EdgeSeqInput};
use treetime::seq::gap_fill::apply_gap_fill;
use treetime_graph::node::GraphNodeKey;
use treetime_io::nwk::NewickDialect;

pub(crate) struct SequenceInputArgs<'a> {
  pub(crate) alignment: &'a AlignmentArgs,
  pub(crate) tree: &'a Path,
  pub(crate) tree_dialect: NewickDialect,
  pub(crate) alphabet_args: &'a AlphabetArgs,
  pub(crate) gap_fill_args: &'a GapFillArgs,
  pub(crate) ignore_missing_alns: bool,
}

pub(crate) struct AncestralReadInputs {
  pub(crate) input: AncestralInput,
  pub(crate) descs: BTreeMap<GraphNodeKey, Option<String>>,
  pub(crate) n_records: usize,
}

pub(crate) fn read_nwk_fasta(
  args: &SequenceInputArgs<'_>,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<AncestralReadInputs, Report> {
  let gap_fill_mode = args.gap_fill_args.effective_gap_fill();
  let alphabet = Alphabet::new(args.alphabet_args.alphabet_name().unwrap_or_default())?;

  cancel.check()?;
  stages.report("Reading input", 0.0, "");

  if args.alignment.alignment.is_empty() {
    return make_error!("--alignment is required: pass one or more FASTA files, or '-' to read standard input");
  }
  let mut aln = read_alignment(&args.alignment.alignment, &alphabet)?;

  for record in &mut aln {
    apply_gap_fill(&mut record.seq, gap_fill_mode, alphabet.gap(), alphabet.unknown());
  }
  let n_records = aln.len();

  cancel.check()?;
  stages.report("Parsing tree", 0.1, "");
  let parse = read_input_tree(args.tree, args.tree_dialect, log)?;

  let names = parse.names();
  let PairedAlignment { mut sequences, descs } =
    pair_alignment(aln, &args.alignment.alignment, &parse.graph, &names, log);
  let alignment_length = sequences.common_length()?;
  complete_alignment_for_leaves(
    &parse.graph,
    &mut sequences.nodes,
    alignment_length,
    &alphabet,
    args.ignore_missing_alns,
    log,
  )?;
  let mask = sequences.mask(alignment_length, &alphabet);

  let edges = parse
    .branch_lengths
    .into_iter()
    .map(|(key, branch_length)| (key, EdgeSeqInput { branch_length }))
    .collect();
  let input = AncestralInput {
    graph: parse.graph,
    nodes: sequences.nodes,
    edges,
    alphabet,
    mask,
  };
  Ok(AncestralReadInputs {
    input,
    descs,
    n_records,
  })
}
