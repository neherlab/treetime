use crate::commands::shared::alignment::{AlignmentArgs, read_alignment, sequence_descriptions};
use crate::commands::shared::alphabet::AlphabetArgs;
use crate::commands::shared::gap_fill::GapFillArgs;
use eyre::Report;
use std::collections::BTreeMap;
use std::path::Path;
use treetime::alphabet::alphabet::Alphabet;
use treetime::ancestral::attach::complete_alignment_for_leaves;
use treetime::ancestral::mask::create_mask;
use treetime::cancel::Cancel;
use treetime::make_error;
use treetime::progress::{LogSink, StageSink};
use treetime::seq::alignment::{AncestralInput, EdgeSeqInput, get_common_length, node_seq_inputs};
use treetime::seq::gap_fill::apply_gap_fill;
use treetime_io::nwk::nwk_read_file;
use treetime_primitives::AlignmentRecord;

pub(crate) struct SequenceInputArgs<'a> {
  pub(crate) alignment: &'a AlignmentArgs,
  pub(crate) tree: &'a Path,
  pub(crate) alphabet_args: &'a AlphabetArgs,
  pub(crate) gap_fill_args: &'a GapFillArgs,
  pub(crate) ignore_missing_alns: bool,
}

pub(crate) struct AncestralReadInputs {
  pub(crate) input: AncestralInput,
  pub(crate) descs: BTreeMap<String, Option<String>>,
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

  let descs = sequence_descriptions(&aln);

  cancel.check()?;
  stages.report("Parsing tree", 0.1, "");
  let parse = nwk_read_file(args.tree)?;

  let names = parse.names();
  let aln = aln.into_iter().map(AlignmentRecord::from).collect();
  let aln = complete_alignment_for_leaves(&parse.graph, aln, &alphabet, args.ignore_missing_alns, &names, log)?;
  let alignment_length = get_common_length(&aln)?;
  let mask = create_mask(&aln, alignment_length, &alphabet);

  let graph = parse.graph;
  let nodes = node_seq_inputs(&graph, &names, aln);
  let edges = parse
    .branch_lengths
    .into_iter()
    .map(|(key, branch_length)| (key, EdgeSeqInput { branch_length }))
    .collect();
  let input = AncestralInput {
    graph,
    nodes,
    edges,
    alphabet,
    mask,
  };
  Ok(AncestralReadInputs { input, descs })
}
