use crate::commands::homoplasy::args::{TreetimeHomoplasyArgs, homoplasy_ancestral_params};
use crate::commands::homoplasy::drms::DrmTable;
use crate::commands::homoplasy::report::render_homoplasy_report;
use crate::commands::homoplasy::result::HomoplasyResult;
use crate::commands::homoplasy::summary::{ResultContext, homoplasy_result};
use crate::commands::shared::ancestral_trees::{AncestralOutputMaps, AncestralTrees, write_ancestral_trees};
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use crate::commands::shared::sequence_inputs::{AncestralReadInputs, SequenceInputArgs, read_nwk_fasta};
use app_output::mutation_filter::UnknownMutationFilter;
use app_output::output_plan::{CommandKind, OutputSelection, ResolvedOutputs};
use eyre::Report;
use std::io::Write;
use treetime::ancestral;
use treetime::branch_lengths::branch_lengths_or_zero;
use treetime::cancel::Cancel;
use treetime::homoplasy;
use treetime::homoplasy::pipeline::{HomoplasyInput, HomoplasyParams};
use treetime::partition::marginal::sample::SampleMode;
use treetime::progress::{LogSink, StageSink};
use treetime::progress_info;
use treetime_utils::io::file::write_file_with;
use treetime_utils::io::json::{JsonPretty, json_write_file};

pub fn run_homoplasy(
  args: &TreetimeHomoplasyArgs,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<HomoplasyResult, Report> {
  let resolved = args.resolve_outputs()?;
  let drms = args.drms.as_deref().map(DrmTable::read_file).transpose()?;

  let sequence_args = SequenceInputArgs {
    alignment: &args.alignment,
    tree: &args.tree,
    alphabet_args: &args.alphabet_args,
    gap_fill_args: &args.gap_fill_args,
    ignore_missing_alns: args.ignore_missing_alns,
  };
  let AncestralReadInputs { mut input, descs } = read_nwk_fasta(&sequence_args, cancel, stages, log)?;
  progress_info!(
    log,
    "Read {} sequences of length {} and a tree with {} leaves",
    descs.len(),
    input.mask.len(),
    input.graph.num_leaves()
  );
  for edge in input.edges.values_mut() {
    edge.branch_length = edge.branch_length.map(|branch_length| branch_length * args.rescale);
  }
  let names = input.names();
  let branch_lengths = input.branch_lengths();
  let alphabet = input.alphabet.clone();
  let topology_order = args.topology_order.resolve_topology_order(&input.graph, &names, None)?;

  let random_step = (args.sample_from_profile != SampleMode::Argmax).then_some("Sampling from the profile");
  let seed = args.seed_args.resolve(random_step, log);
  let params = homoplasy_ancestral_params(args, seed);
  let output = ancestral::pipeline::run(&params, input, None, cancel, stages, log).map_err(|err| err.into_report())?;
  let mut graph = output.graph;
  let raw_mutations = output.edge_mutations;
  let bridged_mutations =
    UnknownMutationFilter::hiding_unknown(output.ambiguous_char).reported_edge_mutations(&graph, raw_mutations.clone())?;

  cancel.check()?;
  stages.report("Counting recurrent mutations", 0.8, "");
  let statistics = homoplasy::pipeline::run(
    &HomoplasyParams {
      constant_sites: args.constant_sites,
    },
    &HomoplasyInput {
      graph: &graph,
      bridged_mutations: &bridged_mutations,
      raw_mutations: &raw_mutations,
      branch_lengths: &branch_lengths_or_zero(&branch_lengths),
      alphabet: &alphabet,
      sequence_length: output.sequence_length,
    },
  )
  .map_err(|err| err.into_report())?;
  let result = homoplasy_result(
    &statistics,
    &ResultContext {
      names: &names,
      drms: drms.as_ref(),
      zero_based: args.zero_based,
    },
  );
  let report = render_homoplasy_report(&result, args.num_mut, args.detailed);

  stages.report("Writing output", 0.9, "");
  write_homoplasy_outputs(&resolved, &result, &report)?;
  topology_order.apply(&mut graph, &names, &branch_lengths)?;
  let maps = AncestralOutputMaps {
    root_sequence: output.root_sequence,
    edge_mutations: bridged_mutations,
  };
  let trees = AncestralTrees {
    graph: &graph,
    names: &names,
    branch_lengths: &branch_lengths,
    maps: &maps,
    amino_acids: None,
  };
  write_ancestral_trees(&trees, None, &resolved, CommandKind::Homoplasy, log)?;

  for line in report.lines() {
    progress_info!(log, "{line}");
  }
  stages.report("Done", 1.0, "");
  Ok(result)
}

fn write_homoplasy_outputs(resolved: &ResolvedOutputs, result: &HomoplasyResult, report: &str) -> Result<(), Report> {
  if let Some(path) = resolved.path(OutputSelection::HomoplasyStats) {
    json_write_file(path, result, JsonPretty(true))?;
  }
  if let Some(path) = resolved.path(OutputSelection::HomoplasyReport) {
    write_file_with(path, |writer| Ok(writer.write_all(report.as_bytes())?))?;
  }
  Ok(())
}
