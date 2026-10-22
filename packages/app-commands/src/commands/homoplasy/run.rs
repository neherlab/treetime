use crate::commands::homoplasy::args::{TreetimeHomoplasyArgs, homoplasy_ancestral_params};
use crate::commands::homoplasy::drms::DrmTable;
use crate::commands::homoplasy::report::render_homoplasy_report;
use crate::commands::homoplasy::result::HomoplasyResult;
use crate::commands::homoplasy::summary::{ResultContext, homoplasy_result};
use crate::commands::shared::ancestral_trees::{
  AncestralOutputMaps, AncestralPartitionTrees, AncestralTrees, write_ancestral_trees,
};
use crate::commands::shared::partition::log_partitions;
use crate::commands::shared::seed::SeedArgs;
use crate::commands::shared::sequence_inputs::{
  AlignmentSource, TreeAndSources, ancestral_input, read_tree_and_sources,
};
use app_output::command_outputs::CommandKind;
use app_output::mutation_filter::UnknownMutationFilter;
use app_output::output_plan::{OutputSelection, ResolvedOutputs};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use std::io::Write;
use treetime::ancestral;
use treetime::branch_lengths::branch_lengths_or_zero;
use treetime::cancel::Cancel;
use treetime::homoplasy;
use treetime::homoplasy::pipeline::{HomoplasyInput, HomoplasyParams};
use treetime::partition::marginal::sample::SampleMode;
use treetime::progress::{LogSink, StageSink};
use treetime::progress_info;
use treetime::seq::mutation::PositionBase;
use treetime_io::nwk::NwkParse;
use treetime_utils::io::file::write_file_with;
use treetime_utils::io::json::{JsonPretty, json_write_file};
use treetime_utils::{make_internal_error, make_internal_report};

pub fn run_homoplasy(
  args: &TreetimeHomoplasyArgs,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<HomoplasyResult, Report> {
  let resolved = args.outputs.resolve()?;
  let drms = args.drms.as_deref().map(DrmTable::read_file).transpose()?;

  let TreeAndSources {
    tree,
    sources,
    partitions,
  } = read_tree_and_sources(
    &args.tree,
    args.newick_dialect,
    &args.sequences,
    &args.model,
    cancel,
    stages,
    log,
  )?;
  log_partitions(&partitions, log);
  let [AlignmentSource { source, n_records, .. }] = <[AlignmentSource; 1]>::try_from(sources)
    .map_err(|sources| make_internal_report!("Expected one alignment source, found {}", sources.len()))?;
  progress_info!(
    log,
    "Read {} sequences of length {} and a tree with {} leaves",
    n_records,
    source.length,
    tree.graph.num_leaves()
  );
  let NwkParse {
    graph,
    names,
    branch_lengths,
    ..
  } = tree;
  let mut input = ancestral_input(graph, branch_lengths, vec![source]);
  for edge in input.edges.values_mut() {
    edge.branch_length = edge.branch_length.map(|branch_length| branch_length * args.rescale);
  }
  let branch_lengths = input.branch_lengths();
  let topology_order = args
    .topology_order
    .resolve(&input.graph, &names, None, args.newick_dialect)?;

  let random_step = (args.sample_from_profile != SampleMode::Argmax).then_some("Sampling from the profile");
  let seed = SeedArgs::resolve(args.seed, random_step, log);
  let params = homoplasy_ancestral_params(args, seed);
  let output = ancestral::pipeline::run(&params, input, None, cancel, stages, log).map_err(|err| err.into_report())?;
  let mut graph = output.graph;
  let result = output.results.into_iter().exactly_one().map_err(|error| {
    make_internal_report!("Homoplasy reconstruction is expected to return one partition result: {error}")
  })?;
  let Some(sequences) = result.sequences else {
    return make_internal_error!("Homoplasy reconstruction returned no reconstructed sequences");
  };
  let alphabet = result.alphabet;
  let partition_number = result.number;
  let raw_mutations = sequences.edge_mutations;
  let bridged_mutations =
    UnknownMutationFilter::hiding_unknown(alphabet.unknown()).reported_edge_mutations(&graph, raw_mutations.clone())?;

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
      sequence_length: result.length,
    },
  )
  .map_err(|err| err.into_report())?;
  let result = homoplasy_result(
    &statistics,
    &ResultContext {
      names: &names,
      drms: drms.as_ref(),
      zero_based: args.positions == PositionBase::Zero,
    },
  );
  let report = render_homoplasy_report(&result, args.num_mut, args.detailed);

  stages.report("Writing output", 0.9, "");
  write_homoplasy_outputs(&resolved, &result, &report)?;
  topology_order.apply(&mut graph, &names, &branch_lengths)?;
  let maps = AncestralOutputMaps {
    root_sequence: sequences.root_sequence,
    edge_mutations: bridged_mutations,
  };
  let trees = AncestralTrees {
    graph: &graph,
    names: &names,
    branch_lengths: &branch_lengths,
    amino_acids: None,
    comment_positions: args.positions,
  };
  let partition = AncestralPartitionTrees {
    number: partition_number,
    alphabet: &alphabet,
    maps,
    node_sequences: None,
  };
  write_ancestral_trees(&trees, vec![partition], &resolved, CommandKind::Homoplasy, log)?;

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
    write_file_with(path, |writer| {
      writer.write_all(report.as_bytes()).wrap_err("When writing the report")
    })?;
  }
  Ok(())
}
