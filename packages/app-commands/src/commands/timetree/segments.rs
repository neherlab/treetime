use crate::commands::shared::dates_input::read_input_dates;
use crate::commands::shared::gtr_output::{partition_models, write_gtr_output};
use crate::commands::shared::partition::log_partitions;
use crate::commands::shared::segment_outputs::{segment_graphs, write_whole_graph_outputs};
use crate::commands::shared::sequence_inputs::{TreeAndSources, read_tree_input_and_sources};
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use app_output::annotated_graph::{AnnotatedGraph, Divergence, TreeDates, TreeSequences};
use app_output::command_outputs::CommandKind;
use app_output::dates_table::write_graph_dates_table_file;
use app_output::mutation_filter::UnknownMutationFilter;
use app_output::output_plan::{OutputSelection, ResolvedOutputs, output_unavailable};
use app_output::tree_output::write_tree_outputs_per_alignment;
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;
use treetime::cancel::Cancel;
use treetime::partition::split::SequenceInput;
use treetime::progress::{LogSink, StageSink};
use treetime::seq::mutation::Mutation;
use treetime::timetree::segments::{
  SegmentTimetreeInput, SegmentTimetreeOutput, SegmentTimetreeParams, run_segment_timetree,
};
use treetime::{progress_info, progress_warn};
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::node::GraphNodeKey;
use treetime_graph::scope::EdgeScopes;
use treetime_io::nwk::NwkParse;
use treetime_primitives::Seq;
use treetime_utils::make_internal_error;

const UNAVAILABLE: [OutputSelection; 6] = [
  OutputSelection::NodeData,
  OutputSelection::ClockModel,
  OutputSelection::Rtt,
  OutputSelection::Tracelog,
  OutputSelection::Coalescent,
  OutputSelection::CoalescentJson,
];

pub(crate) fn run_segment_timetree_command(
  args: &TreetimeTimetreeArgs,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let resolved = args.outputs.resolve()?;
  let Some(clock_rate) = args.clock.clock_rate else {
    return make_internal_error!("A timetree of segment trees needs a fixed clock rate");
  };
  progress_info!(
    log,
    "A graph of segment trees keeps its roots and polytomies and runs without the clock filter"
  );
  let (
    TreeAndSources {
      tree,
      sources,
      partitions,
    },
    scopes,
  ) = read_tree_input_and_sources(
    &args.tree_input,
    args.newick_dialect,
    &args.sequences,
    &args.sequence_model.model,
    cancel,
    stages,
    log,
  )?;
  let Some(scopes) = scopes else {
    return make_internal_error!("Segment trees were read without their edge scopes");
  };
  let NwkParse {
    graph,
    names,
    branch_lengths,
    ..
  } = tree;
  log_partitions(&partitions, log);
  let dates = read_input_dates(&args.dates, &graph, &names, log)?;
  let params = SegmentTimetreeParams {
    dense: args.sequence_model.dense,
    indel_policy: args.branch_lengths.indel_policy,
    partition_rates: args.sequence_model.partition_rates,
    clock_rate,
    confidence: args.dating.confidence,
    max_grid_points: args.dating.max_grid_points,
    cross_parent_messages: args.cross_parent_messages,
    refinement_rounds: args.dating.max_iter,
  };
  let input = SegmentTimetreeInput {
    graph,
    names,
    sequences: SequenceInput {
      alignments: sources.into_iter().map(|source| source.source).collect(),
      partitions,
    },
    dates,
    branch_lengths,
    scopes: Arc::clone(&scopes),
  };
  let output = run_segment_timetree(&params, input, cancel, stages, log).map_err(|err| err.into_report())?;

  stages.report("Writing output", 0.95, "");
  write_gtr_output(
    &resolved,
    &partition_models(&output.partitions, &output.alignments, Some(clock_rate))?,
    "no GTR model was fitted",
    log,
  )?;
  for selection in UNAVAILABLE {
    if let Some(file) = resolved.non_tree_outputs.get(&selection) {
      output_unavailable(
        selection,
        file,
        "a timetree of segment trees does not write it yet",
        log,
      )?;
    }
  }
  if !resolved.alignment_paths(OutputSelection::Sequences).is_empty() {
    progress_warn!(
      log,
      "Not writing reconstructed sequences: a timetree of segment trees does not write them yet"
    );
  }
  write_segment_outputs(args, &resolved, output, &scopes, log)?;
  stages.report("Done", 1.0, "");
  Ok(())
}

fn write_segment_outputs(
  args: &TreetimeTimetreeArgs,
  resolved: &ResolvedOutputs,
  output: SegmentTimetreeOutput,
  scopes: &Arc<EdgeScopes>,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let SegmentTimetreeOutput {
    graph,
    names,
    branch_lengths,
    node_dates,
    date_branch_lengths,
    confidence_intervals,
    alignments: mut results,
    ..
  } = output;
  let segment_graphs = segment_graphs(&graph, Some(scopes))?;
  let maps = results
    .iter_mut()
    .map(|result| {
      let Some(sequences) = result.sequences.take() else {
        return make_internal_error!("The timetree of segment trees returned no reconstructed sequences");
      };
      let filter = UnknownMutationFilter::new(result.alphabet.unknown(), args.tip_sequences.report_ambiguous);
      let result_graph = segment_graphs
        .as_ref()
        .and_then(|graphs| graphs.get(&result.number))
        .unwrap_or(&graph);
      Ok((
        result.number,
        SegmentMaps {
          root_sequence: sequences.root_sequence,
          edge_mutations: filter.reported_edge_mutations(result_graph, sequences.edge_mutations)?,
        },
      ))
    })
    .collect::<Result<BTreeMap<_, _>, Report>>()?;
  let sequences = results
    .iter()
    .map(|result| {
      let maps = &maps[&result.number];
      let sequences = TreeSequences {
        alphabet: &result.alphabet,
        root_sequence: &maps.root_sequence,
        edge_mutations: &maps.edge_mutations,
        amino_acids: None,
        comment_positions: args.positions,
      };
      (result.number, sequences)
    })
    .collect();
  let confidence: Option<NodeMap<[f64; 2]>> = confidence_intervals.map(|intervals| {
    intervals
      .iter()
      .map(|interval| (interval.key, [interval.lower, interval.upper]))
      .collect()
  });
  let excluded = BTreeSet::<GraphNodeKey>::new();
  let annotated = AnnotatedGraph {
    graph: &graph,
    names: &names,
    divergence_branch_lengths: &branch_lengths,
    time_branch_lengths: Some(&date_branch_lengths),
    divergence: Divergence::CumulativeBranchLength,
    sequences: None,
    dates: Some(TreeDates {
      num_date: &node_dates,
      confidence: confidence.as_ref(),
      excluded: &excluded,
      input_dates: None,
    }),
    traits: None,
  };
  write_whole_graph_outputs(&annotated, resolved, &[], log)?;
  if let Some(path) = resolved.path(OutputSelection::Dates) {
    write_graph_dates_table_file(&annotated, path)?;
    progress_info!(log, "Wrote node dates to {path}", path = path.display());
  }
  write_tree_outputs_per_alignment(
    &annotated,
    &sequences,
    segment_graphs.as_ref(),
    resolved,
    CommandKind::Timetree,
    log,
  )
}

struct SegmentMaps {
  root_sequence: Seq,
  edge_mutations: EdgeMap<Vec<Mutation>>,
}
