use crate::commands::shared::gtr_output::{partition_models, write_gtr_output};
use crate::commands::shared::output_args::DivergenceUnits;
use crate::commands::shared::output_args::ResolveOutputs;
use crate::commands::shared::reconstructed_fasta::NucFastaWriters;
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use crate::commands::timetree::branch_length_mode::settle_branch_length_outputs;
use crate::commands::timetree::coalescent::CoalescentArg;
use crate::commands::timetree::initialization::load_input_data;
use crate::commands::timetree::trace::TimetreeTraceSink;
use app_output::annotated_graph::{AnnotatedGraph, Divergence, TreeDates, TreeSequences};
use app_output::augur_node_data_refine::{RefineRun, total_mutation_counts, write_augur_node_data_refine};
use app_output::command_outputs::CommandKind;
use app_output::mutation_filter::UnknownMutationFilter;
use app_output::output_format::FileFormat;
use app_output::output_plan::{OutputSelection, ResolvedOutputs, TreeWriteKind, output_unavailable};
use app_output::table_output::table_write_file;
use app_output::tree_output::{write_run_tree_outputs, write_tree_outputs_per_partition};
use eyre::{Report, WrapErr};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;
use treetime::alphabet::alphabet::Alphabet;
use treetime::cancel::Cancel;
use treetime::clock::divergence::root_to_node_divergences;
use treetime::clock::find_best_root::params::RerootSpec;
use treetime::make_internal_error;
use treetime::partition::marginal::posterior_states::AmbiguityThreshold;
use treetime::partition::result::PartitionSequences;
use treetime::progress::{LogSink, StageSink};
use treetime::progress_info;
use treetime::seq::mutation::Mutation;
use treetime::seq::sink::{SeqItem, SeqSink, SeqTrack};
use treetime::timetree::coalescent::CoalescentOutput;
use treetime::timetree::coalescent_timescale::CoalescentMode;
use treetime::timetree::confidence::NodeConfidenceInterval;
use treetime::timetree::params::{BranchLengthMode, TimetreeParams};
use treetime::timetree::pipeline::{self, TimetreeInput, TimetreeOutput};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_io::dates_csv::DateConstraint;
use treetime_primitives::Seq;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::io::json::{JsonPretty, json_write_file};

pub fn run_timetree_estimation(
  args: &TreetimeTimetreeArgs,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<(), Report> {
  cancel.check()?;
  stages.report("Loading input", 0.0, "");

  let input_data = load_input_data(args, log)?;
  let input_leaf_order = input_data.input_leaf_order.clone();
  let branch_length_mode = input_data.branch_length_mode.mode;

  let mut resolved = args.resolve_outputs()?;
  settle_branch_length_outputs(&mut resolved, &input_data.branch_length_mode, log)?;
  let ambiguity = &args.ambiguity_outputs;
  ambiguity.warn_unwritten(&resolved, log);
  let reconstructed_nuc_fasta = resolved.partition_paths(OutputSelection::Sequences);
  let mutation_units = matches!(args.branch_lengths.divergence_units, DivergenceUnits::Mutations);
  if mutation_units {
    input_data
      .branch_length_mode
      .require_marginal(&["divergence_units"], "--divergence-units mutations")?;
  }
  let sequence_outputs_requested =
    sequence_outputs_requested(&resolved, mutation_units, !reconstructed_nuc_fasta.is_empty());
  let mut trace_sink = TimetreeTraceSink::new(resolved.path(OutputSelection::Tracelog), stages)?;
  let seed = args
    .seed_args
    .resolve((!args.keep_polytomies).then_some("Polytomy resolution"), log);
  let reroot_spec = args.reroot.spec(&args.rooting, &input_data.graph, &input_data.names)?;
  let params = timetree_params(args, branch_length_mode, reroot_spec, sequence_outputs_requested, seed);

  let (sources, aln_descs): (Vec<_>, BTreeMap<_, _>) = input_data
    .sequences
    .into_iter()
    .enumerate()
    .map(|(index, aligned)| (aligned.source, (index + 1, aligned.descs)))
    .unzip();
  let input = TimetreeInput {
    graph: input_data.graph,
    names: input_data.names,
    sources,
    dates: input_data.dates,
    branch_lengths: input_data.branch_lengths,
  };

  let mut recon_sink = (!reconstructed_nuc_fasta.is_empty())
    .then(|| {
      Ok::<_, Report>(ReconstructedNucSink {
        writers: NucFastaWriters::create(reconstructed_nuc_fasta, aln_descs)?,
        ambiguity: ambiguity.threshold_for(OutputSelection::Sequences),
      })
    })
    .transpose()?;

  let output = pipeline::run(
    &params,
    input,
    Some(&mut trace_sink),
    recon_sink.as_mut().map(|sink| -> &mut dyn SeqSink { sink }),
    cancel,
    stages,
    log,
  )
  .map_err(|err| err.into_report())?;
  trace_sink.finish()?;
  if let Some(sink) = recon_sink {
    for path in sink.writers.finish()? {
      progress_info!(
        log,
        "Wrote reconstructed nucleotide FASTA to {path}",
        path = path.display()
      );
    }
  }

  stages.report("Writing output", 0.95, "");
  progress_info!(log, "### TreeTime: writing outputs");
  write_model_outputs(&resolved, &output, log)?;
  let tree_inputs = TreeOutputInputs {
    input_leaf_order,
    report_ambiguous: args.tip_sequences.report_ambiguous,
    mutation_units,
  };
  write_result_outputs(args, &resolved, output, tree_inputs, log)?;

  stages.report("Done", 1.0, "");
  Ok(())
}

fn timetree_params(
  args: &TreetimeTimetreeArgs,
  branch_length_mode: BranchLengthMode,
  reroot_spec: RerootSpec,
  sequence_outputs_requested: bool,
  seed: u64,
) -> TimetreeParams {
  TimetreeParams {
    model: args.model.clone(),
    alphabet: args.alphabet_args.alphabet_name(),
    dense: args.representation.dense(),
    branch_length_mode,
    indel_policy: args.indel_policy,
    partition_rates: args.partition_rates_args.partition_rates(),
    sequence_length: args.sequence_length,
    clock_rate: args.clock_rate,
    clock_std_dev: args.clock_std_dev,
    keep_root: args.rooting.keep_root,
    reroot_spec,
    allow_negative_rate: args.clock_filter.allow_negative_rate,
    clock_filter: args.clock_filter.clock_filter,
    max_iter: args.max_iter,
    resolve_polytomies: !args.keep_polytomies,
    relax: args.relax.clone(),
    coalescent: args.coalescent.map_or(CoalescentMode::Disabled, CoalescentArg::mode),
    n_skyline: args.n_skyline(),
    skyline_stiffness: args.skyline_stiffness(),
    coalescent_confidence: args.coalescent_confidence(),
    gen_per_year: args.gen_per_year(),
    confidence: args.confidence,
    include_leaves: args.include_leaves,
    report_ambiguous: args.tip_sequences.report_ambiguous,
    impute_missing_data: args.impute_missing_data,
    sequence_outputs_requested,
    seed,
    max_grid_points: args.max_grid_points,
  }
}

fn write_model_outputs(resolved: &ResolvedOutputs, output: &TimetreeOutput, log: &dyn LogSink) -> Result<(), Report> {
  if let Some(file) = resolved.non_tree_outputs.get(&OutputSelection::ConfidenceTsv) {
    match output.confidence_intervals.as_ref() {
      Some(intervals) => {
        table_write_file(&file.path, intervals).wrap_err("Failed to write confidence intervals")?;
        progress_info!(log, "Wrote confidence intervals to {}", file.path.display());
      },
      None => output_unavailable(
        OutputSelection::ConfidenceTsv,
        file,
        "no confidence intervals were computed (use --confidence)",
        log,
      )?,
    }
  }

  for selection in [OutputSelection::Coalescent, OutputSelection::CoalescentJson] {
    let Some(file) = resolved.non_tree_outputs.get(&selection) else {
      continue;
    };
    match output.coalescent.as_ref() {
      Some(coalescent) => {
        write_coalescent(selection, coalescent, &file.path)?;
        progress_info!(log, "Wrote coalescent output to {}", file.path.display());
      },
      None => output_unavailable(selection, file, "no coalescent prior was set (use --coalescent)", log)?,
    }
  }

  if let Some(path) = resolved.path(OutputSelection::ClockModel) {
    json_write_file(path, &output.clock_model, JsonPretty(true))?;
  }
  if let Some(path) = resolved.path(OutputSelection::Rtt) {
    table_write_file(path, &output.clock_regression)?;
  }

  write_gtr_output(
    resolved,
    partition_models(&output.results, Some(output.clock_model.clock_rate())),
    "no GTR model was fitted (provide an alignment with --alignment)",
    log,
  )
}

fn write_coalescent(selection: OutputSelection, coalescent: &CoalescentOutput, path: &Path) -> Result<(), Report> {
  match selection.file_format() {
    FileFormat::Table => table_write_file(path, coalescent.rows()),
    FileFormat::Chart | FileFormat::Fixed => json_write_file(path, coalescent, JsonPretty(true)),
  }
}

fn write_result_outputs(
  args: &TreetimeTimetreeArgs,
  resolved: &ResolvedOutputs,
  output: TimetreeOutput,
  inputs: TreeOutputInputs,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let TimetreeOutput {
    mut graph,
    names,
    node_dates,
    divergences,
    outliers,
    bad_branches,
    branch_lengths,
    date_branch_lengths,
    clock_model,
    confidence_intervals,
    dates,
    mut results,
    ..
  } = output;
  let maps = results
    .iter_mut()
    .filter_map(|result| {
      let sequences = result.sequences.take()?;
      Some(
        timetree_output_maps(&graph, sequences, inputs.report_ambiguous, &result.alphabet)
          .map(|maps| (result.number, maps)),
      )
    })
    .collect::<Result<BTreeMap<_, _>, Report>>()?;
  if inputs.mutation_units && maps.is_empty() {
    return make_internal_error!("Mutation divergence units were requested, but no ancestral reconstruction ran");
  }
  let mutation_counts = inputs
    .mutation_units
    .then(|| total_mutation_counts(maps.values().map(|maps| &maps.edge_mutation_counts)));

  let topology_order = args.topology_order.resolve_topology_order(
    &graph,
    &names,
    Some(inputs.input_leaf_order),
    args.newick_dialect.dialect(),
  )?;
  topology_order.apply(&mut graph, &names, &branch_lengths)?;

  let sequences = results
    .iter()
    .filter_map(|result| {
      let maps = maps.get(&result.number)?;
      let sequences = TreeSequences {
        alphabet: &result.alphabet,
        root_sequence: &maps.root_sequence,
        edge_mutations: &maps.edge_mutations,
        amino_acids: None,
        comment_positions: args.zero_based.positions(),
      };
      Some((result.number, sequences))
    })
    .collect();
  let trees = TimetreeTrees {
    graph: &graph,
    names: &names,
    branch_lengths: &branch_lengths,
    date_branch_lengths: &date_branch_lengths,
    divergences: &divergences,
    sequences: &sequences,
    mutation_counts: mutation_counts.as_ref(),
    node_dates: &node_dates,
    confidence_intervals: confidence_intervals.as_deref(),
    outliers: &outliers,
    bad_branches: &bad_branches,
    input_dates: dates.as_ref(),
  };
  let run = RefineRun {
    alignments: &args.alignment.alignment,
    mutation_counts: mutation_counts.as_ref(),
    input_tree: Some(args.tree.as_path()),
    clock_model: Some(&clock_model),
  };
  write_timetree_trees(&trees, &run, resolved, log)
}

struct TreeOutputInputs {
  input_leaf_order: Vec<GraphNodeKey>,
  report_ambiguous: bool,
  mutation_units: bool,
}

fn write_timetree_trees(
  trees: &TimetreeTrees<'_>,
  run: &RefineRun<'_>,
  resolved: &ResolvedOutputs,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let mutation_divergences = trees
    .mutation_counts
    .map(|counts| mutation_divergences(trees.graph, counts))
    .transpose()?;
  let date_confidence: Option<BTreeMap<GraphNodeKey, [f64; 2]>> = trees.confidence_intervals.map(|intervals| {
    intervals
      .iter()
      .map(|interval| (interval.key, [interval.lower, interval.upper]))
      .collect()
  });
  let excluded: BTreeSet<GraphNodeKey> = trees
    .graph
    .get_nodes()
    .map(|node| node.key())
    .filter(|key| trees.outliers.contains(key) || trees.bad_branches[key])
    .collect();
  let annotated = AnnotatedGraph {
    graph: trees.graph,
    names: trees.names,
    divergence_branch_lengths: trees.branch_lengths,
    time_branch_lengths: Some(trees.date_branch_lengths),
    divergence: Divergence::Values(mutation_divergences.as_ref().unwrap_or(trees.divergences)),
    sequences: None,
    dates: Some(TreeDates {
      num_date: trees.node_dates,
      confidence: date_confidence.as_ref(),
      excluded: &excluded,
      input_dates: trees.input_dates,
    }),
    traits: None,
  };
  let tree = write_run_tree_outputs(&annotated, resolved, CommandKind::Timetree, log)?;
  if let (Some(tree), Some(path)) = (&tree, resolved.path(OutputSelection::NodeData)) {
    write_augur_node_data_refine(tree, run, path)?;
    progress_info!(log, "Wrote augur node data JSON to {path}", path = path.display());
  }
  write_tree_outputs_per_partition(&annotated, trees.sequences, resolved, CommandKind::Timetree, log)
}

struct TimetreeTrees<'a> {
  graph: &'a Graph,
  names: &'a NodeNames,
  branch_lengths: &'a BTreeMap<GraphEdgeKey, Option<f64>>,
  date_branch_lengths: &'a BTreeMap<GraphEdgeKey, Option<f64>>,
  divergences: &'a BTreeMap<GraphNodeKey, f64>,
  sequences: &'a BTreeMap<usize, TreeSequences<'a>>,
  mutation_counts: Option<&'a BTreeMap<GraphEdgeKey, usize>>,
  node_dates: &'a BTreeMap<GraphNodeKey, Option<f64>>,
  confidence_intervals: Option<&'a [NodeConfidenceInterval]>,
  outliers: &'a BTreeSet<GraphNodeKey>,
  bad_branches: &'a BTreeMap<GraphNodeKey, bool>,
  input_dates: Option<&'a BTreeMap<GraphNodeKey, DateConstraint>>,
}

fn mutation_divergences(
  graph: &Graph,
  mutation_counts: &BTreeMap<GraphEdgeKey, usize>,
) -> Result<BTreeMap<GraphNodeKey, f64>, Report> {
  root_to_node_divergences(graph, |edge_key| {
    mutation_counts
      .get(&edge_key)
      .copied()
      .unwrap_or_default()
      .to_f64_exact()
  })
}

struct ReconstructedNucSink {
  writers: NucFastaWriters,
  ambiguity: Option<AmbiguityThreshold>,
}

impl SeqSink for ReconstructedNucSink {
  fn emit(&mut self, item: SeqItem<'_>) -> Result<(), Report> {
    match item.track {
      SeqTrack::Nuc(_) if !item.emitted => Ok(()),
      SeqTrack::Nuc(partition) => self
        .writers
        .write(partition, item.key, item.name, &item.seq_with_ambiguity()),
      SeqTrack::Aa(cds) => {
        make_internal_error!("Timetree reconstructed-nucleotide FASTA sink received an amino-acid track '{cds}'")
      },
    }
  }

  fn ambiguity(&self) -> Option<AmbiguityThreshold> {
    self.ambiguity
  }
}

fn timetree_output_maps(
  graph: &Graph,
  sequences: PartitionSequences,
  report_ambiguous: bool,
  alphabet: &Alphabet,
) -> Result<TimetreeOutputMaps, Report> {
  let PartitionSequences {
    root_sequence,
    edge_mutations,
    edge_mutation_counts,
  } = sequences;
  Ok(TimetreeOutputMaps {
    root_sequence,
    edge_mutations: UnknownMutationFilter::new(alphabet.unknown(), report_ambiguous)
      .reported_edge_mutations(graph, edge_mutations)?,
    edge_mutation_counts,
  })
}

fn sequence_outputs_requested(resolved: &ResolvedOutputs, mutation_units: bool, writes_fasta: bool) -> bool {
  mutation_units
    || writes_fasta
    || resolved
      .tree_write_kinds()
      .any(|kind| !matches!(kind, TreeWriteKind::GraphJson | TreeWriteKind::Dot))
}

struct TimetreeOutputMaps {
  root_sequence: Seq,
  edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
  edge_mutation_counts: BTreeMap<GraphEdgeKey, usize>,
}
