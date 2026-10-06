use crate::commands::shared::alignment::sequence_descriptions;
use crate::commands::shared::gtr_output::write_gtr_output;
use crate::commands::shared::output_args::DivergenceUnits;
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use crate::commands::timetree::initialization::load_input_data;
use crate::commands::timetree::trace::TimetreeTraceSink;
use app_output::annotated_graph::{AnnotatedGraph, Divergence, TreeDates, TreeSequences};
use app_output::augur_node_data_refine::{RefineRun, write_augur_node_data_refine};
use app_output::mutation_filter::UnknownMutationFilter;
use app_output::output_plan::{CommandKind, OutputSelection, ResolvedOutputs, TreeWriteKind, output_unavailable};
use app_output::table_output::table_write_file;
use app_output::tree_output::{tree_view_for_outputs, write_graph_outputs, write_tree_outputs};
use eyre::{Report, WrapErr};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use treetime::alphabet::alphabet::Alphabet;
use treetime::cancel::Cancel;
use treetime::clock::divergence::root_to_node_divergences;
use treetime::optimize::params::BranchLengthMode;
use treetime::progress::{LogSink, StageSink};
use treetime::progress_info;
use treetime::seq::mutation::Mutation;
use treetime::seq::sink::{SeqItem, SeqSink, SeqTrack};
use treetime::timetree::coalescent::CoalescentOutput;
use treetime::timetree::confidence::NodeConfidenceInterval;
use treetime::timetree::params::TimetreeParams;
use treetime::timetree::pipeline::{self, TimetreeInput, TimetreeOutput, TimetreeSequences};
use treetime::{make_error, make_internal_error};
use treetime_graph::assign_node_names::assign_node_names;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::dates_csv::DateConstraint;
use treetime_io::fasta::FastaWriter;
use treetime_primitives::{AlignmentRecord, Seq};
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
  let parse_names = input_data.names;

  let resolved = args.resolve_outputs()?;
  let reconstructed_nuc_fasta = reconstructed_nuc_fasta_path(&resolved, args.branch_length_mode, log)?;
  let mutation_units = matches!(args.divergence_units, DivergenceUnits::Mutations);
  if mutation_units && args.branch_length_mode == BranchLengthMode::Input {
    return make_error!(
      "--divergence-units=mutations requires ancestral reconstruction; \
       incompatible with --branch-length-mode=input"
    );
  }
  let sequence_outputs_requested =
    sequence_outputs_requested(&resolved, mutation_units, reconstructed_nuc_fasta.is_some());
  let mut trace_sink = TimetreeTraceSink::new(resolved.path(OutputSelection::Tracelog), stages)?;
  let seed = args
    .seed_args
    .resolve(args.resolve_polytomies.then_some("Polytomy resolution"), log);
  let params = timetree_params(args, sequence_outputs_requested, seed);

  let aln_descs = sequence_descriptions(input_data.aln.iter().flatten());
  let alphabet = input_data.alphabet.clone();
  let unknown = alphabet.unknown();
  let input = TimetreeInput {
    graph: input_data.graph,
    names: parse_names.clone(),
    alphabet: input_data.alphabet,
    sequences: input_data
      .aln
      .map(|records| records.into_iter().map(AlignmentRecord::from).collect()),
    dates: input_data.dates,
    branch_lengths: input_data.branch_lengths,
  };

  let mut recon_sink = reconstructed_nuc_fasta
    .as_ref()
    .map(|path| {
      Ok::<_, Report>(ReconstructedNucSink::new(
        FastaWriter::create(path)?,
        parse_names,
        aln_descs,
      ))
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
  if let (Some(sink), Some(path)) = (recon_sink, &reconstructed_nuc_fasta) {
    sink.writer.finish()?;
    progress_info!(
      log,
      "Wrote reconstructed nucleotide FASTA to {path}",
      path = path.display()
    );
  }

  stages.report("Writing output", 0.95, "");
  progress_info!(log, "### TreeTime: writing outputs");
  write_model_outputs(&resolved, &output, log)?;
  let tree_inputs = TreeOutputInputs {
    alphabet,
    input_leaf_order,
    filter: UnknownMutationFilter::new(unknown, args.report_ambiguous),
    mutation_units,
  };
  write_result_outputs(args, &resolved, output, tree_inputs, log)?;

  stages.report("Done", 1.0, "");
  Ok(())
}

fn timetree_params(args: &TreetimeTimetreeArgs, sequence_outputs_requested: bool, seed: u64) -> TimetreeParams {
  TimetreeParams {
    model: args.model_args.model_name(),
    dense: args.dense,
    branch_length_mode: args.branch_length_mode,
    no_indels: args.no_indels,
    sequence_length: args.sequence_length,
    clock_rate: args.clock_rate,
    clock_std_dev: args.clock_std_dev,
    keep_root: args.keep_root,
    reroot_spec: args.reroot.spec(),
    allow_negative_rate: args.allow_negative_rate,
    clock_filter: args.clock_filter,
    covariation: args.covariation,
    tip_slack: args.tip_slack,
    max_iter: args.max_iter,
    resolve_polytomies: args.resolve_polytomies,
    relax: args.relax.clone(),
    coalescent: args.coalescent,
    coalescent_opt: args.coalescent_opt,
    coalescent_skyline: args.coalescent_skyline,
    skyline_n_points: args.skyline_n_points,
    skyline_stiffness: args.skyline_stiffness,
    coalescent_confidence: args.coalescent_confidence,
    gen_per_year: args.gen_per_year,
    n_branches_posterior: args.n_branches_posterior,
    time_marginal: args.time_marginal,
    confidence: args.confidence,
    include_leaves: args.include_leaves,
    report_ambiguous: args.report_ambiguous,
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
        table_write_file(OutputSelection::ConfidenceTsv, &file.path, intervals)
          .wrap_err("Failed to write confidence intervals")?;
        progress_info!(log, "Wrote confidence intervals to {}", file.path.display());
      },
      None => output_unavailable(
        OutputSelection::ConfidenceTsv,
        file,
        "no confidence intervals were computed (use --time-marginal)",
        log,
      )?,
    }
  }

  for selection in [
    OutputSelection::CoalescentTsv,
    OutputSelection::CoalescentCsv,
    OutputSelection::CoalescentJson,
  ] {
    let Some(file) = resolved.non_tree_outputs.get(&selection) else {
      continue;
    };
    match output.coalescent.as_ref() {
      Some(coalescent) => {
        write_coalescent(selection, coalescent, &file.path)?;
        progress_info!(log, "Wrote coalescent output to {}", file.path.display());
      },
      None => output_unavailable(
        selection,
        file,
        "no coalescent model was set (use --coalescent, --coalescent-opt, or --coalescent-skyline)",
        log,
      )?,
    }
  }

  if let Some(path) = resolved.path(OutputSelection::ClockModel) {
    json_write_file(path, &output.clock_model, JsonPretty(true))?;
  }
  if let Some(path) = resolved.path(OutputSelection::ClockCsv) {
    table_write_file(OutputSelection::ClockCsv, path, &output.clock_regression)?;
  }

  write_gtr_output(
    resolved,
    output.gtr.as_ref().zip(output.model_name),
    "no GTR model was fitted (provide an alignment with --alignment)",
    log,
  )
}

fn write_coalescent(selection: OutputSelection, coalescent: &CoalescentOutput, path: &Path) -> Result<(), Report> {
  match selection.table_format() {
    Some(_) => table_write_file(selection, path, coalescent.rows()),
    None => json_write_file(path, coalescent, JsonPretty(true)),
  }
}

pub(crate) fn reconstructed_nuc_fasta_path(
  resolved: &ResolvedOutputs,
  branch_length_mode: BranchLengthMode,
  log: &dyn LogSink,
) -> Result<Option<PathBuf>, Report> {
  let Some(file) = resolved.non_tree_outputs.get(&OutputSelection::ReconstructedNucFasta) else {
    return Ok(None);
  };
  match branch_length_mode {
    BranchLengthMode::Marginal => Ok(Some(file.path.clone())),
    BranchLengthMode::Input => {
      output_unavailable(
        OutputSelection::ReconstructedNucFasta,
        file,
        "--branch-length-mode=input reconstructs no ancestral sequences",
        log,
      )?;
      Ok(None)
    },
  }
}

struct TreeOutputInputs {
  alphabet: Alphabet,
  input_leaf_order: Vec<String>,
  filter: UnknownMutationFilter,
  mutation_units: bool,
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
    sequences,
    ..
  } = output;
  let (maps, mutation_counts) = timetree_output_maps(&graph, sequences, inputs.filter, inputs.mutation_units)?;

  let topology_order = args
    .topology_order
    .resolve_topology_order(&graph, &names, Some(inputs.input_leaf_order))?;
  topology_order.apply(&mut graph, &names, &branch_lengths)?;

  let trees = TimetreeTrees {
    graph: &graph,
    alphabet: &inputs.alphabet,
    names: &names,
    branch_lengths: &branch_lengths,
    date_branch_lengths: &date_branch_lengths,
    divergences: &divergences,
    maps: &maps,
    mutation_counts: mutation_counts.as_ref(),
    node_dates: &node_dates,
    confidence_intervals: confidence_intervals.as_deref(),
    outliers: &outliers,
    bad_branches: &bad_branches,
    input_dates: dates.as_ref(),
  };
  let run = RefineRun {
    alignment: args.alignment.alignment.first().map(PathBuf::as_path),
    input_tree: Some(args.tree.as_path()),
    clock_model: Some(&clock_model),
  };
  write_timetree_trees(&trees, &run, resolved, log)
}

struct TimetreeTrees<'a> {
  graph: &'a Graph,
  alphabet: &'a Alphabet,
  names: &'a BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &'a BTreeMap<GraphEdgeKey, Option<f64>>,
  date_branch_lengths: &'a BTreeMap<GraphEdgeKey, Option<f64>>,
  divergences: &'a BTreeMap<GraphNodeKey, f64>,
  maps: &'a TimetreeOutputMaps,
  mutation_counts: Option<&'a BTreeMap<GraphEdgeKey, usize>>,
  node_dates: &'a BTreeMap<GraphNodeKey, Option<f64>>,
  confidence_intervals: Option<&'a [NodeConfidenceInterval]>,
  outliers: &'a BTreeSet<GraphNodeKey>,
  bad_branches: &'a BTreeMap<GraphNodeKey, bool>,
  input_dates: Option<&'a BTreeMap<GraphNodeKey, DateConstraint>>,
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
    sequences: trees.maps.root_sequence.as_ref().map(|root_sequence| TreeSequences {
      alphabet: trees.alphabet,
      root_sequence,
      edge_mutations: &trees.maps.edge_mutations,
      mutation_counts: trees.mutation_counts,
      amino_acids: None,
    }),
    dates: Some(TreeDates {
      num_date: trees.node_dates,
      confidence: date_confidence.as_ref(),
      excluded: &excluded,
      input_dates: trees.input_dates,
    }),
    traits: None,
  };
  write_graph_outputs(&annotated, &resolved.tree_outputs)?;
  let Some(tree) = tree_view_for_outputs(&annotated, resolved)? else {
    return Ok(());
  };
  write_tree_outputs(&tree, &resolved.tree_outputs, CommandKind::Timetree, log)?;
  if let Some(path) = resolved.path(OutputSelection::AugurNodeData) {
    write_augur_node_data_refine(&tree, run, path)?;
    progress_info!(log, "Wrote augur node data JSON to {path}", path = path.display());
  }
  Ok(())
}

#[expect(
  clippy::as_conversions,
  reason = "a mutation count is far below 2^53, so the conversion to f64 is exact"
)]
fn mutation_divergences(
  graph: &Graph,
  mutation_counts: &BTreeMap<GraphEdgeKey, usize>,
) -> Result<BTreeMap<GraphNodeKey, f64>, Report> {
  root_to_node_divergences(graph, |edge_key| {
    mutation_counts.get(&edge_key).copied().unwrap_or_default() as f64
  })
}

struct ReconstructedNucSink {
  writer: FastaWriter,
  parse_names: BTreeMap<GraphNodeKey, Option<String>>,
  aln_descs: BTreeMap<String, Option<String>>,
  names: BTreeMap<GraphNodeKey, Option<String>>,
  descs: BTreeMap<GraphNodeKey, Option<String>>,
}

impl ReconstructedNucSink {
  fn new(
    writer: FastaWriter,
    parse_names: BTreeMap<GraphNodeKey, Option<String>>,
    aln_descs: BTreeMap<String, Option<String>>,
  ) -> Self {
    Self {
      writer,
      parse_names,
      aln_descs,
      names: BTreeMap::new(),
      descs: BTreeMap::new(),
    }
  }
}

impl SeqSink for ReconstructedNucSink {
  fn on_topology(&mut self, graph: &Graph) -> Result<(), Report> {
    self.names = assign_node_names(self.parse_names.clone(), graph)?.names;
    self.descs = self
      .names
      .iter()
      .map(|(&key, name)| {
        let desc = name
          .as_deref()
          .and_then(|name| self.aln_descs.get(name).cloned())
          .flatten();
        (key, desc)
      })
      .collect();
    Ok(())
  }

  fn emit(&mut self, item: SeqItem<'_>) -> Result<(), Report> {
    match item.track {
      SeqTrack::Nuc if !item.emitted => Ok(()),
      SeqTrack::Nuc => {
        let name = self.names[&item.key].clone().unwrap_or_default();
        let desc = self.descs[&item.key].clone();
        self.writer.write(&name, desc.as_deref(), item.seq)
      },
      SeqTrack::Aa(cds) => {
        make_internal_error!("Timetree reconstructed-nucleotide FASTA sink received an amino-acid track '{cds}'")
      },
    }
  }
}

fn timetree_output_maps(
  graph: &Graph,
  sequences: Option<TimetreeSequences>,
  filter: UnknownMutationFilter,
  mutation_units: bool,
) -> Result<(TimetreeOutputMaps, Option<BTreeMap<GraphEdgeKey, usize>>), Report> {
  let Some(TimetreeSequences {
    root_sequence,
    edge_mutations,
    edge_mutation_counts,
  }) = sequences
  else {
    if mutation_units {
      return make_internal_error!("Mutation divergence units were requested, but no ancestral reconstruction ran");
    }
    let maps = TimetreeOutputMaps {
      root_sequence: None,
      edge_mutations: BTreeMap::new(),
    };
    return Ok((maps, None));
  };
  let maps = TimetreeOutputMaps {
    root_sequence: Some(root_sequence),
    edge_mutations: filter.reported_edge_mutations(graph, edge_mutations)?,
  };
  Ok((maps, mutation_units.then_some(edge_mutation_counts)))
}

fn sequence_outputs_requested(resolved: &ResolvedOutputs, mutation_units: bool, writes_fasta: bool) -> bool {
  mutation_units
    || writes_fasta
    || resolved
      .tree_outputs
      .keys()
      .any(|kind| !matches!(kind, TreeWriteKind::GraphJson | TreeWriteKind::Dot))
}

struct TimetreeOutputMaps {
  root_sequence: Option<Seq>,
  edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
}
