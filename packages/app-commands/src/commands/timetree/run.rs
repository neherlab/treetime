use crate::commands::shared::alignment::sequence_descriptions;
use crate::commands::shared::output_args::DivergenceUnits;
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use crate::commands::timetree::initialization::load_input_data;
use crate::commands::timetree::trace::TimetreeTraceSink;
use app_output::DateCommentProvider;
use app_output::EdgeMutationCommentProvider;
use app_output::augur_node_data::write_augur_node_data_json;
use app_output::clock_model::write_clock_model;
use app_output::coalescent::{write_coalescent_delimited, write_coalescent_json};
use app_output::confidence::write_confidence_intervals_file;
use app_output::gtr::write_gtr_json;
use app_output::mutation_filter::UnknownMutationFilter;
use app_output::output_plan::OutputSelection;
use app_output::rtt::write_clock_regression_result_csv;
use app_output::timetree_tree_output::write_timetree_tree_outputs;
use app_output::{TimetreeEdgeOut, TimetreeNodeOut, TimetreeOutputMaps};
use eyre::{Report, WrapErr};
use log::debug;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use treetime::cancel::Cancel;
use treetime::gtr::get_gtr::GtrOutput;
use treetime::make_error;
use treetime::optimize::params::BranchLengthMode;
use treetime::progress::{LogSink, StageSink};
use treetime::seq::sink::{SeqItem, SeqSink, SeqTrack};
use treetime::timetree::coalescent::CoalescentOutput;
use treetime::timetree::params::TimetreeParams;
use treetime::timetree::pipeline::{self, TimetreeInput};
use treetime::{progress_info, progress_warn};
use treetime_graph::assign_node_names::assign_node_names;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::FastaWriter;
use treetime_io::graph::TreeWriteKind;
use treetime_io::nwk::CommentProviders;
use treetime_primitives::AlignmentRecord;
use treetime_utils::io::file::create_file_or_stdout;

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
  let confidences = input_data.confidences;
  let parse_names = input_data.names;

  let resolved = args.resolve_outputs()?;
  let reconstructed_nuc_fasta = resolved
    .non_tree_outputs
    .get(&OutputSelection::ReconstructedNucFasta)
    .cloned();
  if reconstructed_nuc_fasta.is_some() && args.branch_length_mode == BranchLengthMode::Input {
    return make_error!(
      "Reconstructed sequence output requires ancestral reconstruction; \
       incompatible with --branch-length-mode=input"
    );
  }
  let mutation_units = matches!(args.divergence_units, DivergenceUnits::Mutations);
  if mutation_units && args.branch_length_mode == BranchLengthMode::Input {
    return make_error!(
      "--divergence-units=mutations requires ancestral reconstruction; \
       incompatible with --branch-length-mode=input"
    );
  }
  let sequence_outputs_requested = reconstructed_nuc_fasta.is_some()
    || mutation_units
    || resolved
      .tree_outputs
      .keys()
      .any(|kind| !matches!(kind, TreeWriteKind::GraphJson | TreeWriteKind::Dot));
  let mut trace_sink = TimetreeTraceSink::new(
    resolved
      .non_tree_outputs
      .get(&OutputSelection::Tracelog)
      .map(PathBuf::as_path),
    stages,
  )?;

  let params = TimetreeParams {
    model: args.model_args.model_name(),
    alphabet_name: args.alphabet_args.alphabet_name().unwrap_or_default(),
    dense: args.dense,
    gap_fill: args.gap_fill_args.effective_gap_fill(),
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
    keep_polytomies: args.keep_polytomies,
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
    impute_missing_data: args.impute_missing_data,
    sequence_outputs_requested,
    zero_based: args.zero_based,
    seed: args.seed,
  };

  let aln_descs = sequence_descriptions(input_data.aln.iter().flatten());
  let unknown = input_data.alphabet.unknown();

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
        FastaWriter::new(create_file_or_stdout(path)?),
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

  let pipeline::TimetreeOutput {
    mut graph,
    names,
    node_dates,
    divergences,
    outliers,
    bad_branches,
    branch_lengths,
    date_branch_lengths,
    clock_model,
    clock_regression,
    confidence_intervals,
    gtr,
    model_name,
    coalescent,
    dates,
    root_sequence,
    edge_mutations,
    edge_mutation_counts,
  } = output;

  let mutation_counts = mutation_units.then_some(edge_mutation_counts);
  let maps = TimetreeOutputMaps {
    root_sequence,
    edge_mutations: UnknownMutationFilter::new(unknown, args.report_ambiguous)
      .reported_edge_mutations(&graph, edge_mutations)?,
  };

  stages.report("Writing output", 0.95, "");
  progress_info!(log, "### TreeTime: writing outputs");

  let topology_order = args
    .topology_order
    .resolve_topology_order(&graph, &names, Some(input_leaf_order))?;
  topology_order.apply(&mut graph, &names, &branch_lengths)?;

  let nodes = timetree_node_outputs(
    &graph,
    &names,
    &confidences,
    &node_dates,
    &divergences,
    &outliers,
    &bad_branches,
  );
  let edges = timetree_edge_outputs(&branch_lengths, &date_branch_lengths);

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::ConfidenceTsv) {
    match confidence_intervals.as_ref() {
      Some(intervals) => {
        write_confidence_intervals_file(intervals, path).wrap_err("Failed to write confidence intervals")?;
        progress_info!(log, "Wrote confidence intervals to {path}", path = path.display());
      },
      None if args.output_confidence_tsv.is_some() => {
        return make_error!(
          "Confidence output requested but no confidence intervals were computed. \
           Use --time-marginal to enable confidence interval computation."
        );
      },
      None => progress_warn!(
        log,
        "Skipping confidence-interval output: no confidence intervals were computed (use --time-marginal)"
      ),
    }
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::CoalescentTsv) {
    write_coalescent_output(
      coalescent.as_ref(),
      path,
      args.output_coalescent_tsv.is_some(),
      |output, path| write_coalescent_delimited(output, path, b'\t'),
      log,
    )?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::CoalescentCsv) {
    write_coalescent_output(
      coalescent.as_ref(),
      path,
      args.output_coalescent_csv.is_some(),
      |output, path| write_coalescent_delimited(output, path, b','),
      log,
    )?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::CoalescentJson) {
    write_coalescent_output(
      coalescent.as_ref(),
      path,
      args.output_coalescent_json.is_some(),
      |output, path| write_coalescent_json(output, path),
      log,
    )?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::ClockModel) {
    write_clock_model(&clock_model, path)?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::ClockCsv) {
    write_clock_regression_result_csv(&clock_regression, path, b',')?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) {
    match (gtr.as_ref(), model_name) {
      (Some(gtr), Some(model_name)) => {
        let gtr_output = GtrOutput::builder().gtr(gtr).model_name(model_name).build();
        write_gtr_json(&gtr_output, path)?;
      },
      _ if args.output_gtr.is_some() => {
        return make_error!("GTR output requested but no GTR model was fitted. Provide sequence alignment input.");
      },
      _ => progress_warn!(
        log,
        "Skipping GTR output: no GTR model was fitted (provide sequence alignment input)"
      ),
    }
  }

  if !resolved.tree_outputs.is_empty() {
    let date_times: BTreeMap<GraphNodeKey, f64> = nodes
      .iter()
      .filter_map(|(key, out)| out.time.map(|time| (*key, time)))
      .collect();
    let date_provider = DateCommentProvider::new(&date_times);
    let mutation_provider = maps
      .root_sequence
      .is_some()
      .then(|| EdgeMutationCommentProvider::new(&maps.edge_mutations, &graph));
    let providers = mutation_provider
      .iter()
      .fold(CommentProviders::new(), |providers, provider| providers.with(provider))
      .with(&date_provider);
    write_timetree_tree_outputs(
      &graph,
      &nodes,
      &edges,
      &maps,
      confidence_intervals.as_deref(),
      mutation_counts.as_ref(),
      &resolved.tree_outputs,
      &providers,
    )?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::AugurNodeData) {
    let alignment = args.alignment.alignment.first().map(PathBuf::as_path);
    write_augur_node_data_json(
      &graph,
      &nodes,
      &edges,
      &clock_model,
      confidence_intervals.as_deref(),
      dates.as_ref(),
      alignment,
      Some(args.tree.as_path()),
      mutation_counts.as_ref(),
      path,
    )?;
    progress_info!(log, "Wrote augur node data JSON to {path}", path = path.display());
  }

  stages.report("Done", 1.0, "");
  Ok(())
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
    self.names = assign_node_names(self.parse_names.clone(), graph)?;
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
        self.writer.write(&name, &desc, item.seq)
      },
      SeqTrack::Aa(cds) => treetime_utils::make_internal_error!(
        "Timetree reconstructed-nucleotide FASTA sink received an amino-acid track '{cds}'"
      ),
    }
  }
}

fn timetree_node_outputs(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  confidences: &BTreeMap<GraphNodeKey, Option<f64>>,
  node_dates: &BTreeMap<GraphNodeKey, Option<f64>>,
  divergences: &BTreeMap<GraphNodeKey, f64>,
  outliers: &BTreeSet<GraphNodeKey>,
  bad_branches: &BTreeMap<GraphNodeKey, bool>,
) -> BTreeMap<GraphNodeKey, TimetreeNodeOut> {
  graph
    .get_nodes()
    .map(|node| {
      let key = node.key();
      let out = TimetreeNodeOut {
        name: names[&key].clone(),
        confidence: confidences.get(&key).copied().flatten(),
        time: node_dates[&key],
        div: divergences[&key],
        is_outlier: outliers.contains(&key),
        bad_branch: bad_branches[&key],
      };
      (key, out)
    })
    .collect()
}

fn timetree_edge_outputs(
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  date_branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> BTreeMap<GraphEdgeKey, TimetreeEdgeOut> {
  date_branch_lengths
    .iter()
    .map(|(&key, &date_branch_length)| {
      let out = TimetreeEdgeOut {
        branch_length: branch_lengths[&key],
        date_branch_length,
      };
      (key, out)
    })
    .collect()
}

fn write_coalescent_output(
  coalescent: Option<&CoalescentOutput>,
  path: &Path,
  explicit: bool,
  write: impl FnOnce(&CoalescentOutput, &Path) -> Result<(), Report>,
  log: &dyn LogSink,
) -> Result<(), Report> {
  match coalescent {
    Some(output) => {
      write(output, path)?;
      progress_info!(log, "Wrote coalescent output to {path}", path = path.display());
    },
    None if explicit => {
      return make_error!(
        "Coalescent output requested but no coalescent model was set. \
         Use --coalescent, --coalescent-opt, or --coalescent-skyline."
      );
    },
    None => debug!(
      "Skipping coalescent output: no coalescent model was set \
       (use --coalescent, --coalescent-opt, or --coalescent-skyline)"
    ),
  }
  Ok(())
}
