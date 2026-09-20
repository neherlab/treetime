use crate::commands::shared::output_args::DivergenceUnits;
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use crate::commands::timetree::initialization::load_input_data;
use app_output::DateCommentProvider;
use app_output::EdgeMutationCommentProvider;
use app_output::augur_node_data::write_augur_node_data_json;
use app_output::coalescent::{write_coalescent_delimited, write_coalescent_json};
use app_output::confidence::write_confidence_intervals_file;
use app_output::output_plan::OutputSelection;
use app_output::timetree_trace::TraceCsvSink;
use app_output::timetree_tree_output::write_timetree_tree_outputs;
use app_output::{TimetreeEdgeOut, TimetreeNodeOut, TimetreeOutputMaps, TimetreeResult};
use eyre::{Report, WrapErr};
use log::{debug, info, warn};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use treetime::clock::clock_output::write_clock_model;
use treetime::clock::clock_state::ClockState;
use treetime::gtr::get_gtr::{GtrOutput, write_gtr_json};
use treetime::make_error;
use treetime::partition::timetree::partition::PartitionTimetree;
use treetime::seq::div::compute_edge_mutation_counts;
use treetime::seq::mutation::MutationTrack;
use treetime::seq::sink::{SeqItem, SeqSink, SeqTrack};
use treetime::timetree::coalescent::CoalescentOutput;
use treetime::timetree::convergence::optimizer::TraceSink;
use treetime::timetree::pipeline::{self, TimetreeInput, TimetreeParams};
use treetime::timetree::timetree_state::TimetreeState;
use treetime_graph::assign_node_names::assign_node_names;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::FastaWriter;
use treetime_io::nwk::CommentProviders;
use treetime_primitives::AlignmentRecord;
use treetime_utils::io::file::create_file_or_stdout;

pub fn run_timetree_estimation(
  args: &TreetimeTimetreeArgs,
  cancel: &dyn treetime::cancel::Cancel,
  progress: &dyn treetime::progress::ProgressSink,
) -> Result<TimetreeResult, Report> {
  cancel.check()?;
  progress.report("Loading input", 0.0, "");

  let input_data = load_input_data(args)?;
  let input_leaf_order = input_data.input_leaf_order.clone();
  let confidences = input_data.confidences;
  let parse_names = input_data.names;

  let resolved = args.resolve_outputs()?;
  let trace_sink: Option<Box<dyn TraceSink>> = match resolved.non_tree_outputs.get(&OutputSelection::Tracelog) {
    Some(path) => Some(Box::new(TraceCsvSink::new(create_file_or_stdout(path)?)?)),
    None => None,
  };

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
    include_leaves: args.include_leaves || args.reconstruct_tip_states,
    impute_missing_data: args.impute_missing_data || args.reconstruct_tip_states,
    report_ambiguous: args.report_ambiguous,
    zero_based: args.zero_based,
    seed: args.seed,
  };

  let aln_descs = input_data
    .aln
    .iter()
    .flatten()
    .fold(BTreeMap::new(), |mut descs, record| {
      descs
        .entry(record.seq_name.clone())
        .or_insert_with(|| record.desc.clone());
      descs
    });

  let input = TimetreeInput {
    graph: input_data.graph,
    alphabet: input_data.alphabet,
    sequences: input_data
      .aln
      .map(|records| records.into_iter().map(AlignmentRecord::from).collect()),
    dates: input_data.dates,
    branch_lengths: input_data.branch_lengths,
  };

  let reconstructed_nuc_fasta = resolved
    .non_tree_outputs
    .get(&OutputSelection::ReconstructedNucFasta)
    .cloned();
  let recon_sink: Option<Box<dyn SeqSink>> = match &reconstructed_nuc_fasta {
    Some(path) => Some(Box::new(ReconstructedNucSink::new(
      FastaWriter::new(create_file_or_stdout(path)?),
      parse_names.clone(),
      aln_descs.clone(),
    ))),
    None => None,
  };

  let mut output = pipeline::run(&params, input, &parse_names, trace_sink, recon_sink, cancel, progress)
    .map_err(|err| err.into_report())?;
  if let Some(path) = &reconstructed_nuc_fasta {
    info!("Wrote reconstructed nucleotide FASTA to {path}", path = path.display());
  }

  let names = std::mem::take(&mut output.names);
  let branch_lengths_opt = std::mem::take(&mut output.branch_lengths);

  let descs: BTreeMap<GraphNodeKey, Option<String>> = names
    .iter()
    .map(|(&key, name)| {
      let desc = name.as_deref().and_then(|name| aln_descs.get(name).cloned()).flatten();
      (key, desc)
    })
    .collect();

  let mutation_counts = match args.divergence_units {
    DivergenceUnits::Mutations => {
      if output.partitions.is_empty() {
        return make_error!(
          "--divergence-units=mutations requires ancestral reconstruction; \
           incompatible with --branch-length-mode=input"
        );
      }
      let partition = &output.partitions[0];
      let edge_subs = output
        .graph
        .get_edges()
        .map(|edge| {
          let edge_key = edge.key();
          Ok((edge_key, partition.edge_subs(&output.graph, edge_key)?))
        })
        .collect::<Result<BTreeMap<_, _>, Report>>()?;
      Some(compute_edge_mutation_counts(&output.graph, &edge_subs))
    },
    DivergenceUnits::MutationsPerSite => None,
  };

  let pipeline::TimetreeOutput {
    mut graph,
    clock_model,
    confidence_intervals,
    partitions,
    dates,
    gtr,
    model_name,
    coalescent,
    rate_susceptibility_dates,
    clock_branch_lengths,
    clock_state,
    timetree_state,
    ..
  } = output;
  let maps = gather_timetree_output_maps(&graph, &partitions)?;

  progress.report("Writing output", 0.95, "");
  info!("### TreeTime: writing outputs");

  let topology_order = args
    .topology_order
    .resolve_topology_order(&graph, &names, Some(input_leaf_order))?;
  topology_order.apply(&mut graph, &names, &branch_lengths_opt)?;

  let (nodes, edges) = gather_timetree_outputs(
    &graph,
    &clock_state,
    &timetree_state,
    &rate_susceptibility_dates,
    &clock_branch_lengths,
    &names,
    &descs,
    &branch_lengths_opt,
    &confidences,
  );

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::ConfidenceTsv) {
    match confidence_intervals.as_ref() {
      Some(intervals) => {
        write_confidence_intervals_file(intervals, path).wrap_err("Failed to write confidence intervals")?;
        info!("Wrote confidence intervals to {path}", path = path.display());
      },
      None if args.output_confidence_tsv.is_some() => {
        return make_error!(
          "Confidence output requested but no confidence intervals were computed. \
           Use --time-marginal to enable confidence interval computation."
        );
      },
      None => warn!("Skipping confidence-interval output: no confidence intervals were computed (use --time-marginal)"),
    }
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::CoalescentTsv) {
    write_coalescent_output(
      coalescent.as_ref(),
      path,
      args.output_coalescent_tsv.is_some(),
      |output, path| write_coalescent_delimited(output, path, b'\t'),
    )?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::CoalescentCsv) {
    write_coalescent_output(
      coalescent.as_ref(),
      path,
      args.output_coalescent_csv.is_some(),
      |output, path| write_coalescent_delimited(output, path, b','),
    )?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::CoalescentJson) {
    write_coalescent_output(
      coalescent.as_ref(),
      path,
      args.output_coalescent_json.is_some(),
      |output, path| write_coalescent_json(output, path),
    )?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::ClockModel) {
    write_clock_model(&clock_model, path)?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) {
    match (gtr.as_ref(), model_name) {
      (Some(gtr), Some(model_name)) => {
        let gtr_output = GtrOutput::new(gtr, model_name);
        write_gtr_json(&gtr_output, path)?;
      },
      _ if args.output_gtr.is_some() => {
        return make_error!("GTR output requested but no GTR model was fitted. Provide sequence alignment input.");
      },
      _ => warn!("Skipping GTR output: no GTR model was fitted (provide sequence alignment input)"),
    }
  }

  if !resolved.tree_outputs.is_empty() {
    let date_times: BTreeMap<GraphNodeKey, f64> = nodes
      .iter()
      .filter_map(|(key, out)| out.time.map(|time| (*key, time)))
      .collect();
    let date_provider = DateCommentProvider::new(&date_times);
    if maps.root_sequence.is_some() {
      let provider = EdgeMutationCommentProvider::new(&maps.edge_mutations, &graph);
      let providers = CommentProviders::new().with(&provider).with(&date_provider);
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
    } else {
      let providers = CommentProviders::new().with(&date_provider);
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
      args.tree.as_deref(),
      mutation_counts.as_ref(),
      path,
    )?;
    info!("Wrote augur node data JSON to {path}", path = path.display());
  }

  if args.plot_rtt.is_some() {
    return make_error!("--plot-rtt is not yet implemented");
  }

  if args.plot_tree.is_some() {
    return make_error!("--plot-tree is not yet implemented");
  }

  progress.report("Done", 1.0, "");
  Ok(TimetreeResult { graph, nodes, edges })
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

fn gather_timetree_outputs(
  graph: &Graph,
  clock_state: &ClockState,
  timetree_state: &TimetreeState,
  rate_susceptibility_dates: &BTreeMap<GraphNodeKey, [f64; 3]>,
  clock_branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  descs: &BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  confidences: &BTreeMap<GraphNodeKey, Option<f64>>,
) -> (
  BTreeMap<GraphNodeKey, TimetreeNodeOut>,
  BTreeMap<GraphEdgeKey, TimetreeEdgeOut>,
) {
  let nodes = graph
    .get_nodes()
    .map(|node| {
      let key = node.key();
      let clock = clock_state.node(key);
      let out = TimetreeNodeOut {
        name: names[&key].clone(),
        desc: descs[&key].clone(),
        confidence: confidences.get(&key).copied().flatten(),
        time: timetree_state.node(key).time,
        div: clock.div,
        is_outlier: clock.is_outlier,
        bad_branch: timetree_state.node(key).bad_branch,
        rate_susceptibility_dates: rate_susceptibility_dates.get(&key).copied(),
      };
      (key, out)
    })
    .collect();

  let edges = graph
    .get_edges()
    .map(|edge| {
      let key = edge.key();
      let edge_state = timetree_state.edge(key);
      let out = TimetreeEdgeOut {
        branch_length: branch_lengths[&key],
        time_length: edge_state.time_length,
        clock_branch_length: clock_branch_lengths.get(&key).copied(),
        gamma: edge_state.gamma,
      };
      (key, out)
    })
    .collect();

  (nodes, edges)
}

pub(crate) fn gather_timetree_output_maps(
  graph: &Graph,
  partitions: &[PartitionTimetree],
) -> Result<TimetreeOutputMaps, Report> {
  let Some(partition) = partitions.first() else {
    return Ok(TimetreeOutputMaps::default());
  };
  let root_sequence = Some(partition.root_sequence(graph)?);
  let edge_mutations = graph
    .get_edges()
    .map(|edge| {
      let key = edge.key();
      Ok((key, partition.edge_mutations(graph, key, &MutationTrack::Nucleotide)?))
    })
    .collect::<Result<BTreeMap<_, _>, Report>>()?;
  Ok(TimetreeOutputMaps {
    root_sequence,
    edge_mutations,
  })
}

fn write_coalescent_output(
  coalescent: Option<&CoalescentOutput>,
  path: &Path,
  explicit: bool,
  write: impl FnOnce(&CoalescentOutput, &Path) -> Result<(), Report>,
) -> Result<(), Report> {
  match coalescent {
    Some(output) => {
      write(output, path)?;
      info!("Wrote coalescent output to {path}", path = path.display());
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
