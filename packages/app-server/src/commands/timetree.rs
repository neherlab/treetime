use crate::commands::support::{default_topology_order, output_plan, reroot_spec};
use app_output::DateCommentProvider;
use app_output::EdgeMutationCommentProvider;
use app_output::augur_node_data::write_augur_node_data_json;
use app_output::coalescent::write_coalescent_delimited;
use app_output::output_plan::{CommandKind, OutputSelection};
use app_output::timetree_trace::TraceCsvSink;
use app_output::timetree_tree_output::write_timetree_tree_outputs;
use app_output::{TimetreeEdgeOut, TimetreeNodeOut, TimetreeOutputMaps, TimetreeResult};
use eyre::{Report, WrapErr};
use log::{debug, info, warn};
use serde::Deserialize;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use treetime::alphabet::alphabet::{Alphabet, AlphabetName};
use treetime::ancestral::params::MethodAncestral;
use treetime::cancel::Cancel;
use treetime::clock::clock_output::write_clock_model;
use treetime::clock::clock_state::ClockState;
use treetime::clock::date_constraints::load_date_constraints;
use treetime::clock::find_best_root::params::RerootMethod;
use treetime::gtr::get_gtr::{GtrModelName, GtrOutput, write_gtr_json};
use treetime::make_error;
use treetime::optimize::params::BranchLengthMode;
use treetime::partition::timetree::partition::PartitionTimetree;
use treetime::progress::ProgressSink;
use treetime::seq::gap_fill::{GapFill, apply_gap_fill};
use treetime::seq::mutation::MutationTrack;
use treetime::seq::sink::{SeqItem, SeqSink, SeqTrack};
use treetime::timetree::coalescent::CoalescentOutput;
use treetime::timetree::convergence::optimizer::TraceSink;
use treetime::timetree::params::TimeMarginalMode;
use treetime::timetree::pipeline::{self, TimetreeInput, TimetreeParams};
use treetime::timetree::timetree_state::TimetreeState;
use treetime_graph::assign_node_names::assign_node_names;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::csv::{default_metadata_delimiters, default_name_candidates};
use treetime_io::dates_csv::{DatesMap, read_dates};
use treetime_io::fasta::{FastaRecord, FastaWriter, read_many_fasta_path};
use treetime_io::nwk::{CommentProviders, nwk_read_file};
use treetime_primitives::AlignmentRecord;
use treetime_utils::io::file::create_file_or_stdout;
use utoipa::ToSchema;

pub(crate) fn run_timetree(
  args: &TimetreeArgs,
  cancel: &dyn Cancel,
  progress: &dyn ProgressSink,
) -> Result<TimetreeResult, Report> {
  cancel.check()?;
  progress.report("Loading input", 0.0, "");

  let input_data = load_input_data(args)?;
  let input_leaf_order = input_data.input_leaf_order.clone();
  let confidences = input_data.confidences;
  let parse_names = input_data.names;

  let tracelog_override: BTreeMap<OutputSelection, PathBuf> = args
    .tracelog
    .as_ref()
    .map(|path| BTreeMap::from([(OutputSelection::Tracelog, PathBuf::from(path))]))
    .unwrap_or_default();
  let resolved = output_plan(CommandKind::Timetree, Path::new(&args.outdir), tracelog_override)?;
  let trace_sink: Option<Box<dyn TraceSink>> = match resolved.non_tree_outputs.get(&OutputSelection::Tracelog) {
    Some(path) => Some(Box::new(TraceCsvSink::new(create_file_or_stdout(path)?)?)),
    None => None,
  };

  let params = args.params();

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

  let mutation_counts: Option<BTreeMap<GraphEdgeKey, usize>> = None;

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

  default_topology_order().apply(&mut graph, &names, &branch_lengths_opt)?;

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

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::CoalescentTsv) {
    write_coalescent_output(coalescent.as_ref(), path, |output, path| {
      write_coalescent_delimited(output, path, b'\t')
    })?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::ClockModel) {
    write_clock_model(&clock_model, path)?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) {
    match (gtr.as_ref(), model_name) {
      (Some(gtr), Some(model_name)) => {
        let gtr_output = GtrOutput::builder().gtr(gtr).model_name(model_name).build();
        write_gtr_json(&gtr_output, path)?;
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
    let alignment = args.input_fastas.first().map(Path::new);
    write_augur_node_data_json(
      &graph,
      &nodes,
      &edges,
      &clock_model,
      confidence_intervals.as_deref(),
      dates.as_ref(),
      alignment,
      args.tree.as_deref().map(Path::new),
      mutation_counts.as_ref(),
      path,
    )?;
    info!("Wrote augur node data JSON to {path}", path = path.display());
  }

  progress.report("Done", 1.0, "");
  Ok(TimetreeResult { graph, nodes, edges })
}

fn load_input_data(args: &TimetreeArgs) -> Result<InputData, Report> {
  let nwk_parsed = if let Some(tree_path) = &args.tree {
    nwk_read_file(Path::new(tree_path)).wrap_err("Failed to load tree from file")?
  } else {
    return make_error!("Tree inference from alignment not yet implemented");
  };
  let confidences = nwk_parsed.confidences();
  let names = nwk_parsed.names();
  let graph = nwk_parsed.graph;
  let branch_lengths = nwk_parsed.branch_lengths;
  let input_leaf_order = graph
    .get_leaves()
    .map(|leaf| {
      let key = leaf.key();
      names[&key]
        .clone()
        .ok_or_else(|| treetime::make_report!("Leaf node {key} has no name"))
    })
    .collect::<Result<Vec<_>, _>>()?;

  let alphabet = Alphabet::new(args.alphabet)?;

  let aln = if !args.input_fastas.is_empty() {
    let paths: Vec<PathBuf> = args.input_fastas.iter().map(PathBuf::from).collect();
    let mut records = read_many_fasta_path(&paths, &alphabet)?;
    let gap_fill_mode = args.effective_gap_fill();
    for record in &mut records {
      apply_gap_fill(&mut record.seq, gap_fill_mode, alphabet.gap(), alphabet.unknown());
    }
    Some(records)
  } else if args.branch_length_mode != BranchLengthMode::Input {
    return make_error!(
      "Alignment required when branch_length_mode is not 'input'. \
       Provide FASTA files or use branch_length_mode=input"
    );
  } else {
    None
  };

  let dates = if let Some(dates_path) = &args.dates {
    let id_columns = args
      .name_column
      .clone()
      .map_or_else(default_name_candidates, |col| vec![col]);
    let dates = read_dates(
      Path::new(dates_path),
      &default_metadata_delimiters(),
      &id_columns,
      &None,
      &args.date_column,
    )
    .wrap_err("When reading dates")?;
    load_date_constraints(&dates, &graph, &names).wrap_err("Failed to load date constraints")?;
    Some(dates)
  } else {
    None
  };

  Ok(InputData {
    graph,
    confidences,
    names,
    branch_lengths,
    input_leaf_order,
    alphabet,
    aln,
    dates,
  })
}

/// Time-tree estimation request (openapi subset).
#[derive(Debug, SmartDefault, Deserialize, ToSchema)]
#[serde(default)]
pub(crate) struct TimetreeArgs {
  input_fastas: Vec<String>,
  tree: Option<String>,
  vcf_reference: Option<String>,
  dates: Option<String>,
  name_column: Option<String>,
  date_column: Option<String>,
  sequence_length: Option<usize>,
  clock_rate: Option<f64>,
  clock_std_dev: Option<f64>,
  #[default(BranchLengthMode::default())]
  #[schema(value_type = String)]
  branch_length_mode: BranchLengthMode,
  #[default(TimeMarginalMode::default())]
  #[schema(value_type = String)]
  time_marginal: TimeMarginalMode,
  confidence: bool,
  keep_polytomies: bool,
  resolve_polytomies: bool,
  relax: Vec<f64>,
  #[default = 2]
  max_iter: usize,
  coalescent: Option<f64>,
  coalescent_opt: bool,
  coalescent_skyline: bool,
  #[default = 20]
  skyline_n_points: usize,
  #[default = 2.0]
  skyline_stiffness: f64,
  #[default = 2.0]
  coalescent_confidence: f64,
  #[default = 50.0]
  gen_per_year: f64,
  n_branches_posterior: Option<usize>,
  tip_labels: bool,
  no_tip_labels: bool,
  clock_filter: f64,
  n_iqd: Option<f64>,
  #[schema(value_type = Option<String>)]
  reroot: Option<RerootMethod>,
  reroot_tips: Vec<String>,
  keep_root: bool,
  allow_negative_rate: bool,
  tip_slack: Option<f64>,
  covariation: bool,
  #[default(GtrModelName::default())]
  #[schema(value_type = String)]
  gtr: GtrModelName,
  gtr_params: Vec<String>,
  #[default(MethodAncestral::default())]
  #[schema(value_type = String)]
  method_anc: MethodAncestral,
  #[default(AlphabetName::default())]
  #[schema(value_type = String)]
  alphabet: AlphabetName,
  dense: Option<bool>,
  aa: bool,
  #[default(GapFill::default())]
  #[schema(value_type = String)]
  gap_fill: GapFill,
  keep_overhangs: bool,
  zero_based: bool,
  include_leaves: bool,
  impute_missing_data: bool,
  reconstruct_tip_states: bool,
  report_ambiguous: bool,
  no_indels: bool,
  outdir: String,
  tracelog: Option<String>,
  seed: Option<u64>,
}

impl TimetreeArgs {
  fn effective_gap_fill(&self) -> GapFill {
    if self.keep_overhangs {
      GapFill::None
    } else {
      self.gap_fill
    }
  }

  fn params(&self) -> TimetreeParams {
    TimetreeParams {
      model: self.gtr,
      alphabet_name: self.alphabet,
      dense: self.dense,
      gap_fill: self.effective_gap_fill(),
      branch_length_mode: self.branch_length_mode,
      no_indels: self.no_indels,
      sequence_length: self.sequence_length,
      clock_rate: self.clock_rate,
      clock_std_dev: self.clock_std_dev,
      keep_root: self.keep_root,
      reroot_spec: reroot_spec(self.reroot, &self.reroot_tips),
      allow_negative_rate: self.allow_negative_rate,
      clock_filter: self.clock_filter,
      covariation: self.covariation,
      tip_slack: self.tip_slack,
      max_iter: self.max_iter,
      resolve_polytomies: self.resolve_polytomies,
      keep_polytomies: self.keep_polytomies,
      relax: self.relax.clone(),
      coalescent: self.coalescent,
      coalescent_opt: self.coalescent_opt,
      coalescent_skyline: self.coalescent_skyline,
      skyline_n_points: self.skyline_n_points,
      skyline_stiffness: self.skyline_stiffness,
      coalescent_confidence: self.coalescent_confidence,
      gen_per_year: self.gen_per_year,
      n_branches_posterior: self.n_branches_posterior,
      time_marginal: self.time_marginal,
      confidence: self.confidence,
      include_leaves: self.include_leaves || self.reconstruct_tip_states,
      impute_missing_data: self.impute_missing_data || self.reconstruct_tip_states,
      report_ambiguous: self.report_ambiguous,
      zero_based: self.zero_based,
      seed: self.seed,
    }
  }
}

struct InputData {
  graph: Graph,
  confidences: BTreeMap<GraphNodeKey, Option<f64>>,
  names: BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  input_leaf_order: Vec<String>,
  alphabet: Alphabet,
  aln: Option<Vec<FastaRecord>>,
  dates: Option<DatesMap>,
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

#[expect(
  clippy::too_many_arguments,
  reason = "each argument is an independent input of this step; a parameter struct would be built only for this call"
)]
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

fn gather_timetree_output_maps(graph: &Graph, partitions: &[PartitionTimetree]) -> Result<TimetreeOutputMaps, Report> {
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
  write: impl FnOnce(&CoalescentOutput, &Path) -> Result<(), Report>,
) -> Result<(), Report> {
  match coalescent {
    Some(output) => {
      write(output, path)?;
      info!("Wrote coalescent output to {path}", path = path.display());
    },
    None => debug!(
      "Skipping coalescent output: no coalescent model was set \
       (use coalescent, coalescent_opt, or coalescent_skyline)"
    ),
  }
  Ok(())
}
