use crate::ancestral::sample::SampleMode;
use crate::clock::clock_output::write_clock_model;
use crate::clock::clock_state::ClockState;
use crate::commands::shared::output::{DivergenceUnits, OutputSelection};
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use crate::commands::shared::tree_output::write_timetree_tree_outputs;
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use crate::commands::timetree::initialization::load_input_data;
use crate::commands::timetree::output::augur_node_data::write_augur_node_data_json;
use crate::commands::timetree::output::coalescent::{
  CoalescentOutput, write_coalescent_delimited, write_coalescent_json,
};
use crate::commands::timetree::output::date_comment::DateCommentProvider;
use crate::commands::timetree::result::{TimetreeEdgeOut, TimetreeNodeOut, TimetreeOutputMaps, TimetreeResult};
use crate::gtr::get_gtr::{GtrOutput, write_gtr_json};
use crate::make_error;
use crate::partition::timetree::marginal::{ancestral_reconstruction_timetree, marginal_update_timetree};
use crate::partition::timetree::partition::PartitionTimetree;
use crate::partition::traits::{EdgeMutationCommentProvider, PartitionBranchOps};
use crate::seq::div::compute_edge_mutation_counts;
use crate::seq::mutation::MutationTrack;
use crate::timetree::confidence::write_confidence_intervals_file;
use crate::timetree::inference::runner::timetree_branch_lengths;
use crate::timetree::pipeline::{self, TimetreeInput, TimetreeParams};
use crate::timetree::timetree_state::TimetreeState;
use eyre::{Report, WrapErr};
use log::{debug, info, warn};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use treetime_graph::assign_node_names::assign_node_names;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::FastaWriter;
use treetime_io::nwk::CommentProviders;
use treetime_utils::io::file::create_file_or_stdout;
use treetime_utils::sync::random::get_random_number_generator;

pub fn run_timetree_estimation(
  args: &TreetimeTimetreeArgs,
  progress: &dyn crate::progress::ProgressSink,
) -> Result<TimetreeResult, Report> {
  progress.check_cancelled()?;
  progress.report("Loading input", 0.0, "");

  let input_data = load_input_data(args)?;
  let input_leaf_order = input_data.input_leaf_order.clone();
  let confidences = input_data.confidences;
  let parse_names = input_data.names;

  // Resolve outputs up front so the tracelog path (which the pipeline writes during the run) is
  // known before the pipeline starts. Topology ordering is resolved separately, after the pipeline.
  let resolved = args.resolve_outputs()?;
  let tracelog: Option<Box<dyn std::io::Write + Send>> = match resolved.non_tree_outputs.get(&OutputSelection::Tracelog)
  {
    Some(path) => Some(Box::new(create_file_or_stdout(path)?)),
    None => None,
  };

  let params = TimetreeParams {
    model: args.model_args.model,
    alphabet_name: args.alphabet_args.alphabet.unwrap_or_default(),
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

  // Description source for the reconstructed-FASTA writer and the node-output gather. Descriptions
  // live only on the input leaf FASTA records; partition init keys them to leaves by matching each
  // leaf name to its record (first record wins on a duplicate name). Capture that same name-keyed
  // map here, before the alignment moves into the pipeline, so the consumers can rebuild a
  // node-keyed description map from the final `names` map.
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
    sequences: input_data.aln,
    dates: input_data.dates,
    branch_lengths: input_data.branch_lengths,
  };

  let mut output = pipeline::run(&params, input, &parse_names, tracelog, progress)?;

  // Name any unnamed internal node before serialization, and capture the resulting node-name map.
  // Rerooting introduces a fresh root node that the load-time naming pass never saw, and polytomy
  // resolution only re-names when it changes the topology, so on a run without polytomy resolution
  // the rerooted root reaches output unnamed. v0 and the `ancestral` command give every internal
  // node a `NODE_<n>` name; assigning one here keeps the reconstructed FASTA, the augur node data,
  // and the tree outputs consistent and matches v0 rather than leaking an empty label or a
  // key-derived placeholder.
  //
  // `assign_node_names` returns `node_names` taken after the write, so this is also the post-mutation
  // re-snapshot: every downstream reader -- the reconstructed FASTA writer, the gather, and the
  // tree/augur output writers -- reads each node's name from this map. The
  // later `marginal_update` and reconstruction touch neither names nor branch lengths, and topology
  // ordering only permutes keys, so the map still describes the final tree at every later point.
  let names = assign_node_names(output.names, &output.graph)?;
  let branch_lengths_opt = output.branch_lengths;

  // Node-keyed descriptions for the reconstructed-FASTA writer and the node-output gather, rebuilt
  // from the name-keyed `aln_descs` captured from the input alignment. A leaf resolves to its FASTA
  // record's description; an internal node (including the fresh root named above, which matches no
  // record) resolves to `None`. This reproduces exactly what partition init wrote onto each leaf by
  // the same name-to-record match.
  let descs: BTreeMap<GraphNodeKey, Option<String>> = names
    .iter()
    .map(|(&key, name)| {
      let desc = name.as_deref().and_then(|name| aln_descs.get(name).cloned()).flatten();
      (key, desc)
    })
    .collect();

  // Reconstructed ancestral-sequence FASTA: the v1 equivalent of v0's `ancestral_sequences.fasta`.
  // The timetree pipeline computes marginal posteriors for branch-length optimization but never
  // materializes the flag-aware per-node sequences, so this reuses the same reconstruction the
  // `ancestral` command runs: `marginal_update` refreshes the posteriors against the final branch
  // lengths, then `ancestral_reconstruction_marginal` writes each node's stored sequence and emits
  // it. `--include-leaves` gates whether tip sequences are emitted; `--impute-missing-data` resolves
  // ambiguous tip states. Reconstruction runs before mutation counting and tree output so every
  // sequence-derived output reflects the same flag-aware states. The pass is opt-in (FASTA requested
  // or a tip-state flag set), so runs that ask for neither leave all other outputs unchanged.
  let reconstructed_nuc_fasta = resolved.non_tree_outputs.get(&OutputSelection::ReconstructedNucFasta);
  if reconstructed_nuc_fasta.is_some() || params.include_leaves || params.impute_missing_data {
    if output.partitions.is_empty() {
      if reconstructed_nuc_fasta.is_some() {
        return make_error!(
          "Reconstructed sequence output requires ancestral reconstruction; \
           incompatible with --branch-length-mode=input"
        );
      }
      warn!(
        "Ignoring tip-state flags (--include-leaves / --impute-missing-data / --reconstruct-tip-states): \
         no ancestral reconstruction was performed under --branch-length-mode=input"
      );
    } else {
      let mut writer = match reconstructed_nuc_fasta {
        Some(path) => Some(FastaWriter::new(create_file_or_stdout(path)?)),
        None => None,
      };
      marginal_update_timetree(
        &output.graph,
        &timetree_branch_lengths(&output.graph, &branch_lengths_opt, &output.clock_branch_lengths),
        &mut output.partitions,
      )?;
      let mut rng = get_random_number_generator(params.seed);
      ancestral_reconstruction_timetree(
        &output.graph,
        params.include_leaves,
        params.impute_missing_data,
        &mut output.partitions,
        SampleMode::Argmax,
        &mut rng,
        |key, seq| match writer.as_mut() {
          Some(writer) => {
            let name = names[&key].clone().unwrap_or_default();
            writer.write(&name, &descs[&key], seq)
          },
          None => Ok(()),
        },
      )?;
      if let Some(path) = reconstructed_nuc_fasta {
        info!("Wrote reconstructed nucleotide FASTA to {path}", path = path.display());
      }
    }
  }

  let mutation_counts = match args.divergence_units {
    DivergenceUnits::Mutations => {
      if output.partitions.is_empty() {
        return make_error!(
          "--divergence-units=mutations requires ancestral reconstruction; \
           incompatible with --branch-length-mode=input"
        );
      }
      Some(compute_edge_mutation_counts(&output.graph, &output.partitions[0])?)
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
  // Gather the per-node/per-edge sequence and mutation values off the pipeline-local partitions into
  // plain value maps the tree writers consume, taking the partition read out of the serialization
  // path. Node and edge keys stay stable through topology ordering, so gathering before it is
  // bit-identical.
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
    // The `date` Newick/Nexus comment, supplied from the gathered node times as a value.
    let date_times: BTreeMap<GraphNodeKey, f64> = nodes
      .iter()
      .filter_map(|(key, out)| out.time.map(|time| (*key, time)))
      .collect();
    let date_provider = DateCommentProvider::new(&date_times);
    // The mutation comment provider now reads the gathered per-edge mutation map rather than the
    // partition; the partition-less case emits only the date comment.
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
        dates.as_ref(),
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
        dates.as_ref(),
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

/// Gather the per-node and per-edge timetree outputs off the final tree into keyed value maps.
///
/// Runs after the pipeline and topology ordering, so it reads the final divergence, estimated time,
/// exclusion flags, and per-edge branch lengths and relaxed-clock rates of the ordered node set. The
/// value states carry the durable results as values; this step surfaces them as a standalone value
/// the output writers consume.
///
/// `rate_susceptibility_dates` carries the per-node date triples the pipeline returns as a value;
/// each node's triple is read from here. `clock_branch_lengths` likewise
/// carries the committed clock branch length per edge as a value; each edge's clock length is read
/// from here. `clock_state` carries each node's divergence and outlier
/// flag as values; both are read from here. `timetree_state` carries each
/// node's committed time as a value; the time is read from here. `names`
/// and `branch_lengths` are the post-mutation node-name and per-edge branch-length maps captured
/// after the final naming pass; each node's name and each edge's branch length is read from them.
/// `descs` is the node-keyed description map rebuilt from the input
/// alignment; each node's description is read from here. `confidences`
/// carries the parse-time input-tree branch support per node; each node's input branch support is
/// read from here, and a node the pipeline created after the parse is
/// absent and reads as `None`.
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
    .iter()
    .map(|node| {
      let node = node.read_arc();
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
    .iter()
    .map(|edge| {
      let key = edge.read_arc().key();
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

/// Gather the per-node nucleotide sequences, root sequence, and per-edge nucleotide mutations the tree
/// writers read off the timetree partition.
pub(crate) fn gather_timetree_output_maps(
  graph: &Graph,
  partitions: &[PartitionTimetree],
) -> Result<TimetreeOutputMaps, Report> {
  let Some(partition) = partitions.first() else {
    return Ok(TimetreeOutputMaps::default());
  };
  let root_sequence = Some(partition.root_sequence(graph)?);
  let node_sequences = graph
    .get_nodes()
    .iter()
    .map(|node| {
      let key = node.read_arc().key();
      (key, partition.node_sequence(key))
    })
    .collect();
  let edge_mutations = graph
    .get_edges()
    .iter()
    .map(|edge| {
      let key = edge.read_arc().key();
      Ok((key, partition.edge_mutations(graph, key, MutationTrack::Nucleotide)?))
    })
    .collect::<Result<BTreeMap<_, _>, Report>>()?;
  Ok(TimetreeOutputMaps {
    root_sequence,
    node_sequences,
    edge_mutations,
  })
}

/// Writes one coalescent output file, or reports the absence of a coalescent.
///
/// A coalescent output exists only when the run inferred one (`--coalescent`, `--coalescent-opt`,
/// or `--coalescent-skyline`). An explicit per-file flag on a run without a coalescent is an error;
/// a file selected only through `--output-all` or `--output-selection` is skipped with a warning.
/// Mirrors the GTR and confidence-interval outputs.
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
    // A coalescent is opt-in, so its absence under a plain `--output-all` run is expected. Unlike
    // GTR and confidence outputs (near-universal in timetree, so worth a warning when missing),
    // report the skip at debug level to avoid warning noise on every non-coalescent run.
    None => debug!(
      "Skipping coalescent output: no coalescent model was set \
       (use --coalescent, --coalescent-opt, or --coalescent-skyline)"
    ),
  }
  Ok(())
}
