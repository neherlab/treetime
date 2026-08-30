use crate::ancestral::marginal::{ancestral_reconstruction_marginal, update_marginal};
use crate::ancestral::sample::SampleMode;
use crate::clock::clock_output::write_clock_model;
use crate::commands::shared::output::{DivergenceUnits, OutputSelection};
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use crate::commands::shared::tree_output::write_timetree_tree_outputs;
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use crate::commands::timetree::initialization::load_input_data;
use crate::commands::timetree::output::augur_node_data::write_augur_node_data_json;
use crate::commands::timetree::output::coalescent::{
  CoalescentOutput, write_coalescent_delimited, write_coalescent_json,
};
use crate::commands::timetree::result::{TimetreeGraphData, TimetreeResult};
use crate::gtr::get_gtr::{GtrOutput, write_gtr_json};
use crate::make_error;
use crate::partition::traits::MutationCommentProvider;
use crate::seq::div::compute_edge_mutation_counts;
use crate::timetree::confidence::write_confidence_intervals_file;
use crate::timetree::pipeline::{self, TimetreeInput, TimetreeParams};
use eyre::{Report, WrapErr};
use log::{debug, info, warn};
use std::path::{Path, PathBuf};
use treetime_graph::assign_node_names::assign_node_names;
use treetime_graph::node::{Described, Named};
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

  let input = TimetreeInput {
    graph: input_data.graph,
    alphabet: input_data.alphabet,
    sequences: input_data.aln,
    dates: input_data.dates,
  };

  let output = pipeline::run(&params, input, tracelog, progress)?;

  // Name any unnamed internal node before serialization. Rerooting introduces a fresh root node
  // that the load-time naming pass never saw, and polytomy resolution only re-names when it changes
  // the topology, so on a run without polytomy resolution the rerooted root reaches output unnamed.
  // v0 and the `ancestral` command give every internal node a `NODE_<n>` name; assigning one here
  // keeps the reconstructed FASTA, the augur node data, and the tree outputs consistent and matches
  // v0 rather than leaking an empty label or a key-derived placeholder.
  assign_node_names(&output.graph)?;

  // Reconstructed ancestral-sequence FASTA: the v1 equivalent of v0's `ancestral_sequences.fasta`.
  // The timetree pipeline computes marginal posteriors for branch-length optimization but never
  // materializes the flag-aware per-node sequences, so this reuses the same reconstruction the
  // `ancestral` command runs: `update_marginal` refreshes the posteriors against the final branch
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
      update_marginal(&output.graph, &output.partitions)?;
      let mut rng = get_random_number_generator(params.seed);
      ancestral_reconstruction_marginal(
        &output.graph,
        params.include_leaves,
        params.impute_missing_data,
        &output.partitions,
        SampleMode::Argmax,
        &mut rng,
        |node, seq| match writer.as_mut() {
          Some(writer) => {
            let name = node.name().map(|n| n.as_ref().to_owned()).unwrap_or_default();
            writer.write(&name, node.desc(), seq)
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
      let guard = output.partitions[0].read_arc();
      Some(compute_edge_mutation_counts(&output.graph, &*guard)?)
    },
    DivergenceUnits::MutationsPerSite => None,
  };

  let pipeline::TimetreeOutput {
    graph,
    clock_model,
    confidence_intervals,
    partitions,
    dates,
    gtr,
    model_name,
    coalescent,
  } = output;
  let mut graph = graph.map_data(TimetreeGraphData::new(
    clock_model,
    confidence_intervals,
    partitions,
    dates,
    gtr,
    model_name,
    mutation_counts,
  ));

  progress.report("Writing output", 0.95, "");
  info!("### TreeTime: writing outputs");

  let topology_order = args
    .topology_order
    .resolve_topology_order(&graph, Some(input_leaf_order))?;
  topology_order.apply(&mut graph)?;

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::ConfidenceTsv) {
    match graph.data().confidence_intervals.as_ref() {
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
    write_clock_model(&graph.data().clock_model, path)?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) {
    match (graph.data().gtr.as_ref(), graph.data().model_name) {
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
    if !graph.data().partitions.is_empty() {
      let guard = graph.data().partitions[0].read_arc();
      let provider = MutationCommentProvider::new(&*guard, &graph);
      let providers = CommentProviders::new().with(&provider);
      write_timetree_tree_outputs(&graph, &resolved.tree_outputs, &providers)?;
    } else {
      write_timetree_tree_outputs(&graph, &resolved.tree_outputs, &CommentProviders::new())?;
    }
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::AugurNodeData) {
    let alignment = args.alignment.alignment.first().map(PathBuf::as_path);
    write_augur_node_data_json(
      &graph,
      &graph.data().clock_model,
      graph.data().confidence_intervals.as_deref(),
      graph.data().dates.as_ref(),
      alignment,
      args.tree.as_deref(),
      graph.data().mutation_counts.as_ref(),
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
  Ok(TimetreeResult { graph })
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
