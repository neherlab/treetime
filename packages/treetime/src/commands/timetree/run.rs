use crate::clock::clock_output::write_clock_model;
use crate::commands::shared::output::{CommandKind, DivergenceUnits, OutputSelection};
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use crate::commands::timetree::initialization::load_input_data;
use crate::commands::timetree::output::augur_node_data::write_augur_node_data_json;
use crate::commands::timetree::output::ir::build_timetree_ir;
use crate::commands::timetree::result::TimetreeResult;
use crate::gtr::get_gtr::{GtrOutput, write_gtr_json};
use crate::make_error;
use crate::partition::traits::MutationCommentProvider;
use crate::seq::div::compute_edge_mutation_counts;
use crate::timetree::confidence::write_confidence_intervals_file;
use crate::timetree::pipeline::{self, TimetreeInput, TimetreeParams};
use eyre::{Report, WrapErr};
use log::{info, warn};
use std::path::PathBuf;
use treetime_io::graph::write_tree_outputs;
use treetime_io::nwk::CommentProviders;
use treetime_utils::io::file::create_file_or_stdout;

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
  let selection: Vec<OutputSelection> = args
    .output_selection
    .iter()
    .copied()
    .map(OutputSelection::from)
    .collect();
  let resolved = args.output.resolve(
    CommandKind::Timetree,
    &selection,
    &[
      (OutputSelection::AugurNodeData, args.output_augur_node_data.as_deref()),
      (OutputSelection::Gtr, args.output_gtr.as_deref()),
      (OutputSelection::ClockModel, args.output_clock_model.as_deref()),
      (OutputSelection::ConfidenceTsv, args.output_confidence_tsv.as_deref()),
      (OutputSelection::Tracelog, args.output_tracelog.as_deref()),
    ],
  )?;

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
    n_skyline: args.n_skyline,
    n_branches_posterior: args.n_branches_posterior,
    time_marginal: args.time_marginal,
    confidence: args.confidence,
    reconstruct_tip_states: args.reconstruct_tip_states,
    report_ambiguous: args.report_ambiguous,
    zero_based: args.zero_based,
  };

  let input = TimetreeInput {
    graph: input_data.graph,
    alphabet: input_data.alphabet,
    sequences: input_data.aln,
    dates: input_data.dates,
  };

  let output = pipeline::run(&params, input, tracelog, progress)?;

  progress.report("Writing output", 0.95, "");
  info!("### TreeTime: writing outputs");

  let topology_order = args
    .topology_order
    .resolve_topology_order(&output.graph, Some(input_leaf_order))?;

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::ConfidenceTsv) {
    match output.confidence_intervals.as_ref() {
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

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::ClockModel) {
    write_clock_model(&output.clock_model, path)?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) {
    match (output.gtr.as_ref(), output.model_name) {
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

  if !resolved.tree_outputs.is_empty() {
    let plan = topology_order.plan(&output.graph)?;
    let ordered = plan.ordered_graph(&output.graph)?;
    // Build the IR from the ordered clone. The clone preserves original node and edge
    // keys, keeping confidence intervals and mutation counts aligned while applying
    // the same topology order to Auspice and the other TreeIR-backed formats.
    let ir = build_timetree_ir(
      &ordered,
      output.confidence_intervals.as_deref(),
      mutation_counts.as_ref(),
    )?;

    if !output.partitions.is_empty() {
      let guard = output.partitions[0].read_arc();
      let provider = MutationCommentProvider::new(&*guard, &output.graph);
      let providers = CommentProviders::new().with(&provider);
      write_tree_outputs(&ordered, &resolved.tree_outputs, &providers, Some(&ir))?;
    } else {
      write_tree_outputs(&ordered, &resolved.tree_outputs, &CommentProviders::new(), Some(&ir))?;
    }
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::AugurNodeData) {
    let alignment = args.alignment.alignment.first().map(PathBuf::as_path);
    write_augur_node_data_json(
      &output.graph,
      &output.clock_model,
      output.confidence_intervals.as_deref(),
      output.dates.as_ref(),
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
  Ok(TimetreeResult {
    graph: output.graph,
    clock_model: output.clock_model,
    confidence_intervals: output.confidence_intervals,
  })
}
