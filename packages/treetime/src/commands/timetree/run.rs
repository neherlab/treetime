use crate::clock::clock_output::write_clock_model;
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use crate::commands::timetree::initialization::load_input_data;
use crate::commands::timetree::output::augur_node_data::write_augur_node_data_json;
use crate::commands::timetree::output::auspice::write_auspice_json;
use crate::commands::timetree::result::TimetreeResult;
use crate::gtr::get_gtr::write_gtr_json;
use crate::make_error;
use crate::partition::traits::MutationCommentProvider;
use crate::timetree::confidence::write_confidence_intervals_file;
use crate::timetree::pipeline::{self, TimetreeInput, TimetreeParams};
use eyre::{Report, WrapErr};
use log::info;
use std::path::PathBuf;
use treetime_io::graph::write_graph_files_with;
use treetime_io::nwk::CommentProviders;
use treetime_utils::io::file::create_file_or_stdout;

pub fn run_timetree_estimation(
  args: &TreetimeTimetreeArgs,
  progress: &dyn crate::progress::ProgressSink,
) -> Result<TimetreeResult, Report> {
  progress.check_cancelled()?;
  progress.report("Loading input", 0.0, "");

  let input_data = load_input_data(args)?;

  let tracelog: Option<Box<dyn std::io::Write + Send>> = if let Some(path) = &args.tracelog {
    Some(Box::new(create_file_or_stdout(path)?))
  } else {
    None
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

  if let Some(intervals) = &output.confidence_intervals {
    let ci_path = args.output.outdir.join("confidence_intervals.tsv");
    write_confidence_intervals_file(intervals, &ci_path).wrap_err("Failed to write confidence intervals")?;
    info!("Wrote confidence intervals to {ci_path}", ci_path = ci_path.display());
  }

  write_outputs(args, &output)?;

  progress.report("Done", 1.0, "");
  Ok(TimetreeResult {
    graph: output.graph,
    clock_model: output.clock_model,
    confidence_intervals: output.confidence_intervals,
  })
}

fn write_outputs(args: &TreetimeTimetreeArgs, output: &pipeline::TimetreeOutput) -> Result<(), Report> {
  if !output.partitions.is_empty() {
    let guard = output.partitions[0].read_arc();
    let provider = MutationCommentProvider::new(&*guard, &output.graph);
    let providers = CommentProviders::new().with(&provider);
    write_graph_files_with(&args.output.outdir, "timetree", &output.graph, &providers)
      .wrap_err("Failed to write tree output")?;
  } else {
    write_graph_files_with(&args.output.outdir, "timetree", &output.graph, &CommentProviders::new())
      .wrap_err("Failed to write tree output")?;
  }

  write_clock_model(&output.clock_model, &args.output.outdir.join("timetree"))?;

  if let (Some(gtr), Some(model_name)) = (&output.gtr, output.model_name) {
    write_gtr_json(gtr, model_name, &args.output.outdir, None)?;
  }

  write_auspice_json(&output.graph, output.confidence_intervals.as_deref(), &args.output.outdir)?;

  let augur_node_data_path = args
    .output_augur_node_data
    .clone()
    .unwrap_or_else(|| args.output.outdir.join("timetree.augur-node-data.json"));
  let alignment = args.alignment.alignment.first().map(PathBuf::as_path);
  write_augur_node_data_json(
    &output.graph,
    &output.clock_model,
    output.confidence_intervals.as_deref(),
    output.dates.as_ref(),
    alignment,
    args.tree.as_deref(),
    &augur_node_data_path,
  )?;
  info!("Wrote augur node data JSON to {path}", path = augur_node_data_path.display());

  if args.plot_rtt.is_some() {
    return make_error!("--plot-rtt is not yet implemented");
  }

  if args.plot_tree.is_some() {
    return make_error!("--plot-tree is not yet implemented");
  }

  Ok(())
}
