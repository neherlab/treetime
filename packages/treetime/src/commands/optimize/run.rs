use crate::alphabet::alphabet::Alphabet;
use crate::commands::optimize::args::TreetimeOptimizeArgs;
use crate::commands::optimize::augur_node_data::write_augur_node_data_json;
use crate::commands::optimize::result::OptimizeResult;
use crate::gtr::get_gtr::write_gtr_json;
use crate::optimize::pipeline::{self, OptimizeInput, OptimizeParams};
use crate::partition::traits::MutationCommentProvider;
use crate::seq::gap_fill::apply_gap_fill;
use eyre::Report;
use log::info;
use std::path::PathBuf;
use treetime_io::fasta::read_many_fasta;
use treetime_io::graph::write_graph_files_with;
use treetime_io::nwk::CommentProviders;
use treetime_io::nwk::nwk_read_file;

pub fn run_optimize(
  args: &TreetimeOptimizeArgs,
  progress: &dyn crate::progress::ProgressSink,
) -> Result<OptimizeResult, Report> {
  progress.check_cancelled()?;
  progress.report("Reading input", 0.0, "");

  let alphabet = Alphabet::new(args.alphabet_args.alphabet.unwrap_or_default())?;
  let gap_fill = args.gap_fill_args.effective_gap_fill();
  let mut aln = read_many_fasta(&args.alignment.alignment, &alphabet)?;
  for record in &mut aln {
    apply_gap_fill(&mut record.seq, gap_fill, alphabet.gap(), alphabet.unknown());
  }
  let graph = nwk_read_file(&args.tree)?;

  let params = OptimizeParams {
    model: args.model_args.model,
    dense: args.dense,
    max_iter: args.max_iter,
    dp: args.dp,
    damping: args.damping,
    opt_method: args.opt_method,
    initial_guess: args.branch_length_initial_guess,
    no_indels: args.no_indels,
  };

  let input = OptimizeInput {
    graph,
    alphabet,
    sequences: aln,
  };

  let output = pipeline::run(&params, input, progress)?;

  progress.report("Writing output", 0.9, "");
  let outdir = &args.output.outdir;

  write_gtr_json(&output.gtr, output.model_name, outdir, None)?;

  let dense = !output.sparse_partitions.is_empty();
  if !output.dense_partitions.is_empty() {
    let guard = output.dense_partitions[0].read_arc();
    let provider = MutationCommentProvider::new(&*guard, &output.graph);
    let providers = CommentProviders::new().with(&provider);
    write_graph_files_with(outdir, "annotated_tree", &output.graph, &providers)?;
  } else if !output.sparse_partitions.is_empty() {
    let guard = output.sparse_partitions[0].read_arc();
    let provider = MutationCommentProvider::new(&*guard, &output.graph);
    let providers = CommentProviders::new().with(&provider);
    write_graph_files_with(outdir, "annotated_tree", &output.graph, &providers)?;
  }

  let augur_node_data_path = args
    .output_augur_node_data
    .clone()
    .unwrap_or_else(|| outdir.join("optimize.augur-node-data.json"));
  let alignment = args.alignment.alignment.first().map(PathBuf::as_path);
  write_augur_node_data_json(&output.graph, alignment, Some(args.tree.as_path()), &augur_node_data_path)?;
  info!("Wrote augur node data JSON to {path}", path = augur_node_data_path.display());

  progress.report("Done", 1.0, "");
  Ok(OptimizeResult { graph: output.graph })
}
