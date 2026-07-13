use crate::alphabet::alphabet::Alphabet;
use crate::commands::optimize::args::TreetimeOptimizeArgs;
use crate::commands::optimize::augur_node_data::write_augur_node_data_json;
use crate::commands::optimize::result::OptimizeResult;
use crate::commands::shared::ir_projection::build_ir_with_mutations;
use crate::commands::shared::output::{CommandKind, DivergenceUnits, OutputSelection};
use crate::gtr::get_gtr::{GtrOutput, write_gtr_json};
use crate::make_error;
use crate::optimize::pipeline::{self, OptimizeInput, OptimizeParams};
use crate::partition::traits::MutationCommentProvider;
use crate::seq::div::compute_edge_mutation_counts;
use crate::seq::gap_fill::apply_gap_fill;
use eyre::Report;
use log::info;
use std::path::PathBuf;
use treetime_io::fasta::read_many_fasta;
use treetime_io::graph::write_tree_outputs;
use treetime_io::nwk::CommentProviders;
use treetime_io::nwk::nwk_read_file;
use treetime_io::tree_ir::types::TreeIrData;

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

  let selection: Vec<OutputSelection> = args
    .output_selection
    .iter()
    .copied()
    .map(OutputSelection::from)
    .collect();
  let resolved = args.output.resolve(
    CommandKind::Optimize,
    &selection,
    &[
      (OutputSelection::AugurNodeData, args.output_augur_node_data.as_deref()),
      (OutputSelection::Gtr, args.output_gtr.as_deref()),
    ],
  )?;

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

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) {
    let gtr_output = GtrOutput::new(&output.gtr, output.model_name);
    write_gtr_json(&gtr_output, path)?;
  }

  if !resolved.tree_outputs.is_empty() {
    let topology_order = args.topology_order.resolve_topology_order(&output.graph, None)?;
    let plan = topology_order.plan(&output.graph)?;
    let ordered = plan.ordered_graph(&output.graph)?;

    if !output.dense_partitions.is_empty() {
      let guard = output.dense_partitions[0].read_arc();
      let provider = MutationCommentProvider::new(&*guard, &output.graph);
      let providers = CommentProviders::new().with(&provider);
      let data = TreeIrData {
        has_mutations: true,
        ..TreeIrData::default()
      };
      let ir = build_ir_with_mutations(&output.graph, &*guard, data)?;
      write_tree_outputs(&ordered, &resolved.tree_outputs, &providers, Some(&ir))?;
    } else if !output.sparse_partitions.is_empty() {
      let guard = output.sparse_partitions[0].read_arc();
      let provider = MutationCommentProvider::new(&*guard, &output.graph);
      let providers = CommentProviders::new().with(&provider);
      let data = TreeIrData {
        has_mutations: true,
        ..TreeIrData::default()
      };
      let ir = build_ir_with_mutations(&output.graph, &*guard, data)?;
      write_tree_outputs(&ordered, &resolved.tree_outputs, &providers, Some(&ir))?;
    } else {
      write_tree_outputs(&ordered, &resolved.tree_outputs, &CommentProviders::new(), None)?;
    }
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::AugurNodeData) {
    let mutation_counts = match args.divergence_units {
      DivergenceUnits::Mutations => {
        let partition: &dyn crate::partition::traits::PartitionBranchOps = if !output.dense_partitions.is_empty() {
          &*output.dense_partitions[0].read_arc()
        } else if !output.sparse_partitions.is_empty() {
          &*output.sparse_partitions[0].read_arc()
        } else {
          return make_error!(
            "--divergence-units=mutations requires ancestral reconstruction but no partitions are available"
          );
        };
        Some(compute_edge_mutation_counts(&output.graph, partition)?)
      },
      DivergenceUnits::MutationsPerSite => None,
    };

    let alignment = args.alignment.alignment.first().map(PathBuf::as_path);
    write_augur_node_data_json(
      &output.graph,
      alignment,
      Some(args.tree.as_path()),
      mutation_counts.as_ref(),
      path,
    )?;
    info!("Wrote augur node data JSON to {path}", path = path.display());
  }

  progress.report("Done", 1.0, "");
  Ok(OptimizeResult { graph: output.graph })
}
