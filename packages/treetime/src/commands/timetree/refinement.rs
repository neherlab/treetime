use crate::commands::ancestral::marginal_unified::run_marginal;
use crate::commands::clock::clock_model::ClockModel;
use crate::commands::clock::clock_regression::{ClockOptions, estimate_clock_model_with_reroot};
use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use crate::commands::timetree::inference::runner::run_timetree;
use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::io::fasta::FastaRecord;
use crate::representation::graph_ancestral::GraphAncestral;
use eyre::{Report, WrapErr};
use log::info;
use parking_lot::RwLock;
use std::sync::Arc;

#[allow(clippy::useless_let_if_seq)]
pub fn run_refinement_iteration(
  args: &TreetimeTimetreeArgs,
  graph: &mut GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll>>],
  aln: Option<&[FastaRecord]>,
  clock_model: &mut ClockModel,
  clock_options: &ClockOptions,
  branch_params: &BranchPointOptimizationParams,
) -> Result<(usize, usize), Report> {
  let mut is_tree_dirty = false;

  if let Some(_coalescent_params) = &args.coalescent {
    todo!("add_coalescent_model not yet implemented");
  }

  if !args.relax.is_empty() {
    todo!("apply_relaxed_clock not yet implemented");
  }

  let n_resolved = if args.resolve_polytomies {
    todo!("resolve_polytomies not yet implemented");
  } else {
    0
  };

  if n_resolved > 0 {
    is_tree_dirty = true;
  }

  if is_tree_dirty {
    info!("Tree structure changed - recomputing timetree then marginal");
    run_timetree(graph, partitions, clock_model).wrap_err("Timetree inference failed")?;

    if aln.is_some() {
      run_marginal(graph, partitions, aln)?;
    }
  } else {
    if aln.is_some() {
      info!("Updating ancestral sequences via marginal reconstruction");
      run_marginal(graph, partitions, aln)?;
    }

    info!("Updating node times via timetree inference");
    run_timetree(graph, partitions, clock_model).wrap_err("Timetree inference failed")?;
  }

  let n_diff = 0;

  *clock_model =
    estimate_clock_model_with_reroot(graph, clock_options, args.clock_rate, args.keep_root, branch_params)
      .wrap_err("Failed to update clock model")?;

  Ok((n_diff, n_resolved))
}
