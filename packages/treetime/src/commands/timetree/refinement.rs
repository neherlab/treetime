use crate::commands::ancestral::marginal::update_marginal;
use crate::commands::clock::clock_model::ClockModel;
use crate::commands::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot};
use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use crate::commands::timetree::inference::runner::run_timetree;
use crate::commands::timetree::optimization::reroot::reroot_tree;
use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::representation::partition::timetree::GraphTimetree;
use crate::representation::payload::timetree::EdgeTimetree;
use crate::representation::payload::timetree::NodeTimetree;
use eyre::{Report, WrapErr};
use log::info;
use parking_lot::RwLock;
use std::sync::Arc;

#[allow(clippy::useless_let_if_seq)]
pub fn run_refinement_iteration(
  args: &TreetimeTimetreeArgs,
  graph: &mut GraphTimetree,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>],
  clock_model: &mut ClockModel,
  clock_params: &ClockParams,
  branch_params: &BranchPointOptimizationParams,
) -> Result<(usize, usize), Report> {
  let mut is_tree_dirty = false;

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
    run_timetree(graph, partitions, clock_model, args.coalescent).wrap_err("Timetree inference failed")?;

    if !partitions.is_empty() {
      update_marginal(graph, partitions)?;
    }
  } else {
    if !partitions.is_empty() {
      info!("Updating ancestral sequences via marginal reconstruction");
      update_marginal(graph, partitions)?;
    }

    info!("Updating node times via timetree inference");
    run_timetree(graph, partitions, clock_model, args.coalescent).wrap_err("Timetree inference failed")?;
  }

  let n_diff = 0;

  *clock_model = if args.keep_root {
    estimate_clock_model_with_reroot(graph, clock_params, args.clock_rate, true, branch_params)
      .wrap_err("Failed to update clock model")?
  } else {
    reroot_tree(graph, partitions, clock_params, args.clock_rate, branch_params)
      .wrap_err("Failed to update clock model with reroot")?
  };

  Ok((n_diff, n_resolved))
}
