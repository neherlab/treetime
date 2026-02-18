use crate::commands::ancestral::marginal::update_marginal;
use crate::commands::clock::clock_model::ClockModel;
use crate::commands::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot};
use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use crate::commands::timetree::inference::runner::run_timetree;
use crate::commands::timetree::optimization::polytomy::{prepare_tree_after_topology_change, resolve_polytomies};
use crate::commands::timetree::optimization::relaxed_clock::apply_relaxed_clock;
use crate::commands::timetree::optimization::reroot::reroot_tree;
use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::representation::partition::timetree::GraphTimetree;
use crate::representation::payload::timetree::EdgeTimetree;
use crate::representation::payload::timetree::NodeTimetree;
use eyre::{Report, WrapErr};
use log::info;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_distribution::Distribution;

#[allow(clippy::useless_let_if_seq)]
pub fn run_refinement_iteration(
  args: &TreetimeTimetreeArgs,
  graph: &mut GraphTimetree,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>],
  clock_model: &mut ClockModel,
  clock_params: &ClockParams,
  branch_params: &BranchPointOptimizationParams,
  coalescent_tc: Option<&Distribution>,
) -> Result<(usize, usize), Report> {
  let mut is_tree_dirty = false;

  if !args.relax.is_empty() {
    let total_length: usize = partitions.iter().map(|p| p.read_arc().get_sequence_length()).sum();
    let one_mutation = 1.0 / total_length as f64;
    info!(
      "Applying relaxed clock with slack={}, coupling={}",
      args.relax.first().copied().unwrap_or(1.0),
      args.relax.get(1).copied().unwrap_or(1.0)
    );
    apply_relaxed_clock(graph, &args.relax, one_mutation);
  }

  let n_resolved = if args.resolve_polytomies {
    let n = resolve_polytomies(graph, partitions).wrap_err("Polytomy resolution failed")?;
    if n > 0 {
      info!("Resolved polytomies, introduced {n} new nodes");
      prepare_tree_after_topology_change(graph).wrap_err("Failed to prepare tree after topology change")?;
    }
    n
  } else {
    0
  };

  if n_resolved > 0 {
    is_tree_dirty = true;
  }

  if is_tree_dirty {
    info!("Tree structure changed - recomputing timetree then marginal");
    run_timetree(graph, partitions, clock_model, coalescent_tc).wrap_err("Timetree inference failed")?;

    if !partitions.is_empty() {
      update_marginal(graph, partitions)?;
    }
  } else {
    if !partitions.is_empty() {
      info!("Updating ancestral sequences via marginal reconstruction");
      update_marginal(graph, partitions)?;
    }

    info!("Updating node times via timetree inference");
    run_timetree(graph, partitions, clock_model, coalescent_tc).wrap_err("Timetree inference failed")?;
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
