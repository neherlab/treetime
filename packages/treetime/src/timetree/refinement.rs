use crate::ancestral::marginal::update_marginal;
use crate::clock::clock_model::ClockModel;
use crate::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot};
use crate::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::partition::timetree::GraphTimetree;
use crate::partition::traits::PartitionTimetreeAll;
use crate::payload::timetree::EdgeTimetree;
use crate::payload::timetree::NodeTimetree;
use crate::timetree::convergence::sequence_changes::{capture_ancestral_states, count_sequence_changes};
use crate::timetree::inference::runner::run_timetree;
use crate::timetree::optimization::polytomy::{prepare_tree_after_topology_change, resolve_polytomies};
use crate::timetree::optimization::relaxed_clock::apply_relaxed_clock;
use eyre::{Report, WrapErr};
use log::info;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_graph::assign_node_names::assign_node_names;

pub struct RefinementParams {
  pub relax: Vec<f64>,
  pub resolve_polytomies: bool,
  pub clock_rate: Option<f64>,
  pub no_indels: bool,
}

pub fn run_refinement_iteration(
  params: &RefinementParams,
  graph: &mut GraphTimetree,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>],
  clock_model: &mut ClockModel,
  clock_params: &ClockParams,
  branch_params: &BranchPointOptimizationParams,
  coalescent_tc: Option<&Distribution>,
) -> Result<(usize, usize), Report> {
  let total_length: usize = partitions.iter().map(|p| p.read_arc().get_sequence_length()).sum();

  if !params.relax.is_empty() {
    if total_length == 0 {
      info!("Skipping relaxed clock: no sequence data (partitions empty or zero-length)");
    } else {
      let one_mutation = 1.0 / total_length as f64;
      info!(
        "Applying relaxed clock with slack={}, coupling={}",
        params.relax.first().copied().unwrap_or(1.0),
        params.relax.get(1).copied().unwrap_or(1.0)
      );
      apply_relaxed_clock(graph, &params.relax, one_mutation, clock_model.clock_rate())?;
    }
  }

  let prev_states = capture_ancestral_states(graph, partitions);

  let n_resolved = if params.resolve_polytomies {
    let zero_branch_slope = clock_model.clock_rate() * total_length as f64;
    let n = resolve_polytomies(graph, partitions, zero_branch_slope, clock_model.clock_rate())
      .wrap_err("Polytomy resolution failed")?;
    if n > 0 {
      info!("Resolved polytomies, introduced {n} new nodes");
      assign_node_names(graph)?;
      prepare_tree_after_topology_change(graph).wrap_err("Failed to prepare tree after topology change")?;

      for partition in partitions {
        partition.write_arc().reconcile_topology(graph);
      }
    }
    n
  } else {
    0
  };

  let is_tree_dirty = n_resolved > 0;

  if is_tree_dirty {
    info!("Tree structure changed - recomputing timetree then marginal");
    run_timetree(graph, partitions, clock_model, coalescent_tc, params.no_indels)
      .wrap_err("Timetree inference failed")?;

    if !partitions.is_empty() {
      update_marginal(graph, partitions)?;
    }
  } else {
    if !partitions.is_empty() {
      info!("Updating ancestral sequences via marginal reconstruction");
      update_marginal(graph, partitions)?;
    }

    info!("Updating node times via timetree inference");
    run_timetree(graph, partitions, clock_model, coalescent_tc, params.no_indels)
      .wrap_err("Timetree inference failed")?;
  }

  let curr_states = capture_ancestral_states(graph, partitions);
  let n_diff = count_sequence_changes(&prev_states, &curr_states);

  *clock_model = estimate_clock_model_with_reroot(
    graph,
    clock_params,
    params.clock_rate,
    true,
    branch_params,
    Some(clock_model.clock_rate()),
  )
  .wrap_err("Failed to update clock model")?;

  Ok((n_diff, n_resolved))
}
