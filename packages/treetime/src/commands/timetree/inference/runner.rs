use crate::commands::clock::clock_model::ClockModel;
use crate::commands::clock::clock_regression::{ClockOptions, estimate_clock_model_with_reroot};
use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::commands::optimize::optimize_unified::OptimizationContribution;
use crate::commands::timetree::inference::backward_pass::propagate_distributions_backward;
use crate::commands::timetree::inference::branch_length_likelihood::compute_branch_length_distribution;
use crate::commands::timetree::inference::forward_pass::propagate_distributions_forward;
use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::distribution::distribution::Distribution;
use crate::graph::edge::GraphEdgeKey;
use crate::representation::graph_ancestral::GraphAncestral;
use eyre::Report;
use parking_lot::RwLock;
use std::sync::Arc;

const BRANCH_GRID_SIZE: usize = 200;

pub fn run_timetree(
  graph: &mut GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll>>],
  keep_root: bool,
) -> Result<(), Report> {
  log::info!("# Running timetree inference");

  log::info!("## Estimating clock rate from root-to-tip regression");
  let clock_model = estimate_clock_model_with_reroot(
    graph,
    &ClockOptions::default(),
    keep_root,
    &BranchPointOptimizationParams::default(),
  )?;
  log::info!("**Estimated clock rate:** {:.6e}", clock_model.clock_rate());
  log::debug!("Clock rate: {:.6e}", clock_model.clock_rate());

  if !partitions.is_empty() {
    log::info!("## Computing branch distributions from partitions");
    compute_branch_distributions(graph, partitions, &clock_model)?;
  } else {
    log::info!("## Creating branch distributions from input lengths");
    create_branch_distributions_from_input_lengths(graph, &clock_model)?;
  }

  log::info!("## Propagating distributions backward");
  propagate_distributions_backward(graph)?;

  log::info!("## Propagating distributions forward");
  propagate_distributions_forward(graph)?;

  log::info!("# Timetree inference completed");
  Ok(())
}

fn compute_branch_distributions(
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll>>],
  clock_model: &ClockModel,
) -> Result<(), Report> {
  // In input branch mode, partitions exist but have no edge data
  // Skip branch distribution computation and use input branch lengths directly
  if partitions.iter().any(|p| p.read_arc().get_sequence_length().is_none()) {
    return Ok(());
  }

  let one_mutation = calculate_one_mutation(partitions);
  let clock_rate = clock_model.clock_rate();

  for edge_ref in graph.get_edges() {
    let edge_key = edge_ref.read_arc().key();
    let mut edge = edge_ref.write_arc().payload().write_arc();
    let branch_length = edge.branch_length.unwrap_or(one_mutation);

    let contributions = collect_contributions(partitions, edge_key)?;
    let distribution = compute_branch_length_distribution(
      &contributions,
      branch_length,
      one_mutation,
      BRANCH_GRID_SIZE,
      clock_rate,
    )?;
    edge.branch_length_distribution = Some(distribution);
  }
  Ok(())
}

fn calculate_one_mutation(partitions: &[Arc<RwLock<dyn PartitionTimetreeAll>>]) -> f64 {
  let total_length: usize = partitions
    .iter()
    .map(|part| part.read_arc().get_sequence_length().unwrap_or(0))
    .sum();
  1.0 / total_length as f64
}

fn collect_contributions(
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll>>],
  edge_key: GraphEdgeKey,
) -> Result<Vec<OptimizationContribution>, Report> {
  partitions
    .iter()
    .map(|partition| partition.read_arc().create_edge_contribution(edge_key))
    .collect()
}

fn create_branch_distributions_from_input_lengths(
  graph: &GraphAncestral,
  clock_model: &ClockModel,
) -> Result<(), Report> {
  let clock_rate = clock_model.clock_rate();

  for edge_ref in graph.get_edges() {
    let mut edge = edge_ref.write_arc().payload().write_arc();

    if let Some(branch_length) = edge.branch_length {
      // Convert branch length (substitutions/site) to time duration (years)
      let time_duration = branch_length / clock_rate;
      let distribution = Distribution::point(time_duration, 1.0);
      edge.branch_length_distribution = Some(Arc::new(distribution));
    }
  }

  Ok(())
}
