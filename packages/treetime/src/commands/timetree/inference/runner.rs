use crate::commands::optimize::optimize_unified::OptimizationContribution;
use crate::commands::timetree::inference::backward_pass::propagate_distributions_backward;
use crate::commands::timetree::inference::branch_length_likelihood::compute_branch_length_distribution;
use crate::commands::timetree::inference::forward_pass::propagate_distributions_forward;
use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::commands::timetree::utils::initialize_node_divergences;
use crate::commands::{
  clock::clock_model::ClockModel, timetree::coalescent::coalescent::compute_coalescent_contributions,
};
use crate::distribution::distribution::Distribution;
use crate::graph::edge::GraphEdgeKey;
use crate::representation::graph_ancestral::GraphAncestral;
use eyre::Report;
use log::{debug, info};
use parking_lot::RwLock;
use std::sync::Arc;

pub const BRANCH_GRID_SIZE: usize = 200;

pub fn run_timetree(
  graph: &mut GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll>>],
  clock_model: &ClockModel,
  coalescent_tc: Option<f64>,
) -> Result<(), Report> {
  info!("# Running timetree inference");

  info!("## Calculating divergence distances");
  initialize_node_divergences(graph);

  info!("## Using clock model");
  info!("**Clock rate:** {:.6e}", clock_model.clock_rate());

  if !partitions.is_empty() {
    info!("## Computing branch distributions from partitions");
    compute_branch_distributions_marginal_mode(graph, partitions, clock_model)?;
  } else {
    info!("## Creating branch distributions from input lengths");
    create_branch_distributions_inpu_mode(graph, clock_model)?;
  }

  let coalescent_contributions = if let Some(tc) = coalescent_tc {
    info!("## Computing coalescent contributions with Tc = {tc:.6e}");
    Some(compute_coalescent_contributions(graph, &Distribution::constant(tc))?)
  } else {
    None
  };

  info!("## Propagating distributions backward");
  propagate_distributions_backward(graph, coalescent_contributions.as_ref())?;

  info!("## Propagating distributions forward");
  propagate_distributions_forward(graph)?;

  info!("# Timetree inference completed");
  Ok(())
}

fn compute_branch_distributions_marginal_mode(
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll>>],
  clock_model: &ClockModel,
) -> Result<(), Report> {
  // In input branch mode, partitions exist but have no edge data
  // Skip branch distribution computation and use input branch lengths directly
  if partitions.iter().any(|p| p.read_arc().get_sequence_length().is_none()) {
    debug!("Skipping branch distribution computation: partitions not initialized");
    return Ok(());
  }

  let one_mutation = calculate_one_mutation(partitions);
  let total_sites: usize = partitions
    .iter()
    .map(|p| p.read_arc().get_sequence_length().unwrap_or(0))
    .sum();

  info!(
    "Computing branch distributions from {} partition(s) with {} total sites",
    partitions.len(),
    total_sites
  );
  debug!("One mutation = {one_mutation:.6e} substitutions/site");

  let clock_rate = clock_model.clock_rate();

  for edge_ref in graph.get_edges() {
    let edge_key = edge_ref.read_arc().key();
    let mut edge = edge_ref.write_arc().payload().write_arc();
    let branch_length = edge.branch_length.unwrap_or(one_mutation);

    debug!("Edge {edge_key:?}: input branch_length = {branch_length:.6e}");

    let contributions = collect_contributions(partitions, edge_key)?;
    let distribution = compute_branch_length_distribution(
      &contributions,
      branch_length,
      one_mutation,
      BRANCH_GRID_SIZE,
      clock_rate,
    )?;

    if let Some(likely_time) = distribution.likely_time() {
      debug!("Edge {edge_key:?}: distribution peak at time = {likely_time:.6e}");
    }

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

fn create_branch_distributions_inpu_mode(graph: &GraphAncestral, clock_model: &ClockModel) -> Result<(), Report> {
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
