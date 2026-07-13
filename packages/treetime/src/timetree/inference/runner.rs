use crate::clock::clock_model::ClockModel;
use crate::coalescent::coalescent::compute_coalescent_contributions;
use crate::optimize::indel::estimate_indel_rate;
use crate::partition::optimization_contribution::OptimizationContribution;
use crate::partition::traits::PartitionTimetreeAll;
use crate::payload::traits::ClockNode;
use crate::payload::traits::{TimetreeEdge, TimetreeNode};
use crate::timetree::inference::backward_pass::propagate_distributions_backward;
use crate::timetree::inference::branch_length_likelihood::compute_branch_length_distribution;
use crate::timetree::inference::forward_pass::propagate_distributions_forward;
use crate::timetree::utils::initialize_node_divergences;
use eyre::Report;
use log::{debug, info};
use parking_lot::RwLock;
use rayon::prelude::*;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};

pub const BRANCH_GRID_SIZE: usize = 1000;

pub fn run_timetree<N, E, P>(
  graph: &mut Graph<N, E, ()>,
  partitions: &[Arc<RwLock<P>>],
  clock_model: &ClockModel,
  coalescent_tc: Option<&Distribution>,
  no_indels: bool,
) -> Result<(), Report>
where
  N: GraphNode + Named + TimetreeNode + ClockNode,
  E: GraphEdge + HasBranchLength + TimetreeEdge,
  P: PartitionTimetreeAll<N, E> + ?Sized,
{
  info!("# Running timetree inference");

  info!("## Calculating divergence distances");
  initialize_node_divergences(graph)?;

  info!("## Using clock model");
  let clock_rate = clock_model.clock_rate();
  info!("**Clock rate:** {clock_rate:.6e}");

  if !partitions.is_empty() {
    info!("## Computing branch distributions from partitions");
    compute_branch_distributions_marginal_mode(graph, partitions, clock_rate, no_indels)?;
  } else {
    info!("## Creating branch distributions from input lengths");
    create_branch_distributions_input_mode(graph, clock_rate)?;
  }

  let coalescent_contributions = if let Some(tc) = coalescent_tc {
    info!("## Computing coalescent contributions");
    Some(compute_coalescent_contributions(graph, tc)?)
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

fn compute_branch_distributions_marginal_mode<N, E, P>(
  graph: &Graph<N, E, ()>,
  partitions: &[Arc<RwLock<P>>],
  clock_rate: f64,
  no_indels: bool,
) -> Result<(), Report>
where
  N: GraphNode + Named + TimetreeNode,
  E: GraphEdge + HasBranchLength + TimetreeEdge,
  P: PartitionTimetreeAll<N, E> + ?Sized,
{
  let one_mutation = calculate_one_mutation(partitions);
  let total_sites: usize = partitions.iter().map(|p| p.read_arc().get_sequence_length()).sum();

  let indel_rate = if no_indels {
    0.0
  } else {
    estimate_indel_rate(graph, partitions)
  };

  info!(
    "Computing branch distributions from {} partition(s) with {} total sites",
    partitions.len(),
    total_sites
  );
  debug!("One mutation = {one_mutation:.6e} substitutions/site");
  debug!("Indel rate = {indel_rate:.6e} indels/(site*time)");

  graph
    .get_edges()
    .par_iter()
    .try_for_each(|edge_ref| -> Result<(), Report> {
      let edge_key = edge_ref.read_arc().key();
      let mut edge = edge_ref.write_arc().payload().write_arc();
      let branch_length = edge.branch_length().unwrap_or(one_mutation);
      let gamma = edge.gamma();

      debug!("Edge {edge_key:?}: input branch_length = {branch_length:.6e}, gamma = {gamma:.4}");

      let contributions = collect_contributions(partitions, edge_key)?;
      let indel_count: usize = if no_indels {
        0
      } else {
        partitions
          .iter()
          .map(|partition| partition.read_arc().edge_indel_count(edge_key))
          .sum()
      };
      let distribution = compute_branch_length_distribution(
        &contributions,
        indel_count,
        indel_rate,
        branch_length,
        one_mutation,
        BRANCH_GRID_SIZE,
        clock_rate,
        gamma,
      )?;

      if let Some(likely_time) = distribution.likely_time() {
        debug!("Edge {edge_key:?}: distribution peak at time = {likely_time:.6e}");
      }

      edge.set_time_length(distribution.likely_time());
      edge.set_branch_length_distribution(Some(distribution));
      Ok(())
    })?;
  Ok(())
}

fn calculate_one_mutation<N, E, P>(partitions: &[Arc<RwLock<P>>]) -> f64
where
  N: GraphNode + Named,
  E: GraphEdge + HasBranchLength,
  P: PartitionTimetreeAll<N, E> + ?Sized,
{
  let total_length: usize = partitions
    .iter()
    .map(|part| part.read_arc().get_sequence_length())
    .sum();
  1.0 / total_length as f64
}

fn collect_contributions<N, E, P>(
  partitions: &[Arc<RwLock<P>>],
  edge_key: GraphEdgeKey,
) -> Result<Vec<OptimizationContribution>, Report>
where
  N: GraphNode + Named,
  E: GraphEdge + HasBranchLength,
  P: PartitionTimetreeAll<N, E> + ?Sized,
{
  partitions
    .iter()
    .map(|partition| partition.read_arc().create_edge_contribution(edge_key))
    .collect()
}

pub(super) fn create_branch_distributions_input_mode<N, E>(
  graph: &Graph<N, E, ()>,
  clock_rate: f64,
) -> Result<(), Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + HasBranchLength + TimetreeEdge,
{
  graph.get_edges().par_iter().for_each(|edge_ref| {
    let mut edge = edge_ref.write_arc().payload().write_arc();

    if let Some(branch_length) = edge.branch_length() {
      // Convert branch length (substitutions/site) to time duration (years)
      // gamma > 1 means faster evolution, so same substitutions correspond to shorter time
      let effective_clock_rate = clock_rate * edge.gamma();
      let time_duration = branch_length / effective_clock_rate;
      let distribution = Distribution::point(time_duration, 1.0);
      edge.set_time_length(Some(time_duration));
      edge.set_branch_length_distribution(Some(Arc::new(distribution)));
    }
  });

  Ok(())
}
