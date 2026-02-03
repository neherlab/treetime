use crate::commands::clock::clock_model::ClockModel;
use crate::commands::clock::clock_traits::ClockNode;
use crate::commands::optimize::optimize_unified::OptimizationContribution;
use crate::commands::timetree::coalescent::coalescent::compute_coalescent_contributions;
use crate::commands::timetree::inference::backward_pass::propagate_distributions_backward;
use crate::commands::timetree::inference::branch_length_likelihood::compute_branch_length_distribution;
use crate::commands::timetree::inference::forward_pass::propagate_distributions_forward;
use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::commands::timetree::timetree_traits::{TimetreeEdge, TimetreeNode};
use crate::commands::timetree::utils::initialize_node_divergences;
use crate::distribution::distribution::Distribution;
use crate::graph::edge::{GraphEdge, GraphEdgeKey, Weighted};
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, Named};
use eyre::Report;
use log::{debug, info};
use parking_lot::RwLock;
use std::sync::Arc;

pub const BRANCH_GRID_SIZE: usize = 200;

pub fn run_timetree<N, E, P>(
  graph: &mut Graph<N, E, ()>,
  partitions: &[Arc<RwLock<P>>],
  clock_model: &ClockModel,
  coalescent_tc: Option<f64>,
) -> Result<(), Report>
where
  N: GraphNode + Named + TimetreeNode + ClockNode,
  E: GraphEdge + Weighted + TimetreeEdge,
  P: PartitionTimetreeAll<N, E> + ?Sized,
{
  info!("# Running timetree inference");

  info!("## Calculating divergence distances");
  initialize_node_divergences(graph);

  info!("## Using clock model");
  let clock_rate = clock_model.clock_rate();
  info!("**Clock rate:** {clock_rate:.6e}");

  if !partitions.is_empty() {
    info!("## Computing branch distributions from partitions");
    compute_branch_distributions_marginal_mode(graph, partitions, clock_rate)?;
  } else {
    info!("## Creating branch distributions from input lengths");
    create_branch_distributions_input_mode(graph, clock_rate)?;
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

fn compute_branch_distributions_marginal_mode<N, E, P>(
  graph: &Graph<N, E, ()>,
  partitions: &[Arc<RwLock<P>>],
  clock_rate: f64,
) -> Result<(), Report>
where
  N: GraphNode + Named + TimetreeNode,
  E: GraphEdge + Weighted + TimetreeEdge,
  P: PartitionTimetreeAll<N, E> + ?Sized,
{
  let one_mutation = calculate_one_mutation(partitions);
  let total_sites: usize = partitions
    .iter()
    .map(|p| p.read_arc().get_sequence_length())
    .sum();

  info!(
    "Computing branch distributions from {} partition(s) with {} total sites",
    partitions.len(),
    total_sites
  );
  debug!("One mutation = {one_mutation:.6e} substitutions/site");

  for edge_ref in graph.get_edges() {
    let edge_key = edge_ref.read_arc().key();
    let mut edge = edge_ref.write_arc().payload().write_arc();
    let branch_length = edge.weight().unwrap_or(one_mutation);

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

    edge.set_time_length(distribution.likely_time());
    edge.set_branch_length_distribution(Some(distribution));
  }
  Ok(())
}

fn calculate_one_mutation<N, E, P>(partitions: &[Arc<RwLock<P>>]) -> f64
where
  N: GraphNode + Named,
  E: GraphEdge + Weighted,
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
  E: GraphEdge + Weighted,
  P: PartitionTimetreeAll<N, E> + ?Sized,
{
  partitions
    .iter()
    .map(|partition| partition.read_arc().create_edge_contribution(edge_key))
    .collect()
}

fn create_branch_distributions_input_mode<N, E>(graph: &Graph<N, E, ()>, clock_rate: f64) -> Result<(), Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + Weighted + TimetreeEdge,
{
  for edge_ref in graph.get_edges() {
    let mut edge = edge_ref.write_arc().payload().write_arc();

    if let Some(branch_length) = edge.weight() {
      // Convert branch length (substitutions/site) to time duration (years)
      let time_duration = branch_length / clock_rate;
      let distribution = Distribution::point(time_duration, 1.0);
      edge.set_time_length(Some(time_duration));
      edge.set_branch_length_distribution(Some(Arc::new(distribution)));
    }
  }

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::graph::edge::TimeLength;
  use crate::io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};
  use crate::representation::edge_timetree::EdgeTimetree;
  use crate::representation::node_timetree::NodeTimetree;
  use approx::assert_abs_diff_eq;

  #[test]
  fn test_create_branch_distributions_input_mode_sets_time_length() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:0.003,B:0.006)AB:0.009,C:0.012)root;")?;
    let clock_rate = 0.001; // 0.001 subs/site/year

    create_branch_distributions_input_mode(&graph, clock_rate)?;

    // Verify each edge has time_length = branch_length / clock_rate
    for edge_ref in graph.get_edges() {
      let edge = edge_ref.read_arc().payload().read_arc();
      let branch_length = edge.weight();
      let time_length = edge.time_length();

      if let Some(bl) = branch_length {
        let expected_time = bl / clock_rate;
        assert_abs_diff_eq!(time_length.unwrap_or(f64::NAN), expected_time, epsilon = 1e-10);
      }
    }

    Ok(())
  }

  #[test]
  fn test_input_mode_newick_output_uses_time_lengths() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:0.003,B:0.006)AB:0.009,C:0.012)root;")?;
    let clock_rate = 0.001;

    create_branch_distributions_input_mode(&graph, clock_rate)?;

    // EdgeTimetree.nwk_weight() returns time_length, so Newick output should show time values
    let newick_output = nwk_write_str(&graph, &NwkWriteOptions::default())?;

    // After conversion: 0.003/0.001=3, 0.006/0.001=6, 0.009/0.001=9, 0.012/0.001=12
    // Use patterns with terminators (comma or close paren) to avoid false matches like :3 matching :30
    assert!(
      newick_output.contains(":3,"),
      "Expected time length 3 in output: {newick_output}"
    );
    assert!(
      newick_output.contains(":6)"),
      "Expected time length 6 in output: {newick_output}"
    );
    assert!(
      newick_output.contains(":9,"),
      "Expected time length 9 in output: {newick_output}"
    );
    assert!(
      newick_output.contains(":12)"),
      "Expected time length 12 in output: {newick_output}"
    );

    Ok(())
  }
}
