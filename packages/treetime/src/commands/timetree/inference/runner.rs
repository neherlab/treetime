use crate::commands::clock::clock_regression::{ClockOptions, estimate_clock_model_with_reroot};
use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::commands::optimize::optimize_unified::OptimizationContribution;
use crate::commands::timetree::inference::branch_length_likelihood::compute_branch_length_distribution;
use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_convolution::distribution_convolution;
use crate::distribution::distribution_function::DistributionFunction;
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::edge::GraphEdgeKey;
use crate::graph::graph::GraphNodeForward;
use crate::representation::graph_ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use eyre::Report;
use parking_lot::RwLock;
use std::cmp::Ordering;
use std::sync::Arc;

const BRANCH_GRID_SIZE: usize = 200;

pub fn run_timetree(
  graph: &mut GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll>>],
  keep_root: bool,
) -> Result<(), Report> {
  log::info!("# Running timetree inference");

  if !partitions.is_empty() {
    log::info!("## Computing branch distributions from partitions");
    compute_branch_distributions(graph, partitions)?;
  } else {
    log::info!("## Estimating clock rate from root-to-tip regression (input branch mode)");
    let clock_rate = estimate_clock_rate_from_root_to_tip(graph, keep_root)?;
    log::info!("**Estimated clock rate:** {clock_rate:.6e}");
    log::debug!("Clock rate: {clock_rate}");

    log::info!("## Creating branch distributions from input lengths");
    create_branch_distributions_from_input_lengths(graph, clock_rate)?;
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
) -> Result<(), Report> {
  // In input branch mode, partitions exist but have no edge data
  // Skip branch distribution computation and use input branch lengths directly
  if partitions.iter().any(|p| p.read_arc().get_sequence_length().is_none()) {
    return Ok(());
  }

  let one_mutation = calculate_one_mutation(partitions);

  for edge_ref in graph.get_edges() {
    let edge_key = edge_ref.read_arc().key();
    let mut edge = edge_ref.write_arc().payload().write_arc();
    let branch_length = edge.branch_length.unwrap_or(one_mutation);

    let contributions = collect_contributions(partitions, edge_key)?;
    let distribution =
      compute_branch_length_distribution(&contributions, branch_length, one_mutation, BRANCH_GRID_SIZE)?;
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

fn estimate_clock_rate_from_root_to_tip(graph: &mut GraphAncestral, keep_root: bool) -> Result<f64, Report> {
  log::debug!("Estimating clock rate using root-to-tip regression");
  let clock_model = estimate_clock_model_with_reroot(
    graph,
    &ClockOptions::default(),
    keep_root,
    &BranchPointOptimizationParams::default(),
  )?;
  log::debug!(
    "Clock model estimated: rate={:.6e}, intercept={:.4}, r_val={:.4}",
    clock_model.clock_rate(),
    clock_model.intercept(),
    clock_model.r_val()
  );
  Ok(clock_model.clock_rate())
}

fn create_branch_distributions_from_input_lengths(graph: &GraphAncestral, clock_rate: f64) -> Result<(), Report> {
  for edge_ref in graph.get_edges() {
    let mut edge = edge_ref.write_arc().payload().write_arc();

    if let Some(branch_length) = edge.branch_length {
      // Convert branch length (substitutions/site) to time duration (years)
      let time_duration = branch_length / clock_rate;

      // Create a delta distribution (point mass) at the time duration
      let distribution = Distribution::point(time_duration, 1.0);
      edge.branch_length_distribution = Some(Arc::new(distribution));
    }
  }

  Ok(())
}

fn propagate_distributions_backward(graph: &GraphAncestral) -> Result<(), Report> {
  graph.par_iter_breadth_first_backward(|mut node| {
    if !node.is_leaf {
      let mut result: Option<Distribution> = None;

      for (child, edge) in &node.children {
        let child = child.read_arc();
        let edge = edge.read_arc();

        if let (Some(branch_dist), Some(child_time_dist)) = (&edge.branch_length_distribution, &child.time_distribution)
        {
          let convolved = convolve_safe(branch_dist, child_time_dist);
          result = Some(multiply_distributions(result, convolved));
        }
      }

      if let Some(dist) = result {
        node.payload.time_distribution = Some(Arc::new(dist));
      }
    }
    GraphTraversalContinuation::Continue
  });
  Ok(())
}

fn convolve_safe(a: &Distribution, b: &Distribution) -> Distribution {
  distribution_convolution(a, b).unwrap_or_else(|_| Distribution::empty())
}

fn multiply_distributions(existing: Option<Distribution>, new: Distribution) -> Distribution {
  if let Some(existing) = existing {
    convolve_safe(&existing, &new)
  } else {
    new
  }
}

fn propagate_distributions_forward(graph: &GraphAncestral) -> Result<(), Report> {
  graph.par_iter_breadth_first_forward(|mut node| {
    extract_time_estimate(&mut node);
    refine_distribution_from_parent(&mut node);
    GraphTraversalContinuation::Continue
  });
  Ok(())
}

fn extract_time_estimate(node: &mut GraphNodeForward<NodeAncestral, EdgeAncestral, ()>) {
  if let Some(time_dist) = &node.payload.time_distribution {
    let time_estimate = match time_dist.as_ref() {
      Distribution::Point(p) => p.t(),
      Distribution::Range(r) => f64::midpoint(r.start(), r.end()),
      Distribution::Function(f) => find_distribution_peak(f),
      Distribution::Empty => 0.0,
    };
    node.payload.time = Some(time_estimate);
  }
}

fn find_distribution_peak(f: &DistributionFunction<f64>) -> f64 {
  let max_idx = f
    .y()
    .iter()
    .enumerate()
    .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(Ordering::Equal))
    .map_or(0, |(idx, _)| idx);
  f.t()[max_idx]
}

fn refine_distribution_from_parent(node: &mut GraphNodeForward<NodeAncestral, EdgeAncestral, ()>) {
  if node.is_root {
    return;
  }

  if let Ok((parent, edge)) = node.get_exactly_one_parent() {
    let parent = parent.read_arc();
    let edge = edge.read_arc();

    if let (Some(parent_time_dist), Some(branch_dist)) = (&parent.time_distribution, &edge.branch_length_distribution) {
      let dist_from_parent = convolve_safe(parent_time_dist, branch_dist);

      let combined = if let Some(subtree_dist) = &node.payload.time_distribution {
        convolve_safe(&dist_from_parent, subtree_dist)
      } else {
        dist_from_parent
      };

      node.payload.time_distribution = Some(Arc::new(combined));
    }
  }
}
