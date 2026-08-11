use crate::payload::clock_set::ClockSet;
use crate::payload::traits::ClockNode;
use crate::payload::traits::TimetreeNode;
use crate::seq::div::{OnlyLeaves, compute_divs};
use eyre::Report;
use ndarray::Array1;
use ordered_float::OrderedFloat;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::{Distribution, DistributionFunction, NegLog};
use treetime_graph::edge::{BranchDistribution, EdgeOptimizeOps, GraphEdge, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};

/// Grid floor as a fraction of one mutation's worth of time. Keeps the first grid point strictly
/// above the hard boundary at `t = 0`, so the divergent `-ln p` there is never stored on the grid.
const MIN_TIME_MUTATION_FRACTION: f64 = 0.01;

pub fn initialize_node_divergences<N, E, D>(graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode + Named + ClockNode,
  E: EdgeOptimizeOps,
  D: Send + Sync,
{
  let divs = compute_divs(graph, OnlyLeaves(false))?;
  for node_ref in graph.get_nodes() {
    let mut node = node_ref.write_arc().payload().write_arc();
    let name = node.name().map(|n| n.as_ref().to_owned());
    if let Some(name) = name {
      if let Some(&div) = divs.get(&name) {
        node.set_div(div);
      }
    }
  }
  Ok(())
}

pub fn initialize_clock_totals_from_time_distributions<N, E, D>(graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode + TimetreeNode + ClockNode,
  E: GraphEdge,
  D: Send + Sync,
{
  for node_ref in graph.get_nodes() {
    let mut node = node_ref.write_arc().payload().write_arc();
    if let Some(dist_arc) = node.time_distribution() {
      if let Some(time) = dist_arc.likely_time() {
        *node.clock_set_mut() = ClockSet::leaf_contribution(Some(time));
      }
    }
  }

  Ok(())
}

pub fn extract_node_times<N, E, D>(graph: &Graph<N, E, D>) -> BTreeMap<String, f64>
where
  N: GraphNode + Named + TimetreeNode,
  E: GraphEdge,
  D: Send + Sync,
{
  graph
    .get_nodes()
    .into_iter()
    .filter_map(|node_ref| {
      let node = node_ref.read_arc();
      let payload = node.payload().read_arc();
      let name = payload.name()?.as_ref().to_owned();
      let time = payload.time()?;
      Some((name, time))
    })
    .collect()
}

/// Construct Poisson branch-length distributions on each edge.
///
/// Replicates v0 Python TreeTime's Poisson branch-length distribution:
/// P(dt) ~ exp(-dt * mu * L) * (dt * mu * L)^(b * L), where:
/// - `mu` = clock rate (substitutions/site/year)
/// - `L` = sequence length
/// - `b` = branch length (substitutions/site)
pub fn create_poisson_branch_distributions<N, E, D>(
  graph: &Graph<N, E, D>,
  mu: f64,
  seq_len: usize,
  n_points: usize,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength + BranchDistribution<Arc<Distribution<NegLog>>>,
  D: Send + Sync,
{
  let seq_len_f64 = seq_len as f64;

  for edge_ref in graph.get_edges() {
    let mut edge = edge_ref.write_arc().payload().write_arc();

    if let Some(branch_length) = edge.branch_length() {
      let expected_time = branch_length / mu;
      let max_time = 3.0 * expected_time.max(1.0);

      // Floor the grid strictly above the hard boundary at `t = 0`. The Poisson density
      // `p(t) ~ (t * mu * L)^(b * L)` vanishes as `t -> 0` for a branch with mutations, so `-ln p`
      // diverges there; gridding from zero would store `+inf`. Start the first grid point a small
      // fraction of one mutation's worth of time above zero instead, as
      // `compute_branch_length_distribution` does, so every stored ordinate is finite.
      let min_time = MIN_TIME_MUTATION_FRACTION / (mu * seq_len_f64);
      let grid = Array1::linspace(min_time, max_time, n_points);

      // Log-likelihood on the grid. Every point is strictly positive, so the density is finite.
      let log_p = grid.mapv(|dt| -dt * mu * seq_len_f64 + branch_length * seq_len_f64 * (dt * mu * seq_len_f64).ln());

      // Negative-log ordinates peak-normalized to `0`: `-ln(p / p_peak) = ln(p_peak) - ln(p)`.
      let log_p_max = log_p.iter().copied().map(OrderedFloat).max().map_or(0.0, |x| x.0);
      let neg_log = log_p.mapv(|value| log_p_max - value);

      let distribution_fn = DistributionFunction::from_range_values((min_time, max_time), neg_log)?;
      let distribution = Distribution::Function(distribution_fn);
      edge.set_branch_length_distribution(Some(Arc::new(distribution)));
    }
  }

  Ok(())
}
