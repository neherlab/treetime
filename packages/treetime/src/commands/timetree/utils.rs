use crate::commands::clock::clock_set::ClockSet;
use crate::commands::clock::clock_traits::ClockNode;
use crate::commands::timetree::timetree_traits::TimetreeNode;
use crate::seq::div::{OnlyLeaves, compute_divs};
use eyre::Report;
use ndarray::Array1;
use ordered_float::OrderedFloat;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::{Distribution, DistributionFunction};
use treetime_graph::edge::{BranchDistribution, EdgeOptimizeOps, GraphEdge, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};

pub fn initialize_node_divergences<N, E, D>(graph: &Graph<N, E, D>)
where
  N: GraphNode + Named + ClockNode,
  E: EdgeOptimizeOps,
  D: Send + Sync,
{
  let divs = compute_divs(graph, OnlyLeaves(false));
  for node_ref in graph.get_nodes() {
    let mut node = node_ref.write_arc().payload().write_arc();
    let name = node.name().map(|n| n.as_ref().to_owned());
    if let Some(name) = name {
      if let Some(&div) = divs.get(&name) {
        node.set_div(div);
      }
    }
  }
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
  E: GraphEdge + HasBranchLength + BranchDistribution<Arc<Distribution>>,
  D: Send + Sync,
{
  let seq_len_f64 = seq_len as f64;

  for edge_ref in graph.get_edges() {
    let mut edge = edge_ref.write_arc().payload().write_arc();

    if let Some(branch_length) = edge.branch_length() {
      let expected_time = branch_length / mu;
      let max_time = 3.0 * expected_time.max(1.0);
      let dx = max_time / (n_points - 1) as f64;

      let y = Array1::from_shape_fn(n_points, |i| {
        let dt = i as f64 * dx;
        if dt < 1e-10 {
          0.0
        } else {
          let log_p = -dt * mu * seq_len_f64 + branch_length * seq_len_f64 * (dt * mu * seq_len_f64).ln();
          log_p.exp()
        }
      });

      let y_max = y.iter().copied().map(OrderedFloat).max().map_or(1.0, |x| x.0);
      let y_normalized = y.mapv(|v| v / y_max);

      let distribution_fn = DistributionFunction::from_start_dx_values(0.0, dx, y_normalized)?;
      let distribution = Distribution::Function(distribution_fn);
      edge.set_branch_length_distribution(Some(Arc::new(distribution)));
    }
  }

  Ok(())
}
