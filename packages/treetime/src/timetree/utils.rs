use crate::clock::clock_state::ClockState;
use crate::seq::div::{OnlyLeaves, compute_divs};
use crate::timetree::timetree_state::TimetreeState;
use eyre::Report;
use ndarray::Array1;
use ordered_float::OrderedFloat;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::{Distribution, DistributionFunction, NegLog};
use treetime_graph::edge::{GraphEdge, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey};

/// Grid floor as a fraction of one mutation's worth of time. Keeps the first grid point strictly
/// above the hard boundary at `t = 0`, so the divergent `-ln p` there is never stored on the grid.
const MIN_TIME_MUTATION_FRACTION: f64 = 0.01;

/// Compute each node's cumulative root divergence and store it in the clock state value.
///
/// The divergence is a durable clock input carried between passes. It lives only in the threaded
/// [`ClockState`], not on the node payload; a node absent from the state (introduced by a topology
/// change since the last rebuild) is inserted with default fields before its divergence is written,
/// so a fresh polytomy or reroot node gets its divergence here rather than a stale zero.
pub fn initialize_node_divergences<N, E, D>(
  graph: &Graph<N, E, D>,
  clock_state: &mut ClockState,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Send + Sync,
{
  let divs = compute_divs(graph, OnlyLeaves(false), branch_lengths, names)?;
  for node_ref in graph.get_nodes() {
    let node = node_ref.read_arc();
    let key = node.key();
    if let Some(name) = &names[&key] {
      if let Some(&div) = divs.get(name) {
        clock_state.nodes.entry(key).or_default().div = div;
      }
    }
  }
  Ok(())
}

pub fn extract_node_times<N, E, D>(
  graph: &Graph<N, E, D>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  state: &TimetreeState,
) -> BTreeMap<String, f64>
where
  N: GraphNode,
  E: GraphEdge,
  D: Send + Sync,
{
  graph
    .get_nodes()
    .into_iter()
    .filter_map(|node_ref| {
      let key = node_ref.read_arc().key();
      let name = names[&key].clone()?;
      let time = state.nodes.get(&key).and_then(|node| node.time)?;
      Some((name, time))
    })
    .collect()
}

/// Build the Poisson branch-length distribution for each edge, keyed by edge.
///
/// Replicates v0 Python TreeTime's Poisson branch-length distribution:
/// P(dt) ~ exp(-dt * mu * L) * (dt * mu * L)^(b * L), where:
/// - `mu` = clock rate (substitutions/site/year)
/// - `L` = sequence length
/// - `b` = branch length (substitutions/site)
///
/// An edge with no branch length is absent from the returned map.
pub fn create_poisson_branch_distributions<N, E, D>(
  graph: &Graph<N, E, D>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  mu: f64,
  seq_len: usize,
  n_points: usize,
) -> Result<BTreeMap<GraphEdgeKey, Arc<Distribution<NegLog>>>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Send + Sync,
{
  let seq_len_f64 = seq_len as f64;

  let mut distributions = BTreeMap::new();
  for edge_ref in graph.get_edges() {
    let edge_key = edge_ref.read_arc().key();

    if let Some(branch_length) = branch_lengths[&edge_key] {
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
      distributions.insert(edge_key, Arc::new(distribution));
    }
  }

  Ok(distributions)
}
