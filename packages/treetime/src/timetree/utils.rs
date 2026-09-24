use crate::clock::clock_state::ClockState;
use crate::seq::div::{OnlyLeaves, compute_divs};
use crate::timetree::timetree_state::TimetreeState;
use eyre::Report;
use ndarray::Array1;
use ordered_float::OrderedFloat;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::{Distribution, DistributionFunction, NegLog};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

const MIN_TIME_MUTATION_FRACTION: f64 = 0.01;

pub(crate) fn initialize_node_divergences(
  graph: &Graph,
  clock_state: &mut ClockState,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<(), Report> {
  let divs = compute_divs(graph, OnlyLeaves(false), branch_lengths, names)?;
  for node_ref in graph.get_nodes() {
    let node = node_ref;
    let key = node.key();
    if let Some(name) = &names[&key] {
      if let Some(&div) = divs.get(name) {
        clock_state.nodes.entry(key).or_default().div = div;
      }
    }
  }
  Ok(())
}

pub(crate) fn extract_node_times(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  state: &TimetreeState,
) -> BTreeMap<String, f64> {
  graph
    .get_nodes()
    .filter_map(|node_ref| {
      let key = node_ref.key();
      let name = names[&key].clone()?;
      let time = state.nodes.get(&key).and_then(|node| node.time)?;
      Some((name, time))
    })
    .collect()
}

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub(crate) fn create_poisson_branch_distributions(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  mu: f64,
  seq_len: usize,
  n_points: usize,
) -> Result<BTreeMap<GraphEdgeKey, Arc<Distribution<NegLog>>>, Report> {
  let seq_len_f64 = seq_len as f64;

  let mut distributions = BTreeMap::new();
  for edge_ref in graph.get_edges() {
    let edge_key = edge_ref.key();

    if let Some(branch_length) = branch_lengths[&edge_key] {
      let expected_time = branch_length / mu;
      let max_time = 3.0 * expected_time.max(1.0);

      let min_time = MIN_TIME_MUTATION_FRACTION / (mu * seq_len_f64);
      let grid = Array1::linspace(min_time, max_time, n_points);

      let log_p = grid.mapv(|dt| -dt * mu * seq_len_f64 + branch_length * seq_len_f64 * (dt * mu * seq_len_f64).ln());

      let log_p_max = log_p.iter().copied().map(OrderedFloat).max().map_or(0.0, |x| x.0);
      let neg_log = log_p.mapv(|value| log_p_max - value);

      let distribution_fn = DistributionFunction::from_range_values((min_time, max_time), neg_log)?;
      let distribution = Distribution::Function(distribution_fn);
      distributions.insert(edge_key, Arc::new(distribution));
    }
  }

  Ok(distributions)
}
