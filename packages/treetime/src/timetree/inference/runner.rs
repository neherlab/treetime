use crate::clock::clock_model::ClockModel;
use crate::clock::clock_state::ClockState;
use crate::clock::date_constraints::DateConstraints;
use crate::coalescent::coalescent::CoalescentModel;
use crate::optimize::gather::{
  gather_timetree_edge_contributions, gather_timetree_edge_indel_counts, timetree_total_sequence_length,
};
use crate::optimize::indel::estimate_indel_rate;
use crate::partition::timetree::partition::PartitionTimetree;
use crate::timetree::inference::backward_pass::propagate_distributions_backward;
use crate::timetree::inference::branch_length_likelihood::compute_branch_length_distribution;
use crate::timetree::inference::forward_pass::propagate_distributions_forward;
use crate::timetree::timetree_state::TimetreeState;
use crate::timetree::utils::initialize_node_divergences;
use eyre::Report;
use log::{debug, info, warn};
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::{Distribution, NegLog};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

pub(crate) const GRID_POINTS: usize = 300;

pub(crate) const EPS: f64 = 5e-4;

#[expect(
  clippy::too_many_arguments,
  reason = "each argument is an independent input of this step; a parameter struct would be built only for this call"
)]
pub(crate) fn run_timetree(
  graph: &Graph,
  constraints: &DateConstraints,
  partitions: &[PartitionTimetree],
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  clock_model: &ClockModel,
  coalescent: Option<&CoalescentModel>,
  no_indels: bool,
  mut state: TimetreeState,
  clock_state: &mut ClockState,
) -> Result<TimetreeState, Report> {
  info!("# Running timetree inference");

  info!("## Calculating divergence distances");
  initialize_node_divergences(graph, clock_state, branch_lengths, names)?;

  state.reseed_from_values(graph);

  info!("## Using clock model");
  let clock_rate = clock_model.clock_rate();
  info!("**Clock rate:** {clock_rate:.6e}");

  if !partitions.is_empty() {
    info!("## Computing branch distributions from partitions");
    compute_branch_distributions_marginal_mode(graph, partitions, branch_lengths, clock_rate, no_indels, &mut state)?;
  } else {
    info!("## Creating branch distributions from input lengths");
    create_branch_distributions_input_mode(graph, branch_lengths, clock_rate, &mut state)?;
  }

  info!("## Propagating distributions backward");
  propagate_distributions_backward(graph, constraints, coalescent, &mut state)?;

  info!("## Propagating distributions forward");
  propagate_distributions_forward(graph, constraints, names, &mut state)?;

  info!("# Timetree inference completed");
  Ok(state)
}

pub(crate) const CLOCK_BRANCH_LENGTH_DAMPING: f64 = 0.5;

pub(crate) fn commit_clock_branch_lengths(
  graph: &Graph,
  clock_rate: f64,
  damping: f64,
  clock_branch_lengths: &mut BTreeMap<GraphEdgeKey, f64>,
  state: &TimetreeState,
) {
  let node_time = |key| state.nodes.get(&key).and_then(|node| node.time);

  let previous_lengths: &BTreeMap<GraphEdgeKey, f64> = clock_branch_lengths;
  let committed: Vec<(GraphEdgeKey, f64, bool)> = graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .filter_map(|edge_ref| {
      let key = edge_ref.key();
      let (Some(parent_time), Some(child_time)) = (node_time(edge_ref.source()), node_time(edge_ref.target())) else {
        return None;
      };

      let duration = child_time - parent_time;
      let fresh = clock_rate * state.edge(key).gamma * duration.max(0.0);
      let value = match previous_lengths.get(&key) {
        Some(previous) => (1.0 - damping) * previous + damping * fresh,
        None => fresh,
      };

      Some((key, value, duration < 0.0))
    })
    .collect();

  let mut inverted = 0_usize;
  for (key, value, is_inverted) in committed {
    clock_branch_lengths.insert(key, value);
    inverted += usize::from(is_inverted);
  }

  if inverted > 0 {
    warn!(
      "Timetree: {inverted} branch(es) run backwards in time, i.e. the child is dated before its \
       parent. Their clock branch lengths were committed as zero. This is expected only where an \
       observed leaf date conflicts with the fitted clock, since the forward pass clamps internal \
       nodes to their parent but leaves leaf dates as given."
    );
  }
}

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
fn compute_branch_distributions_marginal_mode(
  graph: &Graph,
  partitions: &[PartitionTimetree],
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  clock_rate: f64,
  no_indels: bool,
  state: &mut TimetreeState,
) -> Result<(), Report> {
  let total_sites = timetree_total_sequence_length(partitions);
  let one_mutation = 1.0 / total_sites as f64;

  let indel_counts = gather_timetree_edge_indel_counts(graph, partitions);
  let indel_rate = if no_indels {
    0.0
  } else {
    estimate_indel_rate(graph, &indel_counts, branch_lengths)
  };
  let contributions = gather_timetree_edge_contributions(graph, partitions)?;

  info!(
    "Computing branch distributions from {} partition(s) with {} total sites",
    partitions.len(),
    total_sites
  );
  debug!("One mutation = {one_mutation:.6e} substitutions/site");
  debug!("Indel rate = {indel_rate:.6e} indels/(site*time)");

  let edge_states: &TimetreeState = state;
  let distributions: Vec<(GraphEdgeKey, Option<f64>, Arc<Distribution<NegLog>>)> = graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(
      |edge_ref| -> Result<(GraphEdgeKey, Option<f64>, Arc<Distribution<NegLog>>), Report> {
        let edge_key = edge_ref.key();
        let branch_length = branch_lengths[&edge_key].unwrap_or(one_mutation);
        let gamma = edge_states.edge(edge_key).gamma;

        debug!("Edge {edge_key:?}: input branch_length = {branch_length:.6e}, gamma = {gamma:.4}");

        let contributions = &contributions[&edge_key];
        let indel_count: usize = if no_indels { 0 } else { indel_counts[&edge_key] };
        let distribution = compute_branch_length_distribution(
          contributions,
          indel_count,
          indel_rate,
          branch_length,
          one_mutation,
          GRID_POINTS,
          clock_rate,
          gamma,
        )?;

        let time_length = distribution.likely_time();
        if let Some(likely_time) = time_length {
          debug!("Edge {edge_key:?}: distribution peak at time = {likely_time:.6e}");
        }

        Ok((edge_key, time_length, distribution))
      },
    )
    .collect::<Result<Vec<_>, Report>>()?;

  for (key, time_length, distribution) in distributions {
    let entry = state.edge_mut(key);
    entry.time_length = time_length;
    entry.branch_length_distribution = Some(distribution);
  }
  Ok(())
}

pub(super) fn create_branch_distributions_input_mode(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  clock_rate: f64,
  state: &mut TimetreeState,
) -> Result<(), Report> {
  let edge_states: &TimetreeState = state;
  let distributions: Vec<(GraphEdgeKey, f64, Arc<Distribution<NegLog>>)> = graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .filter_map(|edge_ref| {
      let key = edge_ref.key();
      let time_duration = if let Some(branch_length) = branch_lengths[&key] {
        let effective_clock_rate = clock_rate * edge_states.edge(key).gamma;
        Some(branch_length / effective_clock_rate)
      } else {
        edge_states.edge(key).time_length
      };

      time_duration.map(|time_duration| {
        let distribution = Arc::new(Distribution::point(time_duration, 0.0));
        (key, time_duration, distribution)
      })
    })
    .collect();

  for (key, time_duration, distribution) in distributions {
    let entry = state.edge_mut(key);
    entry.time_length = Some(time_duration);
    entry.branch_length_distribution = Some(distribution);
  }

  Ok(())
}

pub(crate) fn timetree_branch_lengths(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  clock_branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> BTreeMap<GraphEdgeKey, f64> {
  graph
    .get_edges()
    .map(|edge| {
      let key = edge.key();
      let branch_length = clock_branch_lengths
        .get(&key)
        .copied()
        .or_else(|| branch_lengths[&key])
        .unwrap_or(0.0);
      (key, branch_length)
    })
    .collect()
}
