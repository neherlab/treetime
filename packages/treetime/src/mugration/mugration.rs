use crate::ancestral::marginal::{marginal_update, profile_branch_lengths};
use crate::constants::MIN_BRANCH_LENGTH_FRACTION;
use crate::gtr::gtr::{GTR, GTRParams};
use crate::gtr::refinement::refine_gtr_iterative;
use crate::make_error;
use crate::mugration::result::{MugrationOutputMaps, MugrationResult, gather_mugration_output_maps};
use crate::partition::marginal::discrete::partition::PartitionMarginalDiscrete;
use crate::partition::storage::discrete::DiscreteStates;
use crate::payload::ancestral::GraphAncestral;
use eyre::Report;
use indexmap::IndexSet;
use itertools::Itertools;
use log::{info, warn};
use ndarray::Array1;
use statrs::statistics::Statistics;
use std::cell::RefCell;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;

#[derive(Debug)]
pub struct WeightCoverageResult {
  pub missing_values: IndexSet<String>,
  pub missing_ratio: f64,
}

pub fn validate_weight_coverage(
  unique_values: &IndexSet<String>,
  weights_keys: &IndexSet<String>,
  missing_data: &str,
  threshold: f64,
) -> Result<WeightCoverageResult, Report> {
  let missing_values: IndexSet<String> = unique_values
    .difference(weights_keys)
    .filter(|&value| value != missing_data)
    .cloned()
    .collect();

  let missing_ratio = missing_values.len() as f64 / unique_values.len().max(1) as f64;

  if missing_ratio > threshold {
    return make_error!(
      "Mugration: too many discrete attributes missing from the weights file. \
       The ratio of missing values {missing_ratio} is greater than the threshold {threshold}."
    );
  }

  Ok(WeightCoverageResult {
    missing_values,
    missing_ratio,
  })
}

pub fn compute_pi_from_weights(states: &DiscreteStates, weights: &BTreeMap<String, f64>) -> Array1<f64> {
  let mean_weight = weights.values().mean();

  let weights_arr: Array1<f64> = states
    .iter()
    .map(|state| *weights.get(state).unwrap_or(&mean_weight))
    .collect();

  let sum = weights_arr.sum();
  weights_arr / sum
}

pub fn compute_pi_uniform(n_states: usize) -> Array1<f64> {
  Array1::from_elem(n_states, 1.0 / n_states as f64)
}

pub fn apply_pseudo_counts(pi: Array1<f64>, pc: Option<f64>) -> Array1<f64> {
  match pc {
    Some(pc_val) => {
      let pi = &pi + pc_val;
      let sum = pi.sum();
      pi / sum
    },
    None => pi,
  }
}

pub fn execute_mugration(
  graph: GraphAncestral,
  confidences: &BTreeMap<GraphNodeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  traits: &BTreeMap<String, String>,
  attribute: &str,
  weights: Option<&BTreeMap<String, f64>>,
  missing_data: &str,
  pc: Option<f64>,
  missing_weights_threshold: f64,
  iterations: usize,
  sampling_bias_correction: Option<f64>,
  smooth_initial_pi: bool,
  filter_uninformative_root: bool,
) -> Result<(MugrationResult, MugrationOutputMaps), Report> {
  let observed_values: IndexSet<String> = traits.values().sorted().cloned().collect();

  let model_values: IndexSet<String> = match weights {
    Some(weights_map) => {
      let weights_keys: IndexSet<String> = weights_map.keys().sorted().cloned().collect();

      let coverage =
        validate_weight_coverage(&observed_values, &weights_keys, missing_data, missing_weights_threshold)?;

      if !coverage.missing_values.is_empty() {
        warn!(
          "Mugration: discrete attributes missing from weights file: {} (ratio: {:.3})",
          coverage.missing_values.iter().join(", "),
          coverage.missing_ratio
        );
      }

      observed_values.union(&weights_keys).cloned().collect()
    },
    None => observed_values,
  };

  let discrete_states = DiscreteStates::from_values(model_values.iter().map(String::as_str), missing_data);
  let n_states = discrete_states.len();

  if n_states < 2 {
    return make_error!(
      "Mugration: only {n_states} discrete attributes provided for mugration. At least 2 are required."
    );
  }

  info!(
    "Mugration: found {n_states} discrete states: {}",
    discrete_states.iter().join(", ")
  );

  let pi = match weights {
    Some(weights_map) => compute_pi_from_weights(&discrete_states, weights_map),
    None => compute_pi_uniform(n_states),
  };

  let fixed_pi = weights.map(|_| pi.clone());

  // v0 builds the initial GTR from the raw equilibrium frequencies and reserves
  // the pseudo-count for infer_gtr regularization. Smoothing the initial pi
  // (a flatter prior for the first reconstruction pass) is opt-in v1 behavior.
  // When enabled it uses the same effective pseudo-count as the refinement path
  // (`pc.unwrap_or(1.0)`), so the two pi-smoothing paths stay consistent.
  let pi = if smooth_initial_pi {
    apply_pseudo_counts(pi, Some(pc.unwrap_or(1.0)))
  } else {
    pi
  };

  let gtr = GTR::new(GTRParams {
    n_states,
    mu: 1.0,
    W: None,
    pi,
  })?;

  let mut partition = PartitionMarginalDiscrete::new(
    gtr,
    discrete_states,
    MIN_BRANCH_LENGTH_FRACTION,
    filter_uninformative_root,
  );
  partition.attach_traits(&graph, traits, names)?;

  let mut partition = RefCell::new(partition);

  let log_lh = marginal_update(
    &graph,
    &profile_branch_lengths(branch_lengths),
    std::slice::from_mut(partition.get_mut()),
  )?;
  info!("Mugration: initial log likelihood = {:.4}", log_lh.value());

  let log_lh = refine_gtr_iterative(
    &graph,
    &partition,
    branch_lengths,
    iterations,
    fixed_pi.as_ref(),
    pc.unwrap_or(1.0),
    sampling_bias_correction,
    true,
  )?;

  let partition = partition.into_inner();

  // Gather the output value maps off the pipeline-local partition before it enters the graph data
  // slot, taking the partition read out of the serialization path. The maps stay a local the command
  // threads to the writers; the partition is dropped at the end of this function.
  let maps = gather_mugration_output_maps(&graph, &partition);
  let result = MugrationResult::new(graph, confidences, names, branch_lengths, &partition, attribute, log_lh);
  Ok((result, maps))
}
