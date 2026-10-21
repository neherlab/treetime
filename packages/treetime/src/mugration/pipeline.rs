use crate::alphabet::state_alphabet::StateAlphabet;
use crate::branch_lengths::branch_lengths_or_zero;
use crate::cancel::Cancel;
use crate::constants::MIN_BRANCH_LENGTH_FRACTION;
use crate::error::OperationError;
use crate::gtr::gtr::{GTR, Mu};
use crate::gtr::infer_gtr::common::InferGtrOptions;
use crate::gtr::refinement::{GtrRefinement, refine_gtr};
use crate::mugration::leaf_observations::leaf_observations;
use crate::partition::marginal::discrete::partition::PartitionMarginalDiscrete;
use crate::partition::marginal::discrete::state_posterior::StatePosterior;
use crate::partition::marginal::shared::update::MarginalUpdate;
use crate::progress::LogSink;
use crate::{make_error, progress_info, progress_warn};
use eyre::Report;
use indexmap::IndexSet;
use itertools::Itertools;
use ndarray::Array1;
use statrs::statistics::Statistics;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_utils::cast::ToF64Exact;

pub fn run(
  params: &MugrationParams,
  input: MugrationInput,
  names: &NodeNames,
  cancel: &dyn Cancel,
  log: &dyn LogSink,
) -> Result<MugrationOutput, OperationError> {
  cancel.check().map_err(OperationError::classify)?;

  let MugrationInput {
    graph,
    traits,
    traits_source,
    observed_values,
    weights,
    branch_lengths,
  } = input;
  let weights = weights.as_ref();

  let observed_values: IndexSet<String> = observed_values.into_iter().collect();

  let model_values: IndexSet<String> = match weights {
    Some(weights_map) => {
      let weights_keys: IndexSet<String> = weights_map.keys().sorted().cloned().collect();

      let coverage = validate_weight_coverage(
        &observed_values,
        &weights_keys,
        &params.missing_data,
        params.missing_weights_threshold,
      )
      .map_err(OperationError::InvalidInput)?;

      if !coverage.missing_values.is_empty() {
        progress_warn!(
          log,
          "Mugration: discrete attributes missing from weights file: {} (ratio: {:.3})",
          coverage.missing_values.iter().join(", "),
          coverage.missing_ratio
        );
      }

      observed_values.union(&weights_keys).cloned().collect()
    },
    None => observed_values,
  };

  let states = StateAlphabet::new(model_values.into_iter().filter(|value| *value != params.missing_data))
    .map_err(OperationError::InvalidInput)?;
  let n_states = states.len();

  progress_info!(
    log,
    "Mugration: found {n_states} discrete states: {}",
    states.iter().join(", ")
  );

  let pi = match weights {
    Some(weights_map) => compute_pi_from_weights(&states, weights_map),
    None => compute_pi_uniform(n_states),
  };

  let fixed_pi = weights.map(|_| pi.clone());

  let pi = if params.smooth_initial_pi {
    apply_pseudo_counts(pi, params.pc)
  } else {
    pi
  };

  let gtr = GTR::builder()
    .n_states(n_states)
    .mu(Mu::FromMatrix(1.0))
    .pi(pi)
    .build()
    .map_err(OperationError::classify)?;

  let observations = leaf_observations(
    &graph,
    &traits,
    names,
    &states,
    &params.missing_data,
    &traits_source,
    log,
  )
  .map_err(OperationError::classify)?;
  let partition = PartitionMarginalDiscrete::new(
    states,
    &graph,
    &observations,
    MIN_BRANCH_LENGTH_FRACTION,
    params.filter_uninformative_root,
  )
  .map_err(OperationError::classify)?;

  let branch_lengths = branch_lengths_or_zero(&branch_lengths);
  let refinement = GtrRefinement {
    iterations: params.iterations,
    inference: InferGtrOptions {
      fixed_pi,
      pc: params.pc,
      ..InferGtrOptions::default()
    },
    sampling_bias_correction: params.sampling_bias_correction,
  };
  let (gtr, MarginalUpdate { node_states, .. }) =
    refine_gtr(&partition, gtr, &(), &refinement, &graph, &branch_lengths, log).map_err(OperationError::classify)?;

  let posteriors = partition
    .state_posteriors(&graph, &node_states)
    .map_err(OperationError::classify)?;

  Ok(MugrationOutput {
    graph,
    gtr,
    states: partition.states,
    posteriors,
  })
}

pub struct MugrationParams {
  pub missing_data: String,
  pub pc: f64,
  pub missing_weights_threshold: f64,
  pub iterations: usize,
  pub sampling_bias_correction: Option<f64>,
  pub smooth_initial_pi: bool,
  pub filter_uninformative_root: bool,
}

pub struct MugrationInput {
  pub graph: Graph,
  pub traits: BTreeMap<GraphNodeKey, String>,
  pub traits_source: String,
  pub observed_values: BTreeSet<String>,
  pub weights: Option<BTreeMap<String, f64>>,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}

#[derive(Debug)]
pub struct MugrationOutput {
  pub graph: Graph,
  pub gtr: GTR,
  pub states: StateAlphabet,
  pub posteriors: BTreeMap<GraphNodeKey, StatePosterior>,
}

pub(crate) fn validate_weight_coverage(
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

  let missing_ratio = missing_values.len().to_f64_exact() / unique_values.len().max(1).to_f64_exact();

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

#[derive(Debug)]
pub struct WeightCoverageResult {
  pub(crate) missing_values: IndexSet<String>,
  pub(crate) missing_ratio: f64,
}

pub(crate) fn compute_pi_from_weights(states: &StateAlphabet, weights: &BTreeMap<String, f64>) -> Array1<f64> {
  let mean_weight = weights.values().mean();

  let weights_arr: Array1<f64> = states
    .iter()
    .map(|state| *weights.get(state).unwrap_or(&mean_weight))
    .collect();

  let sum = weights_arr.sum();
  weights_arr / sum
}

pub(crate) fn compute_pi_uniform(n_states: usize) -> Array1<f64> {
  Array1::from_elem(n_states, 1.0 / n_states.to_f64_exact())
}

pub(crate) fn apply_pseudo_counts(pi: Array1<f64>, pc: f64) -> Array1<f64> {
  let pi = pi + pc;
  let sum = pi.sum();
  pi / sum
}
