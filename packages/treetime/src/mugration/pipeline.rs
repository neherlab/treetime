use crate::alphabet::state_alphabet::StateAlphabet;
use crate::ancestral::plan::fit_marginal;
use crate::branch_lengths::branch_lengths_or_zero;
use crate::cancel::Cancel;
use crate::constants::MIN_BRANCH_LENGTH_FRACTION;
use crate::error::OperationError;
use crate::gtr::gtr::{GTR, Mu};
use crate::gtr::refinement::{GtrRefinement, RefinementRate};
use crate::mugration::leaf_observations::leaf_observations;
use crate::partition::marginal::discrete::partition::PartitionMarginalDiscrete;
use crate::partition::marginal::shared::update::MarginalUpdate;
use crate::partition::storage::dense::DenseNodeState;
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
use treetime_utils::cast::ToF64Exact;

pub fn run(
  params: &MugrationParams,
  input: MugrationInput,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
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
    apply_pseudo_counts(pi, Some(params.pc.unwrap_or(1.0)))
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
    pc: params.pc.unwrap_or(1.0),
    fixed_pi,
    rate: RefinementRate::Optimize {
      scale: params.sampling_bias_correction,
    },
  };
  let (gtr, MarginalUpdate { node_states, .. }) =
    fit_marginal(&partition, gtr, &(), Some(&refinement), &graph, &branch_lengths, log)
      .map_err(OperationError::classify)?;

  let (reconstructed_traits, confidences) = gather_reconstruction_maps(&graph, &partition, &node_states);

  Ok(MugrationOutput {
    graph,
    gtr,
    states: partition.states.clone(),
    n_states: partition.n_states(),
    reconstructed_traits,
    confidences,
  })
}

pub struct MugrationParams {
  pub missing_data: String,
  pub pc: Option<f64>,
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
  pub n_states: usize,
  pub reconstructed_traits: BTreeMap<GraphNodeKey, Option<String>>,
  pub confidences: BTreeMap<GraphNodeKey, Option<Array1<f64>>>,
}

fn gather_reconstruction_maps(
  graph: &Graph,
  partition: &PartitionMarginalDiscrete,
  node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
) -> (
  BTreeMap<GraphNodeKey, Option<String>>,
  BTreeMap<GraphNodeKey, Option<Array1<f64>>>,
) {
  let reconstructed_traits = graph
    .get_nodes()
    .map(|node| {
      let key = node.key();
      (key, partition.get_reconstructed_trait(node_states, key))
    })
    .collect();
  let confidences = graph
    .get_nodes()
    .map(|node| {
      let key = node.key();
      (key, partition.get_confidence(node_states, key))
    })
    .collect();
  (reconstructed_traits, confidences)
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

pub(crate) fn apply_pseudo_counts(pi: Array1<f64>, pc: Option<f64>) -> Array1<f64> {
  match pc {
    Some(pc_val) => {
      let pi = &pi + pc_val;
      let sum = pi.sum();
      pi / sum
    },
    None => pi,
  }
}
