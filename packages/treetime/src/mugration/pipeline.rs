use crate::ancestral::marginal::branch_lengths_or_zero;
use crate::cancel::Cancel;
use crate::constants::MIN_BRANCH_LENGTH_FRACTION;
use crate::error::OperationError;
use crate::gtr::gtr::{GTR, GTRParams};
use crate::gtr::refinement::refine_gtr_model_and_rate;
use crate::partition::marginal::discrete::partition::PartitionMarginalDiscrete;
use crate::partition::marginal::shared::update::{MarginalPasses, MarginalUpdate};
use crate::partition::storage::dense::DenseNodeState;
use crate::partition::storage::discrete::DiscreteStates;
use crate::{make_error, make_report};
use eyre::Report;
use indexmap::IndexSet;
use itertools::Itertools;
use log::{info, warn};
use ndarray::Array1;
use statrs::statistics::Statistics;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

/// Scientific and algorithm policy for discrete-trait mugration inference.
pub struct MugrationParams {
  /// Marker value for missing discrete attributes, excluded from the model alphabet.
  pub missing_data: String,
  /// Pseudo-count regularization for GTR inference. `None` uses the refinement default of `1.0`.
  pub pc: Option<f64>,
  /// Maximum tolerated ratio of observed attributes absent from the weights file.
  pub missing_weights_threshold: f64,
  /// Number of GTR re-estimation iterations.
  pub iterations: usize,
  /// Optional sampling-bias correction applied to the inferred rate.
  pub sampling_bias_correction: Option<f64>,
  /// When set, smooths the initial equilibrium prior with a pseudo-count before the first pass.
  pub smooth_initial_pi: bool,
  /// When set, drops uninformative (uniform-posterior) roots from the equilibrium-frequency estimate.
  pub filter_uninformative_root: bool,
}

/// Parsed domain input for mugration inference.
pub struct MugrationInput {
  /// Tree topology.
  pub graph: Graph,
  /// Observed discrete attribute per leaf, keyed by leaf name.
  pub traits: BTreeMap<String, String>,
  /// Optional per-state weights (equilibrium-frequency prior), keyed by state name.
  pub weights: Option<BTreeMap<String, f64>>,
  /// Raw per-edge branch lengths captured from the Newick parse, keyed by edge id.
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}

/// Aggregate result of mugration inference.
///
/// Holds the reconstructed discrete-trait value maps gathered from the pipeline-local partition
/// before it is dropped, alongside the inferred model. The output writers project every response
/// body, file, and wire form from these value maps; the partition never leaves the core.
#[derive(Debug)]
pub struct MugrationOutput {
  /// Tree topology.
  pub graph: Graph,
  /// The inferred discrete GTR model.
  pub gtr: GTR,
  /// Discrete state names in order.
  pub states: DiscreteStates,
  /// Number of real states (excludes the missing-data marker).
  pub n_states: usize,
  /// Reconstructed discrete trait per node (argmax state name), or `None` when the node has no profile.
  pub reconstructed_traits: BTreeMap<GraphNodeKey, Option<String>>,
  /// Confidence profile per node (raw, unfiltered), or `None` when the node has no profile.
  pub confidences: BTreeMap<GraphNodeKey, Option<Array1<f64>>>,
}

/// Result of the weights-coverage check: the observed attributes absent from the weights file.
#[derive(Debug)]
pub struct WeightCoverageResult {
  pub missing_values: IndexSet<String>,
  pub missing_ratio: f64,
}

pub fn run(
  params: &MugrationParams,
  input: MugrationInput,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  cancel: &dyn Cancel,
) -> Result<MugrationOutput, OperationError> {
  cancel.check()?;

  let MugrationInput {
    graph,
    traits,
    weights,
    branch_lengths,
  } = input;
  let weights = weights.as_ref();

  let observed_values: IndexSet<String> = traits.values().sorted().cloned().collect();

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

  let discrete_states = DiscreteStates::from_values(model_values.iter().map(String::as_str), &params.missing_data);
  let n_states = discrete_states.len();

  if n_states < 2 {
    return Err(OperationError::InvalidInput(make_report!(
      "Mugration: only {n_states} discrete attributes provided for mugration. At least 2 are required."
    )));
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
  let pi = if params.smooth_initial_pi {
    apply_pseudo_counts(pi, Some(params.pc.unwrap_or(1.0)))
  } else {
    pi
  };

  let gtr = GTR::new(GTRParams {
    n_states,
    mu: 1.0,
    W: None,
    pi,
  })?;

  let partition = PartitionMarginalDiscrete::new(
    discrete_states,
    MIN_BRANCH_LENGTH_FRACTION,
    params.filter_uninformative_root,
  );
  let node_states = partition.attach_traits(&graph, &traits, names)?;

  let profile_lengths = branch_lengths_or_zero(&branch_lengths);
  let update = partition.marginal_update(&gtr, &graph, &profile_lengths, node_states)?;
  info!("Mugration: initial log likelihood = {:.4}", update.log_lh.value());

  // The partition is an immutable source; refinement threads the model through as a value and returns
  // the refined model with its own reconstruction result maps. Mugration optimizes the rate.
  let (gtr, MarginalUpdate { node_states, .. }) = refine_gtr_model_and_rate(
    &partition,
    gtr,
    update,
    params.iterations,
    fixed_pi.as_ref(),
    params.pc.unwrap_or(1.0),
    params.sampling_bias_correction,
    &graph,
    &profile_lengths,
  )?;

  // Gather the reconstructed value maps off the pipeline-local partition and its node states before
  // they leave scope, taking the partition read out of the serialization path. The reads are keyed by
  // node key and independent of node ordering, so the maps stay bit-identical regardless of any later
  // topology ordering the adapter applies to the returned graph.
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

/// Gather the per-node reconstructed trait (argmax state name) and confidence profile off the
/// mugration discrete partition. Keyed over every node, so the adapter reads plain value maps instead
/// of the partition during serialization.
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

#[allow(clippy::as_conversions, reason = "count/index numeric cast is exact for the domain range")]
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

#[allow(clippy::as_conversions, reason = "count/index numeric cast is exact for the domain range")]
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
