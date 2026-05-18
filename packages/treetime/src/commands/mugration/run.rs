use crate::ancestral::marginal::update_marginal_mut;
use crate::commands::mugration::args::TreetimeMugrationArgs;
use crate::commands::mugration::comment_provider::PartitionCommentProvider;
use crate::commands::mugration::input::MugrationInput;
use crate::commands::mugration::output::{MugrationGtrOutput, MugrationResult, MugrationTraitsOutput};
use crate::constants::MIN_BRANCH_LENGTH_FRACTION;
use crate::gtr::gtr::{GTR, GTRParams};
use crate::gtr::refinement::refine_gtr_iterative;
use crate::partition::discrete_states::DiscreteStates;
use crate::partition::marginal_discrete::PartitionMarginalDiscrete;
use crate::partition::payload::ancestral::GraphAncestral;
use crate::{make_error, make_report};
use eyre::Report;
use indexmap::IndexSet;
use itertools::Itertools;
use log::{info, warn};
use ndarray::Array1;
use statrs::statistics::Statistics;
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;
use treetime_io::discrete_states_csv::read_discrete_attrs;
use treetime_io::nex::{NexWriteOptions, nex_write_str_with};
use treetime_io::nwk::CommentProviders;
use treetime_io::nwk::nwk_read_file;
use treetime_utils::io::json::{JsonPretty, json_write_file};

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

pub fn parse_mugration_input(args: &TreetimeMugrationArgs) -> Result<MugrationInput, Report> {
  let TreetimeMugrationArgs {
    tree,
    attribute,
    states,
    weights,
    name_column,
    pc,
    missing_data,
    missing_weights_threshold,
    ..
  } = args;

  let tree_path = tree.as_ref().ok_or_else(|| make_report!("Tree file is required"))?;
  let graph: GraphAncestral = nwk_read_file(tree_path)?;

  let (attr_values, _attr_name) =
    read_discrete_attrs::<String>(states, name_column, &Some(attribute.clone()), |s| Ok(s.to_owned()))?;
  let traits: BTreeMap<String, String> = attr_values.into_iter().collect();

  let weights_map = if let Some(weights_filepath) = weights {
    let (map, _) = read_discrete_attrs::<f64>(
      weights_filepath,
      &Some(attribute.clone()),
      &Some("weight".to_owned()),
      |s| Ok(s.parse::<f64>()?),
    )?;
    Some(map.into_iter().collect())
  } else {
    None
  };

  Ok(MugrationInput {
    graph,
    traits,
    attribute: attribute.clone(),
    weights: weights_map,
    missing_data: missing_data.clone(),
    pc: *pc,
    missing_weights_threshold: *missing_weights_threshold,
    iterations: args.iterations,
    sampling_bias_correction: args.sampling_bias_correction,
  })
}

pub fn execute_mugration(input: MugrationInput) -> Result<MugrationResult, Report> {
  let MugrationInput {
    graph,
    traits,
    attribute,
    weights,
    missing_data,
    pc,
    missing_weights_threshold,
    iterations,
    sampling_bias_correction,
  } = input;

  let observed_values: IndexSet<String> = traits.values().sorted().cloned().collect();

  let model_values: IndexSet<String> = match &weights {
    Some(weights_map) => {
      let weights_keys: IndexSet<String> = weights_map.keys().sorted().cloned().collect();

      let coverage = validate_weight_coverage(
        &observed_values,
        &weights_keys,
        &missing_data,
        missing_weights_threshold,
      )?;

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

  let discrete_states = DiscreteStates::from_values(model_values.iter().map(String::as_str), &missing_data);
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

  let pi = match &weights {
    Some(weights_map) => compute_pi_from_weights(&discrete_states, weights_map),
    None => compute_pi_uniform(n_states),
  };

  let fixed_pi = weights.as_ref().map(|_| pi.clone());

  let pi = apply_pseudo_counts(pi, pc);

  let gtr = GTR::new(GTRParams {
    n_states,
    mu: 1.0,
    W: None,
    pi,
  })?;

  let mut partition = PartitionMarginalDiscrete::new(gtr, discrete_states, MIN_BRANCH_LENGTH_FRACTION);
  partition.attach_traits(&graph, &traits)?;

  let log_lh = update_marginal_mut(&graph, &mut partition)?;
  info!("Mugration: initial log likelihood = {log_lh:.4}");

  let log_lh = refine_gtr_iterative(
    &graph,
    &mut partition,
    iterations,
    fixed_pi.as_ref(),
    pc.unwrap_or(1.0),
    sampling_bias_correction,
  )?;

  Ok(MugrationResult::new(graph, partition, &attribute, log_lh))
}

pub fn run_mugration(mugration_args: &TreetimeMugrationArgs) -> Result<(), Report> {
  let outdir = &mugration_args.outdir;
  fs::create_dir_all(outdir)?;

  let input = parse_mugration_input(mugration_args)?;
  let result = execute_mugration(input)?;

  write_annotated_tree(&result.graph, &result.partition, &result.traits, outdir)?;
  write_gtr_json_file(&result.gtr, outdir)?;

  if let Some(confidence_path) = &mugration_args.confidence {
    write_confidence_csv(&result, confidence_path)?;
  }

  info!("Mugration: wrote output to {}", outdir.display());
  Ok(())
}

fn write_annotated_tree(
  graph: &GraphAncestral,
  partition: &PartitionMarginalDiscrete,
  traits: &MugrationTraitsOutput,
  outdir: &Path,
) -> Result<(), Report> {
  let provider = PartitionCommentProvider::new(partition, &traits.attribute);
  let providers = CommentProviders::new().with(&provider);
  let nexus = nex_write_str_with(graph, &NexWriteOptions::default(), &providers)?;
  fs::write(outdir.join("annotated_tree.nexus"), format!("{nexus}\n"))?;

  fs::write(outdir.join("traits.csv"), traits.render_csv())?;

  Ok(())
}

fn write_gtr_json_file(gtr: &MugrationGtrOutput, outdir: &Path) -> Result<(), Report> {
  json_write_file(outdir.join("gtr.json"), gtr, JsonPretty(true))
}

fn write_confidence_csv(result: &MugrationResult, output_path: &Path) -> Result<(), Report> {
  fs::write(output_path, result.confidence.render_csv())?;
  Ok(())
}
