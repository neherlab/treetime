use crate::commands::mugration::args::TreetimeMugrationArgs;
use crate::commands::mugration::discrete_marginal::{attach_traits, run_discrete_marginal};
use crate::gtr::gtr::{GTR, GTRParams};
use crate::make_error;
use crate::representation::discrete_states::DiscreteStates;
use crate::representation::partition::discrete::PartitionDiscrete;
use crate::representation::payload::ancestral::GraphAncestral;
use eyre::Report;
use indexmap::IndexSet;
use itertools::Itertools;
use log::{info, warn};
use ndarray::Array1;
use serde::Serialize;
use statrs::statistics::Statistics;
use std::collections::BTreeMap;
use std::fmt::Write as FmtWrite;
use std::fs;
use std::io::Write;
use std::sync::Arc;
use treetime_graph::node::Named;
use treetime_io::discrete_states_csv::read_discrete_attrs;
use treetime_io::nwk::nwk_read_file;
use treetime_io::nwk::{EdgeToNwk, NodeToNwk, NwkWriteOptions, format_weight};

/// Result of validating weight coverage against unique values.
#[derive(Debug)]
pub struct WeightCoverageResult {
  /// Values present in unique_values but missing from weights.
  pub missing_values: IndexSet<String>,
  /// Ratio of missing values to total unique values.
  pub missing_ratio: f64,
}

/// Validate that weights cover enough of the unique values.
///
/// Returns the set of missing values and the missing ratio.
/// Returns an error if the missing ratio exceeds the threshold.
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

/// Compute equilibrium frequencies (pi) from a weights map.
///
/// For states not present in the weights map, uses the mean weight as fallback.
pub fn compute_pi_from_weights(states: &DiscreteStates, weights: &BTreeMap<String, f64>) -> Array1<f64> {
  let mean_weight = weights.values().mean();

  let weights_arr: Array1<f64> = states
    .iter()
    .map(|state| *weights.get(state).unwrap_or(&mean_weight))
    .collect();

  let sum = weights_arr.sum();
  weights_arr / sum
}

/// Compute uniform equilibrium frequencies (pi) for n states.
pub fn compute_pi_uniform(n_states: usize) -> Array1<f64> {
  Array1::from_elem(n_states, 1.0 / n_states as f64)
}

/// Apply pseudo-counts to equilibrium frequencies and re-normalize.
///
/// If `pc` is None, returns the input unchanged.
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

pub fn run_mugration(mugration_args: &TreetimeMugrationArgs) -> Result<(), Report> {
  let TreetimeMugrationArgs {
    tree,
    attribute,
    states,
    weights,
    name_column,
    confidence,
    pc,
    missing_data,
    missing_weights_threshold,
    outdir,
    ..
  } = mugration_args;

  fs::create_dir_all(outdir)?;

  // Read tree
  let tree_path = tree.as_ref().ok_or_else(|| eyre::eyre!("Tree file is required"))?;
  let graph: GraphAncestral = nwk_read_file(tree_path)?;

  // Read trait values
  let (attr_values, _attr_name) =
    read_discrete_attrs::<String>(states, name_column, &Some(attribute.clone()), |s| Ok(s.to_owned()))?;

  let unique_values: IndexSet<String> = attr_values.values().sorted().cloned().collect();

  // Build DiscreteStates
  let discrete_states = DiscreteStates::from_values(unique_values.iter().map(String::as_str), missing_data);
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

  // Compute equilibrium frequencies (pi)
  let pi = if let Some(weights_filepath) = weights {
    let weights_map = read_discrete_attrs::<f64>(
      weights_filepath,
      &Some(attribute.clone()),
      &Some("weight".to_owned()),
      |s| Ok(s.parse::<f64>()?),
    )?
    .0;

    let weights_keys: IndexSet<String> = weights_map.keys().sorted().cloned().collect();

    let coverage = validate_weight_coverage(&unique_values, &weights_keys, missing_data, *missing_weights_threshold)?;

    if !coverage.missing_values.is_empty() {
      warn!(
        "Mugration: discrete attributes missing from weights file: {} (ratio: {:.3})",
        coverage.missing_values.iter().join(", "),
        coverage.missing_ratio
      );
    }

    compute_pi_from_weights(&discrete_states, &weights_map)
  } else {
    compute_pi_uniform(n_states)
  };

  // Add pseudo-counts if specified
  let pi = apply_pseudo_counts(pi, *pc);

  // Create GTR model
  let gtr = GTR::new(GTRParams {
    n_states,
    mu: 1.0,
    W: None, // uniform rates
    pi,
  })?;

  // Create partition and attach traits
  let mut partition = PartitionDiscrete::new(0, gtr, discrete_states);

  // Convert attr_values to BTreeMap<String, String>
  let traits: BTreeMap<String, String> = attr_values.into_iter().collect();
  attach_traits(&mut partition, &graph, &traits)?;

  // Run discrete marginal reconstruction
  let log_lh = run_discrete_marginal(&graph, &mut partition)?;
  info!("Mugration: total log likelihood = {log_lh:.4}");

  // Write output files
  write_annotated_tree(&graph, &partition, attribute, outdir)?;
  write_gtr_json(&partition, attribute, outdir)?;

  if let Some(confidence_path) = confidence {
    write_confidence_csv(&graph, &partition, confidence_path)?;
  }

  info!("Mugration: wrote output to {}", outdir.display());
  Ok(())
}

fn write_annotated_tree(
  graph: &GraphAncestral,
  partition: &PartitionDiscrete,
  attribute: &str,
  outdir: &std::path::Path,
) -> Result<(), Report> {
  let leaf_names = graph
    .get_leaves()
    .iter()
    .filter_map(|node| {
      let node = node.read_arc();
      let payload = node.payload().read_arc();
      payload.name().map(|name| name.as_ref().to_owned())
    })
    .join(" ");
  let nwk = write_annotated_tree_nwk(graph, partition, attribute)?;
  let nexus = format!(
    r#"#NEXUS
Begin Taxa;
  Dimensions NTax={};
  TaxLabels {};
End;
Begin Trees;
  Tree tree1={nwk};
End;
"#,
    graph.num_leaves(),
    leaf_names
  );
  fs::write(outdir.join("annotated_tree.nexus"), format!("{nexus}\n"))?;

  // Write trait assignments as separate CSV
  let mut file = fs::File::create(outdir.join("traits.csv"))?;
  writeln!(file, "node,{attribute}")?;

  for node in graph.get_nodes() {
    let node_guard = node.read_arc();
    let node_key = node_guard.key();
    let payload = node_guard.payload().read_arc();
    let node_name = payload
      .name()
      .map_or_else(|| format!("node_{}", node_key.0), |n| n.as_ref().to_owned());

    if let Some(trait_value) = partition.get_reconstructed_trait(node_key) {
      writeln!(file, "{node_name},{trait_value}")?;
    }
  }

  Ok(())
}

fn write_annotated_tree_nwk(
  graph: &GraphAncestral,
  partition: &PartitionDiscrete,
  attribute: &str,
) -> Result<String, Report> {
  let roots = graph.get_roots();
  let root = if roots.len() == 1 {
    &roots[0]
  } else if roots.is_empty() {
    return make_error!("When converting graph to mugration-annotated Newick format: No roots found.");
  } else {
    return make_error!(
      "When converting graph to mugration-annotated Newick format: Multiple roots are not supported."
    );
  };

  let mut nwk = String::new();
  let mut stack = vec![(Arc::clone(root), None, 0_usize)];
  while let Some((node, edge, child_visit)) = stack.pop() {
    let children = graph.children_of(&node.read_arc()).into_iter().collect_vec();

    if child_visit < children.len() {
      stack.push((node, edge, child_visit + 1));

      if child_visit == 0 {
        write!(&mut nwk, "(")?;
      } else {
        write!(&mut nwk, ",")?;
      }

      let (child, child_edge) = &children[child_visit];
      stack.push((Arc::clone(child), Some(Arc::clone(child_edge)), 0));
      continue;
    }

    if child_visit > 0 {
      write!(&mut nwk, ")")?;
    }

    let node_key = node.read_arc().key();
    let (name, mut comments) = {
      let node_payload = node.read_arc().payload().read_arc();
      (
        node_payload.nwk_name().map(|node_name| node_name.as_ref().to_owned()),
        node_payload.nwk_comments(),
      )
    };
    if let Some(trait_value) = partition.get_reconstructed_trait(node_key) {
      comments.insert(attribute.to_owned(), trait_value);
    }

    let weight = edge.and_then(|edge| edge.read_arc().payload().read_arc().nwk_weight());

    if let Some(name) = name {
      write!(&mut nwk, "{name}")?;
    }

    if let Some(weight) = weight {
      write!(&mut nwk, ":{}", format_weight(weight, &NwkWriteOptions::default()))?;
    }

    if !comments.is_empty() {
      let comments = comments
        .iter()
        .filter(|(_, value)| !value.is_empty())
        .map(|(key, value)| format!("[&{key}=\"{value}\"]"))
        .join("");
      if !comments.is_empty() {
        write!(&mut nwk, "{comments}")?;
      }
    }
  }

  write!(&mut nwk, ";")?;
  Ok(nwk)
}

#[derive(Serialize)]
struct GTROutput {
  attribute: String,
  n_states: usize,
  states: Vec<String>,
  pi: Vec<f64>,
  mu: f64,
}

fn write_gtr_json(partition: &PartitionDiscrete, attribute: &str, outdir: &std::path::Path) -> Result<(), Report> {
  let output = GTROutput {
    attribute: attribute.to_owned(),
    n_states: partition.n_states(),
    states: partition.states.iter().map(|s| s.to_owned()).collect(),
    pi: partition.gtr.pi.to_vec(),
    mu: partition.gtr.mu,
  };

  let json = serde_json::to_string_pretty(&output)?;
  let mut file = fs::File::create(outdir.join("gtr.json"))?;
  file.write_all(json.as_bytes())?;

  Ok(())
}

fn write_confidence_csv(
  graph: &GraphAncestral,
  partition: &PartitionDiscrete,
  output_path: &std::path::Path,
) -> Result<(), Report> {
  let mut file = fs::File::create(output_path)?;

  // Header
  let state_names: Vec<&str> = partition.states.iter().collect();
  writeln!(file, "node,{}", state_names.join(","))?;

  // Data rows
  for node in graph.get_nodes() {
    let node_guard = node.read_arc();
    let node_key = node_guard.key();
    let payload = node_guard.payload().read_arc();
    let node_name = payload
      .name()
      .map_or_else(|| format!("node_{}", node_key.0), |n| n.as_ref().to_owned());

    if let Some(confidence) = partition.get_confidence(node_key) {
      let probs: Vec<String> = confidence.iter().map(|p| format!("{p:.6}")).collect();
      writeln!(file, "{node_name},{}", probs.join(","))?;
    }
  }

  Ok(())
}
