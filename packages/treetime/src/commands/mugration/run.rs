use crate::commands::mugration::args::TreetimeMugrationArgs;
use crate::commands::mugration::discrete_marginal::{attach_traits, run_discrete_marginal};
use crate::commands::mugration::input::MugrationInput;
use crate::commands::mugration::output::{MugrationGtrOutput, MugrationResult, MugrationTraitsOutput};
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
use statrs::statistics::Statistics;
use std::collections::BTreeMap;
use std::fmt::Write;
use std::fs;
use std::sync::Arc;
use treetime_graph::node::Named;
use treetime_io::discrete_states_csv::read_discrete_attrs;
use treetime_io::json::{JsonPretty, json_write_file};
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

/// Parse CLI arguments into a MugrationInput struct.
///
/// Reads tree, trait values, and optional weights from files and assembles
/// them into an in-memory input suitable for the execution layer.
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

  // Read tree
  let tree_path = tree.as_ref().ok_or_else(|| eyre::eyre!("Tree file is required"))?;
  let graph: GraphAncestral = nwk_read_file(tree_path)?;

  // Read trait values
  let (attr_values, _attr_name) =
    read_discrete_attrs::<String>(states, name_column, &Some(attribute.clone()), |s| Ok(s.to_owned()))?;
  let traits: BTreeMap<String, String> = attr_values.into_iter().collect();

  // Read weights if provided
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
  })
}

/// Execute mugration from parsed input and return structured results.
///
/// This is the reusable execution core that can be called by tests directly
/// without file I/O. The command wrapper handles parsing and output writing.
pub fn execute_mugration(input: MugrationInput) -> Result<MugrationResult, Report> {
  let MugrationInput {
    graph,
    traits,
    attribute,
    weights,
    missing_data,
    pc,
    missing_weights_threshold,
  } = input;

  let unique_values: IndexSet<String> = traits.values().sorted().cloned().collect();

  // Build DiscreteStates
  let discrete_states = DiscreteStates::from_values(unique_values.iter().map(String::as_str), &missing_data);
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
  let pi = match &weights {
    Some(weights_map) => {
      let weights_keys: IndexSet<String> = weights_map.keys().sorted().cloned().collect();

      let coverage = validate_weight_coverage(&unique_values, &weights_keys, &missing_data, missing_weights_threshold)?;

      if !coverage.missing_values.is_empty() {
        warn!(
          "Mugration: discrete attributes missing from weights file: {} (ratio: {:.3})",
          coverage.missing_values.iter().join(", "),
          coverage.missing_ratio
        );
      }

      compute_pi_from_weights(&discrete_states, weights_map)
    },
    None => compute_pi_uniform(n_states),
  };

  // Add pseudo-counts if specified
  let pi = apply_pseudo_counts(pi, pc);

  // Create GTR model
  let gtr = GTR::new(GTRParams {
    n_states,
    mu: 1.0,
    W: None, // uniform rates
    pi,
  })?;

  // Create partition and attach traits
  let mut partition = PartitionDiscrete::new(0, gtr, discrete_states);
  attach_traits(&mut partition, &graph, &traits)?;

  // Run discrete marginal reconstruction
  let log_lh = run_discrete_marginal(&graph, &mut partition)?;
  info!("Mugration: total log likelihood = {log_lh:.4}");

  Ok(MugrationResult::new(graph, partition, &attribute, log_lh))
}

pub fn run_mugration(mugration_args: &TreetimeMugrationArgs) -> Result<(), Report> {
  let outdir = &mugration_args.outdir;
  fs::create_dir_all(outdir)?;

  // Parse input and execute
  let input = parse_mugration_input(mugration_args)?;
  let result = execute_mugration(input)?;

  // Write output files
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
  partition: &PartitionDiscrete,
  traits: &MugrationTraitsOutput,
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
  let nwk = write_annotated_tree_nwk(graph, partition, &traits.attribute)?;
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
  fs::write(outdir.join("traits.csv"), traits.render_csv())?;

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

fn write_gtr_json_file(gtr: &MugrationGtrOutput, outdir: &std::path::Path) -> Result<(), Report> {
  json_write_file(outdir.join("gtr.json"), gtr, JsonPretty(true))
}

fn write_confidence_csv(result: &MugrationResult, output_path: &std::path::Path) -> Result<(), Report> {
  fs::write(output_path, result.confidence.render_csv())?;
  Ok(())
}
