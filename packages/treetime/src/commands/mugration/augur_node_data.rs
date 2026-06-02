use crate::mugration::result::MugrationResult;
use crate::partition::discrete_states::DiscreteStates;
use crate::partition::marginal_discrete::PartitionMarginalDiscrete;
use crate::partition::traits::HasGtr;
use crate::payload::ancestral::GraphAncestral;
use eyre::Report;
use std::collections::BTreeMap;
use std::path::Path;
use treetime_graph::node::{GraphNodeKey, Named};
use treetime_utils::io::json::{JsonPretty, json_write_file};
use util_augur_node_data_json::{
  AugurNodeDataJsonGeneratedBy, AugurNodeDataJsonTraitModel, AugurNodeDataJsonTraits, AugurNodeDataJsonTraitsBranches,
  AugurNodeDataJsonTraitsMeta, AugurNodeDataJsonTraitsNode,
};

pub fn build_augur_node_data_json(result: &MugrationResult) -> Result<AugurNodeDataJsonTraits, Report> {
  let attribute = &result.traits.attribute;
  let partition = &result.partition;
  let graph = &result.graph;

  let models = build_models(attribute, partition);
  let nodes = build_nodes(attribute, graph, partition);
  let branches = build_branches(attribute, graph, partition);

  Ok(AugurNodeDataJsonTraits {
    generated_by: Some(AugurNodeDataJsonGeneratedBy {
      program: "treetime".to_owned(),
      version: env!("CARGO_PKG_VERSION").to_owned(),
    }),
    metadata: AugurNodeDataJsonTraitsMeta {
      models: Some(models),
      other: if branches.is_empty() {
        BTreeMap::new()
      } else {
        let mut other = BTreeMap::new();
        other.insert("branches".to_owned(), serde_json::to_value(branches)?);
        other
      },
    },
    nodes,
  })
}

pub fn write_augur_node_data_json(result: &MugrationResult, path: &Path) -> Result<(), Report> {
  let data = build_augur_node_data_json(result)?;
  json_write_file(path, &data, JsonPretty(true))?;
  Ok(())
}

fn build_models(
  attribute: &str,
  partition: &PartitionMarginalDiscrete,
) -> BTreeMap<String, AugurNodeDataJsonTraitModel> {
  let gtr = partition.gtr();
  let n_states = partition.n_states();

  // Alphabet includes missing data marker "?" at the end (n_states+1 elements)
  let mut alphabet: Vec<String> = partition.states.iter().map(|s| s.to_owned()).collect();
  alphabet.push("?".to_owned());

  // Equilibrium probabilities exclude missing (n_states elements)
  let equilibrium_probabilities: Vec<f64> = (0..n_states).map(|i| gtr.pi[i]).collect();

  // Transition matrix excludes missing (n_states x n_states)
  let transition_matrix: Vec<Vec<f64>> = (0..n_states)
    .map(|i| (0..n_states).map(|j| gtr.W[[i, j]]).collect())
    .collect();

  let mut models = BTreeMap::new();
  models.insert(
    attribute.to_owned(),
    AugurNodeDataJsonTraitModel {
      rate: gtr.mu,
      alphabet,
      equilibrium_probabilities,
      transition_matrix,
      other: BTreeMap::new(),
    },
  );
  models
}

fn build_nodes(
  attribute: &str,
  graph: &GraphAncestral,
  partition: &PartitionMarginalDiscrete,
) -> BTreeMap<String, AugurNodeDataJsonTraitsNode> {
  let confidence_key = format!("{attribute}_confidence");
  let entropy_key = format!("{attribute}_entropy");

  let mut nodes = BTreeMap::new();

  for node in graph.get_nodes() {
    let node_guard = node.read_arc();
    let node_key = node_guard.key();
    let payload = node_guard.payload().read_arc();
    let node_name = payload
      .name()
      .map_or_else(|| format!("node_{}", node_key.0), |n| n.as_ref().to_owned());

    let mut fields = BTreeMap::new();

    if let Some(trait_value) = partition.get_reconstructed_trait(node_key) {
      fields.insert(attribute.to_owned(), serde_json::Value::String(trait_value));
    }

    if let Some(profile) = partition.get_confidence(node_key) {
      let confidence = build_confidence_map(&partition.states, &profile);
      if !confidence.is_empty() {
        fields.insert(confidence_key.clone(), serde_json::to_value(&confidence).unwrap());
      }

      let entropy = compute_entropy(&profile);
      fields.insert(entropy_key.clone(), serde_json::json!(entropy));
    }

    nodes.insert(node_name, AugurNodeDataJsonTraitsNode { fields });
  }

  nodes
}

/// Build confidence map: state -> probability, sorted descending, filtered > 0.001.
fn build_confidence_map(states: &DiscreteStates, profile: &ndarray::Array1<f64>) -> BTreeMap<String, f64> {
  let mut pairs: Vec<(&str, f64)> = states.iter().zip(profile.iter()).map(|(s, &p)| (s, p)).collect();
  pairs.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
  pairs
    .into_iter()
    .filter(|(_, p)| *p > 0.001)
    .map(|(s, p)| (s.to_owned(), p))
    .collect()
}

/// Shannon entropy: -sum(p * ln(p + 1e-12)) over real states (excludes missing).
fn compute_entropy(profile: &ndarray::Array1<f64>) -> f64 {
  const TINY: f64 = 1e-12;
  -profile.iter().map(|&p| p * (p + TINY).ln()).sum::<f64>()
}

fn build_branches(
  attribute: &str,
  graph: &GraphAncestral,
  partition: &PartitionMarginalDiscrete,
) -> BTreeMap<String, AugurNodeDataJsonTraitsBranches> {
  let root_key = graph.get_exactly_one_root().ok().map(|r| r.read_arc().key());
  let mut branches = BTreeMap::new();

  let parent_traits = build_parent_trait_map(graph, partition);

  for node in graph.get_nodes() {
    let node_guard = node.read_arc();
    let node_key = node_guard.key();
    let payload = node_guard.payload().read_arc();
    let node_name = payload
      .name()
      .map_or_else(|| format!("node_{}", node_key.0), |n| n.as_ref().to_owned());

    let child_trait = partition.get_reconstructed_trait(node_key);

    let label = if Some(node_key) == root_key {
      // Root gets just the state name (no arrow)
      child_trait.clone()
    } else {
      let parent_trait = parent_traits.get(&node_key).and_then(|t| t.as_deref());
      match (parent_trait, child_trait.as_deref()) {
        (Some(parent), Some(child)) if parent != child => Some(format!("{parent} \u{2192} {child}")),
        _ => None,
      }
    };

    if let Some(label) = label {
      let mut labels = BTreeMap::new();
      labels.insert(attribute.to_owned(), label);
      branches.insert(
        node_name,
        AugurNodeDataJsonTraitsBranches {
          labels: Some(labels),
          other: BTreeMap::new(),
        },
      );
    }
  }

  branches
}

fn build_parent_trait_map(
  graph: &GraphAncestral,
  partition: &PartitionMarginalDiscrete,
) -> BTreeMap<GraphNodeKey, Option<String>> {
  let mut map = BTreeMap::new();
  for node in graph.get_nodes() {
    let node_guard = node.read_arc();
    let node_key = node_guard.key();
    let inbound = node_guard.inbound().to_vec();
    if let Some(parent_edge_key) = inbound.first() {
      let parent_node_key = graph.get_source_node_key(*parent_edge_key).ok();
      let parent_trait = parent_node_key.and_then(|k| partition.get_reconstructed_trait(k));
      map.insert(node_key, parent_trait);
    }
  }
  map
}
