use crate::mugration_result::MugrationResult;
use crate::mugration_tree_output::{build_confidence_map, compute_entropy};
use eyre::Report;
use std::collections::BTreeMap;
use std::path::Path;
use treetime::mugration::pipeline::MugrationOutput;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::io::json::{JsonPretty, json_write_file};
use util_augur_node_data_json::{
  AugurNodeDataJsonGeneratedBy, AugurNodeDataJsonTraitModel, AugurNodeDataJsonTraits, AugurNodeDataJsonTraitsBranches,
  AugurNodeDataJsonTraitsMeta, AugurNodeDataJsonTraitsNode,
};

pub fn build_augur_node_data_json(
  result: &MugrationResult,
  output: &MugrationOutput,
) -> Result<AugurNodeDataJsonTraits, Report> {
  let attribute = &result.traits.attribute;
  let graph = &output.graph;
  let names: BTreeMap<GraphNodeKey, Option<String>> = result
    .nodes
    .iter()
    .map(|(key, node)| (*key, node.name.clone()))
    .collect();

  let models = build_models(attribute, output);
  let nodes = build_nodes(attribute, graph, &names, output);
  let branches = build_branches(attribute, graph, &names, output);

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

pub fn write_augur_node_data_json(
  result: &MugrationResult,
  output: &MugrationOutput,
  path: &Path,
) -> Result<(), Report> {
  let data = build_augur_node_data_json(result, output)?;
  json_write_file(path, &data, JsonPretty(true))?;
  Ok(())
}

fn build_models(attribute: &str, output: &MugrationOutput) -> BTreeMap<String, AugurNodeDataJsonTraitModel> {
  let gtr = &output.gtr;
  let n_states = output.n_states;

  // Alphabet includes missing data marker "?" at the end (n_states+1 elements)
  let mut alphabet: Vec<String> = output.states.iter().map(|s| s.to_owned()).collect();
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
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  output: &MugrationOutput,
) -> BTreeMap<String, AugurNodeDataJsonTraitsNode> {
  let confidence_key = format!("{attribute}_confidence");
  let entropy_key = format!("{attribute}_entropy");

  let mut nodes = BTreeMap::new();

  for node in graph.get_nodes() {
    let node_guard = node;
    let node_key = node_guard.key();
    let node_name = names[&node_key]
      .as_deref()
      .map_or_else(|| format!("node_{}", node_key.0), str::to_owned);

    let mut fields = BTreeMap::new();

    if let Some(trait_value) = output.reconstructed_traits[&node_key].clone() {
      fields.insert(attribute.to_owned(), serde_json::Value::String(trait_value));
    }

    if let Some(profile) = output.confidences[&node_key].as_ref() {
      let confidence = build_confidence_map(&output.states, profile);
      if !confidence.is_empty() {
        fields.insert(confidence_key.clone(), serde_json::to_value(&confidence).unwrap());
      }

      let entropy = compute_entropy(profile);
      fields.insert(entropy_key.clone(), serde_json::json!(entropy));
    }

    nodes.insert(node_name, AugurNodeDataJsonTraitsNode { fields });
  }

  nodes
}

fn build_branches(
  attribute: &str,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  output: &MugrationOutput,
) -> BTreeMap<String, AugurNodeDataJsonTraitsBranches> {
  let root_key = graph.get_exactly_one_root().ok().map(|r| r.key());
  let mut branches = BTreeMap::new();

  let parent_traits = build_parent_trait_map(graph, output);

  for node in graph.get_nodes() {
    let node_guard = node;
    let node_key = node_guard.key();
    let node_name = names[&node_key]
      .as_deref()
      .map_or_else(|| format!("node_{}", node_key.0), str::to_owned);

    let child_trait = output.reconstructed_traits[&node_key].clone();

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

fn build_parent_trait_map(graph: &Graph, output: &MugrationOutput) -> BTreeMap<GraphNodeKey, Option<String>> {
  let mut map = BTreeMap::new();
  for node in graph.get_nodes() {
    let node_guard = node;
    let node_key = node_guard.key();
    let inbound = node_guard.inbound().to_vec();
    if let Some(parent_edge_key) = inbound.first() {
      let parent_node_key = graph.get_source_node_key(*parent_edge_key).ok();
      let parent_trait = parent_node_key.and_then(|k| output.reconstructed_traits[&k].clone());
      map.insert(node_key, parent_trait);
    }
  }
  map
}
