use crate::annotated_graph::{AnnotatedTreeView, TreeTraits};
use crate::auspice::ensure_finite;
use crate::trait_profile::{build_confidence_map, compute_entropy};
use deser_value::{Value, to_value};
use eyre::Report;
use itertools::izip;
use maplit::btreemap;
use std::collections::BTreeMap;
use std::path::Path;
use treetime::gtr::gtr::GTR;
use treetime_graph::assign_node_names::node_name_or_key;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::io::json::{JsonPretty, json_write_file};
use treetime_utils::make_internal_report;
use util_augur_node_data_json::{
  AugurNodeDataJsonGeneratedBy, AugurNodeDataJsonTraitModel, AugurNodeDataJsonTraits, AugurNodeDataJsonTraitsBranches,
  AugurNodeDataJsonTraitsMeta, AugurNodeDataJsonTraitsNode,
};

pub fn write_augur_node_data_traits(tree: &AnnotatedTreeView<'_>, model: &GTR, path: &Path) -> Result<(), Report> {
  json_write_file(path, &build_augur_node_data_traits(tree, model)?, JsonPretty(true))
}

pub fn build_augur_node_data_traits(
  tree: &AnnotatedTreeView<'_>,
  model: &GTR,
) -> Result<AugurNodeDataJsonTraits, Report> {
  let traits = tree
    .graph()
    .traits
    .as_ref()
    .ok_or_else(|| make_internal_report!("Trait node data requires reconstructed traits"))?;
  let branches = trait_branches(tree, traits);
  let mut other = BTreeMap::new();
  if !branches.is_empty() {
    other.insert("branches".to_owned(), to_value(&branches)?);
  }
  Ok(AugurNodeDataJsonTraits {
    generated_by: Some(AugurNodeDataJsonGeneratedBy {
      program: "treetime".to_owned(),
      version: env!("CARGO_PKG_VERSION").to_owned(),
    }),
    metadata: AugurNodeDataJsonTraitsMeta {
      models: Some(trait_models(traits, model)),
      other,
    },
    nodes: trait_nodes(tree, traits)?,
  })
}

fn trait_models(traits: &TreeTraits<'_>, model: &GTR) -> BTreeMap<String, AugurNodeDataJsonTraitModel> {
  let alphabet = traits.states.iter().chain(["?"]).map(str::to_owned).collect();
  let model = AugurNodeDataJsonTraitModel {
    rate: model.mu,
    alphabet,
    equilibrium_probabilities: model.pi.to_vec(),
    transition_matrix: model.W.outer_iter().map(|row| row.to_vec()).collect(),
    other: BTreeMap::new(),
  };
  btreemap! { traits.attribute.to_owned() => model }
}

fn trait_nodes(
  tree: &AnnotatedTreeView<'_>,
  traits: &TreeTraits<'_>,
) -> Result<BTreeMap<String, AugurNodeDataJsonTraitsNode>, Report> {
  let graph = tree.graph();
  let confidence_key = format!("{}_confidence", traits.attribute);
  let entropy_key = format!("{}_entropy", traits.attribute);
  let mut nodes = BTreeMap::new();
  for node in graph.graph.get_nodes() {
    let key = node.key();
    let name = node_name_or_key(key, graph.names[&key].as_deref());
    let mut fields = BTreeMap::new();
    if let Some(value) = &traits.values[&key] {
      fields.insert(traits.attribute.to_owned(), Value::from(value.as_str()));
    }
    if let Some(profile) = &traits.profiles[&key] {
      for (state, probability) in izip!(traits.states.iter(), profile) {
        ensure_finite(*probability, &name, &format!("trait state '{state}' probability"))?;
      }
      let entropy = compute_entropy(profile);
      ensure_finite(entropy, &name, "trait entropy")?;
      let confidence = build_confidence_map(traits.states, profile);
      if !confidence.is_empty() {
        fields.insert(confidence_key.clone(), to_value(&confidence)?);
      }
      fields.insert(entropy_key.clone(), Value::from(entropy));
    }
    nodes.insert(name, AugurNodeDataJsonTraitsNode { fields });
  }
  Ok(nodes)
}

fn trait_branches(
  tree: &AnnotatedTreeView<'_>,
  traits: &TreeTraits<'_>,
) -> BTreeMap<String, AugurNodeDataJsonTraitsBranches> {
  let graph = tree.graph();
  graph
    .graph
    .get_nodes()
    .filter_map(|node| {
      let key = node.key();
      let label = branch_label(tree, traits, key)?;
      let branch = AugurNodeDataJsonTraitsBranches {
        labels: Some(btreemap! { traits.attribute.to_owned() => label }),
        other: BTreeMap::new(),
      };
      Some((node_name_or_key(key, graph.names[&key].as_deref()), branch))
    })
    .collect()
}

fn branch_label(tree: &AnnotatedTreeView<'_>, traits: &TreeTraits<'_>, key: GraphNodeKey) -> Option<String> {
  let child = traits.values[&key].as_ref()?;
  let Some((parent_key, _)) = tree.tree().parent(key) else {
    return Some(child.clone());
  };
  let parent = traits.values[&parent_key].as_ref()?;
  (parent != child).then(|| format!("{parent} \u{2192} {child}"))
}
