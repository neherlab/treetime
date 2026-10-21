use crate::annotated_graph::{AnnotatedTreeView, TreeTraits};
use crate::auspice::ensure_finite;
use crate::trait_profile::{build_confidence_map, trait_change};
use deser_value::{Value, to_value};
use eyre::Report;
use itertools::{Itertools, izip};
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
  let branches = trait_branches(tree, traits)?;
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
    rate: model.mu(),
    alphabet,
    equilibrium_probabilities: model.pi().to_vec(),
    transition_matrix: model.W().outer_iter().map(|row| row.to_vec()).collect(),
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
    let posterior = &traits.posteriors[&key];
    let mut fields = BTreeMap::new();
    fields.insert(traits.attribute.to_owned(), Value::from(traits.state_name(key)?));
    for (state, probability) in izip!(traits.states.iter(), posterior.profile()) {
      ensure_finite(*probability, &name, &format!("trait state '{state}' probability"))?;
    }
    let entropy = posterior.entropy();
    ensure_finite(entropy, &name, "trait entropy")?;
    let confidence = build_confidence_map(traits.states, posterior.profile());
    if !confidence.is_empty() {
      fields.insert(confidence_key.clone(), to_value(&confidence)?);
    }
    fields.insert(entropy_key.clone(), Value::from(entropy));
    nodes.insert(name, AugurNodeDataJsonTraitsNode { fields });
  }
  Ok(nodes)
}

fn trait_branches(
  tree: &AnnotatedTreeView<'_>,
  traits: &TreeTraits<'_>,
) -> Result<BTreeMap<String, AugurNodeDataJsonTraitsBranches>, Report> {
  let graph = tree.graph();
  graph
    .graph
    .get_nodes()
    .map(|node| {
      let key = node.key();
      let branch = branch_label(tree, traits, key)?.map(|label| {
        let branch = AugurNodeDataJsonTraitsBranches {
          labels: Some(btreemap! { traits.attribute.to_owned() => label }),
          other: BTreeMap::new(),
        };
        (node_name_or_key(key, graph.names[&key].as_deref()), branch)
      });
      Ok(branch)
    })
    .flatten_ok()
    .collect()
}

fn branch_label(
  tree: &AnnotatedTreeView<'_>,
  traits: &TreeTraits<'_>,
  key: GraphNodeKey,
) -> Result<Option<String>, Report> {
  match tree.tree().parent(key) {
    Some((parent_key, _)) => trait_change(traits, parent_key, key),
    None => Ok(Some(traits.state_name(key)?.to_owned())),
  }
}
