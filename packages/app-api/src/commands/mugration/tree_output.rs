use crate::commands::mugration::augur_node_data::{build_confidence_map, compute_entropy};
use crate::commands::shared::tree_output::{
  TraitValue, auspice_data, auspice_from_graph, auspice_node, coloring, cumulative_branch_length_from, ensure_finite,
  finite_number, generation_date, mutation_free_mat, node_name_value, write_tree_outputs,
};
use eyre::Report;
use maplit::btreemap;
use serde_json::json;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime::mugration::result::{MugrationNodeOut, MugrationOutputMaps};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::auspice_types::{AuspiceTree, AuspiceTreeBranchAttrsLabels};
use treetime_io::graph::TreeWriteKind;
use treetime_io::nwk::CommentProviders;
use treetime_io::usher_mat::UsherTree;

pub fn write_mugration_tree_outputs(
  graph: &Graph,
  nodes: &BTreeMap<GraphNodeKey, MugrationNodeOut>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  maps: &MugrationOutputMaps,
  attribute: &str,
  outputs: &BTreeMap<TreeWriteKind, PathBuf>,
  providers: &CommentProviders,
) -> Result<(), Report> {
  let updated = generation_date();
  let names: BTreeMap<GraphNodeKey, Option<String>> =
    nodes.iter().map(|(key, node)| (*key, node.name.clone())).collect();
  write_tree_outputs(
    graph,
    &names,
    branch_lengths,
    branch_lengths,
    outputs,
    providers,
    "mugration",
    || mugration_to_auspice(graph, nodes, branch_lengths, maps, attribute, &updated),
    || mugration_to_mat(graph, &names, branch_lengths),
  )
}

pub(crate) fn mugration_to_auspice(
  graph: &Graph,
  nodes: &BTreeMap<GraphNodeKey, MugrationNodeOut>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  maps: &MugrationOutputMaps,
  attribute: &str,
  updated: &str,
) -> Result<AuspiceTree, Report> {
  let data = auspice_data(
    "TreeTime mugration analysis",
    updated,
    vec![coloring(attribute, attribute, "categorical")],
    vec![attribute.to_owned()],
    Some(attribute.to_owned()),
    None,
    None,
    false,
  );
  auspice_from_graph(graph, data, |context| {
    let name = node_name_value(context.node_key, nodes[&context.node_key].name.as_deref());
    let traits = mugration_traits(graph, maps, context.node_key, &name, attribute)?;
    Ok(auspice_node(
      name.clone(),
      finite_number(
        cumulative_branch_length_from(graph, branch_lengths, context.node_key)?,
        6,
        "mugration",
        &name,
        "div",
      )?,
      None,
      None,
      None,
      traits,
      BTreeMap::new(),
      mugration_transition_label(graph, maps, context.node_key, attribute)?,
    ))
  })
}

pub(crate) fn mugration_to_mat(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  nwk_weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<UsherTree, Report> {
  mutation_free_mat(graph, names, nwk_weights)
}

fn mugration_traits(
  graph: &Graph,
  maps: &MugrationOutputMaps,
  node_key: GraphNodeKey,
  node_name: &str,
  attribute: &str,
) -> Result<BTreeMap<String, TraitValue>, Report> {
  let Some(value) = maps.reconstructed_traits[&node_key].clone() else {
    return Ok(BTreeMap::new());
  };
  let profile = maps.confidences[&node_key].clone();
  if let Some(profile) = profile.as_ref() {
    for (state, probability) in maps.states.iter().zip(profile) {
      ensure_finite(
        *probability,
        "mugration",
        node_name,
        &format!("trait state '{state}' probability"),
      )?;
    }
  }
  let confidence = profile
    .as_ref()
    .map(|profile| build_confidence_map(&maps.states, profile))
    .unwrap_or_default();
  let entropy = profile.as_ref().map(compute_entropy);
  if let Some(entropy) = entropy {
    ensure_finite(entropy, "mugration", node_name, "trait entropy")?;
  }
  Ok(btreemap! {
    attribute.to_owned() => TraitValue { value, confidence, entropy },
  })
}

fn mugration_transition(
  graph: &Graph,
  maps: &MugrationOutputMaps,
  node_key: GraphNodeKey,
) -> Result<Option<(String, String)>, Report> {
  let Some((parent_key, _edge_key)) = graph.node_parent(node_key)? else {
    return Ok(None);
  };
  let parent = maps.reconstructed_traits[&parent_key].clone();
  let child = maps.reconstructed_traits[&node_key].clone();
  Ok(match (parent, child) {
    (Some(parent), Some(child)) if parent != child => Some((parent, child)),
    _ => None,
  })
}

fn mugration_transition_label(
  graph: &Graph,
  maps: &MugrationOutputMaps,
  node_key: GraphNodeKey,
  attribute: &str,
) -> Result<Option<AuspiceTreeBranchAttrsLabels>, Report> {
  let Some((parent, child)) = mugration_transition(graph, maps, node_key)? else {
    return Ok(None);
  };
  Ok(Some(AuspiceTreeBranchAttrsLabels {
    aa: None,
    clade: None,
    other: json!({ attribute.to_owned(): format!("{parent} → {child}") }),
  }))
}
