use crate::mugration_result::MugrationNodeOut;
use crate::tree_output::{
  TraitValue, auspice_data, auspice_from_graph, auspice_node, coloring, cumulative_branch_length_from, ensure_finite,
  finite_number, generation_date, mutation_free_mat, node_name_value, write_tree_outputs,
};
use eyre::Report;
use maplit::btreemap;
use serde_json::json;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime::mugration::pipeline::MugrationOutput;
use treetime::partition::storage::discrete::DiscreteStates;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::auspice_types::{AuspiceTree, AuspiceTreeBranchAttrsLabels};
use treetime_io::graph::TreeWriteKind;
use treetime_io::nwk::CommentProviders;
use treetime_io::usher_mat::UsherTree;

pub fn write_mugration_tree_outputs(
  output: &MugrationOutput,
  nodes: &BTreeMap<GraphNodeKey, MugrationNodeOut>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  attribute: &str,
  outputs: &BTreeMap<TreeWriteKind, PathBuf>,
  providers: &CommentProviders,
) -> Result<(), Report> {
  let graph = &output.graph;
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
    || mugration_to_auspice(graph, nodes, branch_lengths, output, attribute, &updated),
    || mugration_to_mat(graph, &names, branch_lengths),
  )
}

pub(crate) fn mugration_to_auspice(
  graph: &Graph,
  nodes: &BTreeMap<GraphNodeKey, MugrationNodeOut>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  output: &MugrationOutput,
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
    let traits = mugration_traits(graph, output, context.node_key, &name, attribute)?;
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
      mugration_transition_label(graph, output, context.node_key, attribute)?,
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
  output: &MugrationOutput,
  node_key: GraphNodeKey,
  node_name: &str,
  attribute: &str,
) -> Result<BTreeMap<String, TraitValue>, Report> {
  let Some(value) = output.reconstructed_traits[&node_key].clone() else {
    return Ok(BTreeMap::new());
  };
  let profile = output.confidences[&node_key].clone();
  if let Some(profile) = profile.as_ref() {
    for (state, probability) in output.states.iter().zip(profile) {
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
    .map(|profile| build_confidence_map(&output.states, profile))
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
  output: &MugrationOutput,
  node_key: GraphNodeKey,
) -> Result<Option<(String, String)>, Report> {
  let Some((parent_key, _edge_key)) = graph.node_parent(node_key)? else {
    return Ok(None);
  };
  let parent = output.reconstructed_traits[&parent_key].clone();
  let child = output.reconstructed_traits[&node_key].clone();
  Ok(match (parent, child) {
    (Some(parent), Some(child)) if parent != child => Some((parent, child)),
    _ => None,
  })
}

fn mugration_transition_label(
  graph: &Graph,
  output: &MugrationOutput,
  node_key: GraphNodeKey,
  attribute: &str,
) -> Result<Option<AuspiceTreeBranchAttrsLabels>, Report> {
  let Some((parent, child)) = mugration_transition(graph, output, node_key)? else {
    return Ok(None);
  };
  Ok(Some(AuspiceTreeBranchAttrsLabels {
    aa: None,
    clade: None,
    other: json!({ attribute.to_owned(): format!("{parent} → {child}") }),
  }))
}

pub(crate) fn build_confidence_map(states: &DiscreteStates, profile: &ndarray::Array1<f64>) -> BTreeMap<String, f64> {
  let mut pairs: Vec<(&str, f64)> = states.iter().zip(profile.iter()).map(|(s, &p)| (s, p)).collect();
  pairs.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
  pairs
    .into_iter()
    .filter(|(_, p)| *p > 0.001)
    .map(|(s, p)| (s.to_owned(), p))
    .collect()
}

pub(crate) fn compute_entropy(profile: &ndarray::Array1<f64>) -> f64 {
  const TINY: f64 = 1e-12;
  -profile.iter().map(|&p| p * (p + TINY).ln()).sum::<f64>()
}
