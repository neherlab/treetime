use chrono::Utc;
use eyre::{Report, WrapErr};
use maplit::{btreemap, btreeset};
use serde_json::{Value, json};
use std::collections::{BTreeMap, VecDeque};
use std::path::PathBuf;
use treetime::seq::mutation::{Mutation, MutationEvent, MutationTrack, mutation_event_strings};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::auspice::auspice_write_file;
use treetime_io::auspice_types::{
  AuspiceColoring, AuspiceDisplayDefaults, AuspiceGenomeAnnotations, AuspiceNumDate, AuspiceTree,
  AuspiceTreeBranchAttrs, AuspiceTreeBranchAttrsLabels, AuspiceTreeData, AuspiceTreeMeta, AuspiceTreeNode,
  AuspiceTreeNodeAttr, AuspiceTreeNodeAttrs,
};
use treetime_io::graph::TreeWriteKind;
use treetime_io::graphviz::graphviz_write_file;
use treetime_io::nex::{NexWriteOptions, nex_write_file_with};
use treetime_io::nwk::{CommentProviders, NwkWriteOptions, nwk_write_file_with, nwk_write_str};
use treetime_io::usher_mat::{
  UsherMatJsonOptions, UsherMetadata, UsherMutation, UsherMutationList, UsherTree, UsherTreeNode,
  usher_mat_json_write_file, usher_mat_pb_write_file,
};
use treetime_primitives::AsciiChar;
use treetime_utils::io::json::{JsonPretty, json_write_file};
use treetime_utils::{make_error, make_internal_report};

pub(crate) const COLORING_BAD_BRANCH: &str = "bad_branch";
const COLORING_GENOTYPE: &str = "gt";
pub(crate) const COLORING_NUM_DATE: &str = "num_date";
pub(crate) const NUC_TRACK: &str = "nuc";

pub(crate) fn write_tree_outputs<A, M>(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  nwk_weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  graphviz_weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  outputs: &BTreeMap<TreeWriteKind, PathBuf>,
  providers: &CommentProviders,
  command: &str,
  to_auspice: A,
  to_mat: M,
) -> Result<(), Report>
where
  A: Fn() -> Result<AuspiceTree, Report>,
  M: Fn() -> Result<UsherTree, Report>,
{
  for (kind, path) in outputs {
    match kind {
      TreeWriteKind::Nwk(spec) => {
        graph
          .get_exactly_one_root()
          .wrap_err_with(|| format!("When converting {command} graph to Newick"))?;
        nwk_write_file_with(
          path,
          graph,
          names,
          nwk_weights,
          &NwkWriteOptions {
            style: spec.style,
            ..NwkWriteOptions::default()
          },
          providers,
        )?;
      },
      TreeWriteKind::Nexus(spec) => {
        graph
          .get_exactly_one_root()
          .wrap_err_with(|| format!("When converting {command} graph to Nexus"))?;
        nex_write_file_with(
          path,
          graph,
          names,
          nwk_weights,
          &NexWriteOptions {
            style: spec.style,
            ..NexWriteOptions::default()
          },
          providers,
        )?;
      },
      TreeWriteKind::Auspice => auspice_write_file(path, &to_auspice()?)?,
      TreeWriteKind::MatPb => usher_mat_pb_write_file(path, &to_mat()?)?,
      TreeWriteKind::MatJson => {
        usher_mat_json_write_file(path, &to_mat()?, &UsherMatJsonOptions::default())?;
      },
      TreeWriteKind::GraphJson => {
        json_write_file(path, graph, JsonPretty(true))?;
      },
      TreeWriteKind::Dot => graphviz_write_file(path, graph, names, graphviz_weights)?,
    }
  }
  Ok(())
}

#[allow(
  clippy::field_scoped_visibility_modifiers,
  reason = "crate-internal fields are the record interface"
)]
pub(crate) struct GraphNodeContext {
  pub(crate) node_key: GraphNodeKey,
  pub(crate) edge_key: Option<GraphEdgeKey>,
}

#[derive(Clone, Debug, Default, PartialEq)]
#[allow(
  clippy::field_scoped_visibility_modifiers,
  reason = "crate-internal fields are the record interface"
)]
pub(crate) struct TraitValue {
  pub(crate) value: String,
  pub(crate) confidence: BTreeMap<String, f64>,
  pub(crate) entropy: Option<f64>,
}

pub(crate) fn auspice_data(
  title: &str,
  updated: &str,
  mut colorings: Vec<AuspiceColoring>,
  filters: Vec<String>,
  color_by: Option<String>,
  genome_annotations: Option<AuspiceGenomeAnnotations>,
  root_sequences: Option<BTreeMap<String, String>>,
  has_mutations: bool,
) -> AuspiceTreeData {
  if has_mutations {
    colorings.push(coloring(COLORING_GENOTYPE, "Genotype", "categorical"));
  }
  AuspiceTreeData {
    version: Some("v2".to_owned()),
    meta: AuspiceTreeMeta {
      title: Some(title.to_owned()),
      updated: Some(updated.to_owned()),
      panels: vec!["tree".to_owned()],
      genome_annotations,
      colorings,
      filters,
      display_defaults: AuspiceDisplayDefaults {
        color_by,
        ..AuspiceDisplayDefaults::default()
      },
      ..AuspiceTreeMeta::default()
    },
    root_sequence: root_sequences.filter(|sequences| !sequences.is_empty()),
    other: Value::default(),
  }
}

pub(crate) fn auspice_node(
  name: String,
  div: Option<f64>,
  date: Option<f64>,
  date_confidence: Option<[f64; 2]>,
  bad_branch: Option<bool>,
  traits: BTreeMap<String, TraitValue>,
  mutations: BTreeMap<String, Vec<String>>,
  labels: Option<AuspiceTreeBranchAttrsLabels>,
) -> AuspiceTreeNode {
  AuspiceTreeNode {
    name,
    branch_attrs: AuspiceTreeBranchAttrs {
      mutations,
      labels,
      other: Value::default(),
    },
    node_attrs: AuspiceTreeNodeAttrs {
      div,
      num_date: date.map(|value| AuspiceNumDate {
        value,
        confidence: date_confidence,
      }),
      bad_branch: bad_branch.map(|bad| AuspiceTreeNodeAttr::new(if bad { "Yes" } else { "No" })),
      clade_membership: None,
      region: None,
      country: None,
      division: None,
      other: build_trait_attrs(traits),
    },
    children: vec![],
    other: Value::default(),
  }
}

pub(crate) fn sequence_auspice_node(
  name: &str,
  div: Option<f64>,
  confidence: Option<f64>,
  mutations: Vec<Mutation>,
  date: Option<f64>,
  bad_branch: Option<bool>,
) -> Result<AuspiceTreeNode, Report> {
  let mut other = serde_json::Map::new();
  if let Some(confidence) = confidence {
    ensure_finite(confidence, "tree output", name, "input branch support")?;
    other.insert("confidence".to_owned(), json!({ "value": confidence }));
  }
  let mut node = auspice_node(
    name.to_owned(),
    finite_number(div, 6, "tree output", name, "div")?,
    finite_number(date, 3, "tree output", name, "date")?,
    None,
    bad_branch,
    BTreeMap::new(),
    group_mutations(mutations)?,
    None,
  );
  if let Value::Object(target) = &mut node.node_attrs.other {
    target.extend(other);
  }
  Ok(node)
}

pub(crate) fn auspice_from_graph<F>(graph: &Graph, data: AuspiceTreeData, mut convert: F) -> Result<AuspiceTree, Report>
where
  F: FnMut(&GraphNodeContext) -> Result<AuspiceTreeNode, Report>,
{
  let root_key = graph
    .get_exactly_one_root()
    .wrap_err("When converting graph to Auspice v2 JSON")?
    .key();
  let mut node_map = btreemap! {};
  let mut queue = VecDeque::from([(root_key, None)]);
  while let Some((node_key, edge_key)) = queue.pop_front() {
    let converted = convert(&GraphNodeContext { node_key, edge_key })?;
    if converted.node_attrs.div.is_none() && converted.node_attrs.num_date.is_none() {
      return make_error!(
        "Auspice v2 node '{}' requires divergence or numerical date data",
        converted.name
      );
    }
    node_map.insert(node_key, converted);
    let node = graph
      .get_node(node_key)
      .ok_or_else(|| make_internal_report!("Node {node_key} not found in graph"))?;
    for (child_key, child_edge_key) in graph.children_keys_of(node) {
      queue.push_back((child_key, Some(child_edge_key)));
    }
  }
  attach_auspice_children(graph, root_key, &mut node_map)?;
  let tree = node_map
    .remove(&root_key)
    .ok_or_else(|| make_internal_report!("Auspice root node {root_key} was not converted"))?;
  Ok(AuspiceTree { data, tree })
}

fn attach_auspice_children(
  graph: &Graph,
  root_key: GraphNodeKey,
  node_map: &mut BTreeMap<GraphNodeKey, AuspiceTreeNode>,
) -> Result<(), Report> {
  let mut visited = btreeset! {};
  let mut stack = vec![root_key];
  while let Some(key) = stack.pop() {
    let node = graph
      .get_node(key)
      .ok_or_else(|| make_internal_report!("Node {key} not found in graph"))?;
    if visited.contains(&key) {
      let children = graph
        .children_keys_of(node)
        .map(|(child_key, _)| {
          node_map
            .remove(&child_key)
            .ok_or_else(|| make_internal_report!("Auspice child node {child_key} was not converted"))
        })
        .collect::<Result<Vec<_>, _>>()?;
      node_map
        .get_mut(&key)
        .ok_or_else(|| make_internal_report!("Auspice parent node {key} was not converted"))?
        .children = children;
    } else {
      visited.insert(key);
      stack.push(key);
      stack.extend(graph.children_keys_of(node).map(|(child_key, _)| child_key));
    }
  }
  Ok(())
}

pub(crate) fn mutation_free_mat(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  nwk_weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<UsherTree, Report> {
  mat_from_graph(graph, names, nwk_weights, None, |_node_key, _edge_key| Ok(vec![]))
}

pub(crate) fn mat_from_graph<F>(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  nwk_weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  reference: Option<&str>,
  mut edge_mutations: F,
) -> Result<UsherTree, Report>
where
  F: FnMut(GraphNodeKey, GraphEdgeKey) -> Result<Vec<Mutation>, Report>,
{
  graph
    .get_exactly_one_root()
    .wrap_err("When converting graph to UShER MAT")?;
  let mut node_mutations = vec![];
  let mut condensed_nodes = vec![];
  let mut metadata = vec![];
  graph.iter_depth_first_preorder_forward(|node| {
    let name = names[&node.key].clone().unwrap_or_default();
    let mutations = node
      .parent_keys
      .first()
      .map(|(_, edge_key)| edge_mutations(node.key, *edge_key))
      .transpose()?
      .unwrap_or_default();
    let mutations = mutations
      .iter()
      .map(|mutation| mat_mutation(mutation, reference, &name))
      .collect::<Result<Vec<_>, _>>()?;
    node_mutations.push(UsherMutationList { mutation: mutations });
    condensed_nodes.push(UsherTreeNode {
      node_name: name,
      condensed_leaves: vec![],
    });
    metadata.push(UsherMetadata {
      clade_annotations: vec![],
    });
    Ok(())
  })?;
  Ok(UsherTree {
    newick: nwk_write_str(graph, names, nwk_weights, &NwkWriteOptions::default())?,
    node_mutations,
    condensed_nodes,
    metadata,
  })
}

pub(crate) fn mat_mutation(
  mutation: &Mutation,
  reference: Option<&str>,
  node_name: &str,
) -> Result<UsherMutation, Report> {
  if mutation.track != MutationTrack::Nucleotide {
    return make_error!("Node '{node_name}' has an amino-acid mutation that UShER MAT cannot represent");
  }
  let MutationEvent::Substitution(substitution) = &mutation.event else {
    return make_error!("Node '{node_name}' has an insertion or deletion that UShER MAT cannot represent");
  };
  let reference = reference.ok_or_else(|| {
    eyre::eyre!("Node '{node_name}' has nucleotide mutations, but UShER MAT requires a root nucleotide reference")
  })?;
  let position = substitution
    .pos()
    .checked_add(1)
    .ok_or_else(|| eyre::eyre!("Node '{node_name}' mutation coordinate overflow"))?;
  let position = i32::try_from(position).wrap_err_with(|| {
    format!("Node '{node_name}' mutation position {position} exceeds the UShER MAT i32 coordinate range")
  })?;
  let reference_state = reference.as_bytes().get(substitution.pos()).copied().ok_or_else(|| {
    eyre::eyre!(
      "Node '{node_name}' mutation position {} is outside the root nucleotide reference of length {}",
      substitution.pos() + 1,
      reference.len()
    )
  })?;
  let reference_state = AsciiChar::try_new(reference_state)?;
  Ok(UsherMutation {
    position,
    ref_nuc: mat_nucleotide(reference_state, node_name, "root reference")?,
    par_nuc: mat_nucleotide(substitution.reff(), node_name, "parent")?,
    mut_nuc: vec![mat_nucleotide(substitution.qry(), node_name, "child")?],
    chromosome: String::new(),
  })
}

fn mat_nucleotide(nucleotide: AsciiChar, node_name: &str, role: &str) -> Result<i32, Report> {
  match char::from(nucleotide).to_ascii_uppercase() {
    'A' => Ok(0),
    'C' => Ok(1),
    'G' => Ok(2),
    'T' => Ok(3),
    state => {
      make_error!("Node '{node_name}' has {role} nucleotide '{state}', but UShER MAT accepts only A, C, G, or T")
    },
  }
}

pub(crate) fn group_mutations(mutations: Vec<Mutation>) -> Result<BTreeMap<String, Vec<String>>, Report> {
  let mut grouped = BTreeMap::new();
  for mutation in mutations {
    let track = match &mutation.track {
      MutationTrack::Nucleotide => NUC_TRACK.to_owned(),
      MutationTrack::AminoAcid(track) => {
        if track.is_empty()
          || !track
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'*' | b'_' | b'.' | b'(' | b')' | b'-'))
        {
          return make_error!("Auspice v2 cannot represent amino-acid mutation track '{track}'");
        }
        track.clone()
      },
    };
    if matches!(mutation.track, MutationTrack::Nucleotide) && !matches!(mutation.event, MutationEvent::Substitution(_))
    {
      continue;
    }
    grouped
      .entry(track)
      .or_insert_with(Vec::new)
      .extend(mutation_event_strings(&mutation.event)?);
  }
  Ok(grouped)
}

fn build_trait_attrs(traits: BTreeMap<String, TraitValue>) -> Value {
  Value::Object(
    traits
      .into_iter()
      .map(|(attribute, value)| {
        let mut fields = serde_json::Map::new();
        fields.insert("value".to_owned(), json!(value.value));
        if !value.confidence.is_empty() {
          fields.insert(
            "confidence".to_owned(),
            json!(
              value
                .confidence
                .iter()
                .map(|(state, probability)| (state.clone(), format_number(*probability, 3)))
                .collect::<BTreeMap<_, _>>()
            ),
          );
        }
        if let Some(entropy) = value.entropy {
          fields.insert("entropy".to_owned(), json!(format_number(entropy, 3)));
        }
        (attribute, Value::Object(fields))
      })
      .collect(),
  )
}

pub(crate) fn cumulative_branch_length_from(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  mut key: GraphNodeKey,
) -> Result<Option<f64>, Report> {
  let mut total = 0.0;
  while let Some((parent, edge_key)) = graph.node_parent(key)? {
    let Some(length) = branch_lengths[&edge_key] else {
      return Ok(None);
    };
    total += length;
    key = parent;
  }
  Ok(Some(total))
}

pub(crate) fn node_name_value(key: GraphNodeKey, name: Option<&str>) -> String {
  name.map_or_else(|| format!("node_{}", key.as_usize()), str::to_owned)
}

pub(crate) fn ensure_optional_finite(
  value: Option<f64>,
  command: &str,
  node_name: &str,
  field: &str,
) -> Result<(), Report> {
  if let Some(value) = value {
    ensure_finite(value, command, node_name, field)?;
  }
  Ok(())
}

pub(crate) fn ensure_finite(value: f64, command: &str, node_name: &str, field: &str) -> Result<(), Report> {
  if value.is_finite() {
    Ok(())
  } else {
    make_error!("{command} node '{node_name}' has non-finite {field}={value}")
  }
}

pub(crate) fn finite_number(
  value: Option<f64>,
  precision: i32,
  command: &str,
  node_name: &str,
  field: &str,
) -> Result<Option<f64>, Report> {
  ensure_optional_finite(value, command, node_name, field)?;
  Ok(value.map(|value| format_number(value, precision)))
}

#[allow(
  clippy::as_conversions,
  clippy::expect_used,
  reason = "count/index numeric cast is exact for the domain range; expect on a value an upstream invariant guarantees is present"
)]
pub(crate) fn format_number(number: f64, precision: i32) -> f64 {
  if number == 0.0 || !number.is_finite() {
    return number;
  }
  let integral = number.abs().trunc();
  let significand = if integral >= 1.0 {
    integral.log10().floor() as i32 + 1
  } else {
    0
  };
  let significant_figures = (significand + precision).max(1) as usize;
  format!("{number:.*e}", significant_figures - 1)
    .parse()
    .expect("a float formatted in scientific notation must parse back")
}

pub(crate) fn coloring(key: &str, title: &str, type_: &str) -> AuspiceColoring {
  AuspiceColoring {
    key: key.to_owned(),
    title: title.to_owned(),
    type_: type_.to_owned(),
    ..AuspiceColoring::default()
  }
}

pub(crate) fn generation_date() -> String {
  Utc::now().format("%Y-%m-%d").to_string()
}
