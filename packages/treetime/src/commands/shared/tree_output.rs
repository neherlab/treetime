use crate::ancestral::pipeline::AncestralPartition;
use crate::clock::clock_graph::GraphClock;
use crate::commands::ancestral::result::AncestralGraphData;
use crate::commands::clock::run::ClockGraphData;
use crate::commands::mugration::augur_node_data::{build_confidence_map, compute_entropy};
use crate::commands::optimize::result::OptimizeGraphData;
use crate::commands::prune::result::PruneGraphData;
use crate::commands::timetree::result::TimetreeGraphData;
use crate::mugration::result::MugrationGraphData;
use crate::partition::traits::{BranchTopology, PartitionBranchOps};
use crate::payload::ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use crate::payload::timetree::{EdgeTimetree, NodeTimetree};
use crate::seq::mutation::{Mutation, MutationEvent, MutationTrack, mutation_event_strings};
use chrono::Utc;
use eyre::{Report, WrapErr};
use maplit::{btreemap, btreeset};
use parking_lot::RwLock;
use percent_encoding::{AsciiSet, NON_ALPHANUMERIC, utf8_percent_encode};
use serde::Serialize;
use serde_json::{Value, json};
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::path::PathBuf;
use std::sync::Arc;
use treetime_graph::edge::{Edge, GraphEdge, GraphEdgeKey, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey, Named, Node};
use treetime_io::auspice::auspice_write_file;
use treetime_io::auspice_types::{
  AuspiceColoring, AuspiceDisplayDefaults, AuspiceGenomeAnnotationCds, AuspiceGenomeAnnotationNuc,
  AuspiceGenomeAnnotations, AuspiceNumDate, AuspiceTree, AuspiceTreeBranchAttrs, AuspiceTreeBranchAttrsLabels,
  AuspiceTreeData, AuspiceTreeMeta, AuspiceTreeNode, AuspiceTreeNodeAttr, AuspiceTreeNodeAttrs, Segments, StartEnd,
};
use treetime_io::graph::TreeWriteKind;
use treetime_io::graphviz::{EdgeToGraphviz, NodeToGraphviz, graphviz_write_file};
use treetime_io::nex::{NexWriteOptions, nex_write_file_with};
use treetime_io::nwk::{CommentProviders, EdgeToNwk, NodeToNwk, NwkWriteOptions, nwk_write_file_with, nwk_write_str};
use treetime_io::phyloxml::{
  Phyloxml, PhyloxmlClade, PhyloxmlConfidence, PhyloxmlDate, PhyloxmlJsonOptions, PhyloxmlMolSeq, PhyloxmlPhylogeny,
  PhyloxmlProperty, PhyloxmlSequence, phyloxml_json_write_file, phyloxml_write_file,
};
use treetime_io::usher_mat::{
  UsherMatJsonOptions, UsherMetadata, UsherMutation, UsherMutationList, UsherTree, UsherTreeNode,
  usher_mat_json_write_file, usher_mat_pb_write_file,
};
use treetime_primitives::AsciiChar;
use treetime_utils::io::json::{JsonPretty, json_write_file};
use treetime_utils::{make_error, make_internal_report, make_report};
use util_augur_node_data_json::AugurNodeDataJsonAnnotationEntry;

const APPLIES_BRANCH: &str = "parent_branch";
const APPLIES_NODE: &str = "node";
const BRANCH_LENGTH_UNIT: &str = "subs/site";
const COLORING_BAD_BRANCH: &str = "bad_branch";
const COLORING_GENOTYPE: &str = "gt";
const COLORING_NUM_DATE: &str = "num_date";
const DT_BOOLEAN: &str = "xsd:boolean";
const DT_DOUBLE: &str = "xsd:double";
const DT_STRING: &str = "xsd:string";
const NUC_TRACK: &str = "nuc";
const REF_BAD_BRANCH: &str = "treetime:bad_branch";
const REF_DATE_INFERRED: &str = "treetime:date_inferred";
const REF_DIV: &str = "treetime:divergence";
const REF_GAMMA: &str = "treetime:gamma";
const REF_MUTATION: &str = "treetime:mutation";
const REF_TRAIT_CONFIDENCE_PREFIX: &str = "treetime:trait_confidence:";
const REF_TRAIT_ENTROPY_PREFIX: &str = "treetime:trait_entropy:";
const REF_TRAIT_PREFIX: &str = "treetime:trait:";
const REF_TRAIT_TRANSITION_PREFIX: &str = "treetime:trait_transition:";
const TYPE_INPUT_BRANCH_SUPPORT: &str = "treetime:input_branch_support";
const PROPERTY_TOKEN_ENCODE_SET: &AsciiSet = &NON_ALPHANUMERIC.remove(b'-').remove(b'.').remove(b'_').remove(b'~');

pub fn write_ancestral_tree_outputs(
  graph: &GraphAncestral<AncestralGraphData>,
  outputs: &BTreeMap<TreeWriteKind, PathBuf>,
  providers: &CommentProviders,
) -> Result<(), Report> {
  let updated = generation_date();
  write_tree_outputs(
    graph,
    outputs,
    providers,
    "ancestral",
    || ancestral_to_auspice(graph, &updated),
    || ancestral_to_phyloxml(graph),
    || ancestral_to_mat(graph),
  )
}

pub fn write_optimize_tree_outputs(
  graph: &GraphAncestral<OptimizeGraphData>,
  outputs: &BTreeMap<TreeWriteKind, PathBuf>,
  providers: &CommentProviders,
) -> Result<(), Report> {
  let updated = generation_date();
  write_tree_outputs(
    graph,
    outputs,
    providers,
    "optimize",
    || optimize_to_auspice(graph, &updated),
    || optimize_to_phyloxml(graph),
    || optimize_to_mat(graph),
  )
}

pub fn write_prune_tree_outputs(
  graph: &GraphAncestral<PruneGraphData>,
  outputs: &BTreeMap<TreeWriteKind, PathBuf>,
  providers: &CommentProviders,
) -> Result<(), Report> {
  let updated = generation_date();
  write_tree_outputs(
    graph,
    outputs,
    providers,
    "prune",
    || prune_to_auspice(graph, &updated),
    || prune_to_phyloxml(graph),
    || prune_to_mat(graph),
  )
}

pub fn write_clock_tree_outputs(
  graph: &GraphClock<ClockGraphData>,
  outputs: &BTreeMap<TreeWriteKind, PathBuf>,
  providers: &CommentProviders,
) -> Result<(), Report> {
  let updated = generation_date();
  write_tree_outputs(
    graph,
    outputs,
    providers,
    "clock",
    || clock_to_auspice(graph, &updated),
    || clock_to_phyloxml(graph),
    || clock_to_mat(graph),
  )
}

pub fn write_mugration_tree_outputs(
  graph: &GraphAncestral<MugrationGraphData>,
  outputs: &BTreeMap<TreeWriteKind, PathBuf>,
  providers: &CommentProviders,
) -> Result<(), Report> {
  let updated = generation_date();
  write_tree_outputs(
    graph,
    outputs,
    providers,
    "mugration",
    || mugration_to_auspice(graph, &updated),
    || mugration_to_phyloxml(graph),
    || mugration_to_mat(graph),
  )
}

pub fn write_timetree_tree_outputs(
  graph: &Graph<NodeTimetree, EdgeTimetree, TimetreeGraphData>,
  outputs: &BTreeMap<TreeWriteKind, PathBuf>,
  providers: &CommentProviders,
) -> Result<(), Report> {
  let updated = generation_date();
  write_tree_outputs(
    graph,
    outputs,
    providers,
    "timetree",
    || timetree_to_auspice(graph, &updated),
    || timetree_to_phyloxml(graph),
    || timetree_to_mat(graph),
  )
}

fn write_tree_outputs<N, E, D, A, P, M>(
  graph: &Graph<N, E, D>,
  outputs: &BTreeMap<TreeWriteKind, PathBuf>,
  providers: &CommentProviders,
  command: &str,
  to_auspice: A,
  to_phyloxml: P,
  to_mat: M,
) -> Result<(), Report>
where
  N: GraphNode + Named + NodeToNwk + NodeToGraphviz + Serialize,
  E: GraphEdge + HasBranchLength + EdgeToNwk + EdgeToGraphviz + Serialize,
  D: Send + Sync + Serialize,
  A: Fn() -> Result<AuspiceTree, Report>,
  P: Fn() -> Result<Phyloxml, Report>,
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
          &NexWriteOptions {
            style: spec.style,
            ..NexWriteOptions::default()
          },
          providers,
        )?;
      },
      TreeWriteKind::Auspice => auspice_write_file(path, &to_auspice()?)?,
      TreeWriteKind::Phyloxml => phyloxml_write_file(path, &to_phyloxml()?)?,
      TreeWriteKind::PhyloxmlJson => {
        phyloxml_json_write_file(path, &to_phyloxml()?, &PhyloxmlJsonOptions::default())?;
      },
      TreeWriteKind::MatPb => usher_mat_pb_write_file(path, &to_mat()?)?,
      TreeWriteKind::MatJson => {
        usher_mat_json_write_file(path, &to_mat()?, &UsherMatJsonOptions::default())?;
      },
      TreeWriteKind::GraphJson => {
        // Graph JSON is an unstable debug dump of the concrete internal structs. It has no portable schema or input contract.
        json_write_file(path, graph, JsonPretty(true))?;
      },
      TreeWriteKind::Dot => graphviz_write_file(path, graph)?,
    }
  }
  Ok(())
}

pub(crate) fn ancestral_to_auspice(
  graph: &GraphAncestral<AncestralGraphData>,
  updated: &str,
) -> Result<AuspiceTree, Report> {
  let root_sequences = ancestral_root_sequences(graph)?;
  let genome_annotations = ancestral_genome_annotations(graph, &root_sequences)?;
  let data = auspice_data(
    "TreeTime ancestral analysis",
    updated,
    vec![],
    vec![],
    None,
    genome_annotations,
    Some(root_sequences),
    !ancestral_all_mutations(graph)?.is_empty(),
  );
  auspice_from_graph(graph, data, |context| {
    let mutations = ancestral_node_mutations(graph, context.node_key, context.edge_key)?;
    ancestral_auspice_node(context, mutations, None, None)
  })
}

pub(crate) fn optimize_to_auspice(
  graph: &GraphAncestral<OptimizeGraphData>,
  updated: &str,
) -> Result<AuspiceTree, Report> {
  let root_sequences = optimize_root_sequences(graph)?;
  let data = auspice_data(
    "TreeTime optimize analysis",
    updated,
    vec![],
    vec![],
    None,
    None,
    Some(root_sequences),
    optimize_has_mutations(graph),
  );
  auspice_from_graph(graph, data, |context| {
    ancestral_auspice_node(context, optimize_mutations(graph, context.edge_key)?, None, None)
  })
}

pub(crate) fn prune_to_auspice(graph: &GraphAncestral<PruneGraphData>, updated: &str) -> Result<AuspiceTree, Report> {
  let root_sequences = prune_root_sequences(graph)?;
  let data = auspice_data(
    "TreeTime prune analysis",
    updated,
    vec![],
    vec![],
    None,
    None,
    Some(root_sequences),
    !graph.data().partitions.is_empty(),
  );
  auspice_from_graph(graph, data, |context| {
    ancestral_auspice_node(context, prune_mutations(graph, context.edge_key)?, None, None)
  })
}

pub(crate) fn clock_to_auspice(graph: &GraphClock<ClockGraphData>, updated: &str) -> Result<AuspiceTree, Report> {
  let data = auspice_data(
    "TreeTime clock analysis",
    updated,
    vec![
      coloring(COLORING_NUM_DATE, "Date", "continuous"),
      coloring(COLORING_BAD_BRANCH, "Excluded", "categorical"),
    ],
    vec![COLORING_BAD_BRANCH.to_owned()],
    Some(COLORING_BAD_BRANCH.to_owned()),
    None,
    None,
    false,
  );
  auspice_from_graph(graph, data, |context| {
    let name = node_name(context.node_key, context.node);
    Ok(auspice_node(
      name.clone(),
      finite_number(Some(context.node.div), 6, "clock", &name, "div")?,
      finite_number(context.node.time, 3, "clock", &name, "date")?,
      None,
      Some(context.node.bad_branch || context.node.is_outlier),
      BTreeMap::new(),
      BTreeMap::new(),
      None,
    ))
  })
}

pub(crate) fn mugration_to_auspice(
  graph: &GraphAncestral<MugrationGraphData>,
  updated: &str,
) -> Result<AuspiceTree, Report> {
  let attribute = &graph.data().traits.attribute;
  let data = auspice_data(
    "TreeTime mugration analysis",
    updated,
    vec![coloring(attribute, attribute, "categorical")],
    vec![attribute.clone()],
    Some(attribute.clone()),
    None,
    None,
    false,
  );
  auspice_from_graph(graph, data, |context| {
    let name = node_name(context.node_key, context.node);
    let traits = mugration_traits(graph, context.node_key, &name)?;
    Ok(auspice_node(
      name.clone(),
      finite_number(
        cumulative_branch_length(graph, context.node_key)?,
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
      mugration_transition_label(graph, context.node_key)?,
    ))
  })
}

pub(crate) fn timetree_to_auspice(
  graph: &Graph<NodeTimetree, EdgeTimetree, TimetreeGraphData>,
  updated: &str,
) -> Result<AuspiceTree, Report> {
  let root_sequences = timetree_root_sequences(graph)?;
  let data = auspice_data(
    "TreeTime timetree analysis",
    updated,
    vec![
      coloring(COLORING_NUM_DATE, "Date", "continuous"),
      coloring(COLORING_BAD_BRANCH, "Excluded", "categorical"),
    ],
    vec![COLORING_BAD_BRANCH.to_owned()],
    Some(COLORING_BAD_BRANCH.to_owned()),
    None,
    Some(root_sequences),
    !graph.data().partitions.is_empty(),
  );
  auspice_from_graph(graph, data, |context| {
    let name = node_name(context.node_key, context.node);
    let div = timetree_divergence(graph, context.node_key, context.node)?;
    let confidence = timetree_date_confidence(graph, context.node_key, &name)?;
    Ok(auspice_node(
      name.clone(),
      finite_number(Some(div), 6, "timetree", &name, "div")?,
      finite_number(context.node.time, 3, "timetree", &name, "date")?,
      confidence,
      Some(context.node.bad_branch || context.node.is_outlier),
      BTreeMap::new(),
      group_mutations(timetree_mutations(graph, context.edge_key)?)?,
      None,
    ))
  })
}

fn ancestral_auspice_node<D: Send + Sync>(
  context: &GraphNodeContext<NodeAncestral, EdgeAncestral, D>,
  mutations: Vec<Mutation>,
  date: Option<f64>,
  bad_branch: Option<bool>,
) -> Result<AuspiceTreeNode, Report> {
  let name = node_name(context.node_key, context.node);
  let div = cumulative_branch_length(context.graph, context.node_key)?;
  let mut other = serde_json::Map::new();
  if let Some(confidence) = context.node.confidence {
    ensure_finite(confidence, "tree output", &name, "input branch support")?;
    other.insert("confidence".to_owned(), json!({ "value": confidence }));
  }
  let mut node = auspice_node(
    name.clone(),
    finite_number(div, 6, "tree output", &name, "div")?,
    finite_number(date, 3, "tree output", &name, "date")?,
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

fn auspice_data(
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

fn auspice_node(
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

fn auspice_from_graph<N, E, D, F>(
  graph: &Graph<N, E, D>,
  data: AuspiceTreeData,
  mut convert: F,
) -> Result<AuspiceTree, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Send + Sync,
  F: FnMut(&GraphNodeContext<N, E, D>) -> Result<AuspiceTreeNode, Report>,
{
  let root = graph
    .get_exactly_one_root()
    .wrap_err("When converting graph to Auspice v2 JSON")?;
  let mut node_map = btreemap! {};
  let mut queue = VecDeque::from([(Arc::clone(&root), None)]);
  while let Some((current_node, current_edge)) = queue.pop_front() {
    let node_key = current_node.read_arc().key();
    let node = current_node.read_arc().payload().read_arc();
    let edge_key = current_edge
      .as_ref()
      .map(|edge: &Arc<RwLock<Edge<E>>>| edge.read_arc().key());
    let edge = current_edge
      .as_ref()
      .map(|edge: &Arc<RwLock<Edge<E>>>| edge.read_arc().payload().read_arc());
    let converted = convert(&GraphNodeContext {
      node_key,
      node: &node,
      edge_key,
      edge: edge.as_deref(),
      graph,
    })?;
    if converted.node_attrs.div.is_none() && converted.node_attrs.num_date.is_none() {
      return make_error!(
        "Auspice v2 node '{}' requires divergence or numerical date data",
        converted.name
      );
    }
    node_map.insert(node_key, converted);
    for (child, edge) in graph.children_of(&current_node.read_arc()) {
      queue.push_back((child, Some(edge)));
    }
  }
  attach_auspice_children(graph, &root, &mut node_map)?;
  let root_key = root.read_arc().key();
  let tree = node_map
    .remove(&root_key)
    .ok_or_else(|| make_internal_report!("Auspice root node {root_key} was not converted"))?;
  Ok(AuspiceTree { data, tree })
}

fn attach_auspice_children<N, E, D>(
  graph: &Graph<N, E, D>,
  root: &Arc<RwLock<Node<N>>>,
  node_map: &mut BTreeMap<GraphNodeKey, AuspiceTreeNode>,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Send + Sync,
{
  let mut visited = btreeset! {};
  let mut stack = vec![Arc::clone(root)];
  while let Some(node) = stack.pop() {
    let key = node.read_arc().key();
    if visited.contains(&key) {
      let children = graph
        .children_of(&node.read_arc())
        .into_iter()
        .map(|(child, _)| {
          let child_key = child.read_arc().key();
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
      stack.push(Arc::clone(&node));
      stack.extend(graph.children_of(&node.read_arc()).into_iter().map(|(child, _)| child));
    }
  }
  Ok(())
}

struct GraphNodeContext<'a, N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: Send + Sync,
{
  node_key: GraphNodeKey,
  node: &'a N,
  edge_key: Option<GraphEdgeKey>,
  edge: Option<&'a E>,
  graph: &'a Graph<N, E, D>,
}

#[derive(Clone, Debug, Default, PartialEq)]
struct TraitValue {
  value: String,
  confidence: BTreeMap<String, f64>,
  entropy: Option<f64>,
}

pub(crate) fn ancestral_to_phyloxml(graph: &GraphAncestral<AncestralGraphData>) -> Result<Phyloxml, Report> {
  phyloxml_from_graph(graph, "TreeTime ancestral analysis", |context| {
    let mutations = ancestral_node_mutations(graph, context.node_key, context.edge_key)?;
    ancestral_phyloxml_clade(context, mutations, &ancestral_node_sequences(graph, context.node_key))
  })
}

pub(crate) fn optimize_to_phyloxml(graph: &GraphAncestral<OptimizeGraphData>) -> Result<Phyloxml, Report> {
  phyloxml_from_graph(graph, "TreeTime optimize analysis", |context| {
    ancestral_phyloxml_clade(
      context,
      optimize_mutations(graph, context.edge_key)?,
      &optimize_node_sequences(graph, context.node_key),
    )
  })
}

pub(crate) fn prune_to_phyloxml(graph: &GraphAncestral<PruneGraphData>) -> Result<Phyloxml, Report> {
  phyloxml_from_graph(graph, "TreeTime prune analysis", |context| {
    ancestral_phyloxml_clade(
      context,
      prune_mutations(graph, context.edge_key)?,
      &prune_node_sequences(graph, context.node_key),
    )
  })
}

pub(crate) fn clock_to_phyloxml(graph: &GraphClock<ClockGraphData>) -> Result<Phyloxml, Report> {
  phyloxml_from_graph(graph, "TreeTime clock analysis", |context| {
    let name = node_name(context.node_key, context.node);
    ensure_optional_finite(context.node.time, "clock", &name, "date")?;
    ensure_finite(context.node.div, "clock", &name, "divergence")?;
    let property = vec![
      property(REF_DIV, DT_DOUBLE, APPLIES_NODE, &context.node.div.to_string()),
      property(
        REF_BAD_BRANCH,
        DT_BOOLEAN,
        APPLIES_NODE,
        if context.node.bad_branch || context.node.is_outlier {
          "true"
        } else {
          "false"
        },
      ),
    ];
    Ok(PhyloxmlClade {
      name: context.node.name.clone(),
      branch_length_elem: context.edge.and_then(HasBranchLength::branch_length),
      branch_length_attr: None,
      confidence: vec![],
      width: None,
      color: None,
      node_id: None,
      taxonomy: vec![],
      sequence: vec![],
      events: None,
      binary_characters: None,
      distribution: vec![],
      date: context.node.time.map(phyloxml_date),
      reference: vec![],
      property,
      clade: vec![],
      other: BTreeMap::new(),
    })
  })
}

pub(crate) fn mugration_to_phyloxml(graph: &GraphAncestral<MugrationGraphData>) -> Result<Phyloxml, Report> {
  phyloxml_from_graph(graph, "TreeTime mugration analysis", |context| {
    let name = node_name(context.node_key, context.node);
    let div = cumulative_branch_length(graph, context.node_key)?;
    ensure_optional_finite(div, "mugration", &name, "divergence")?;
    let traits = mugration_traits(graph, context.node_key, &name)?;
    let mut properties = div
      .map(|div| vec![property(REF_DIV, DT_DOUBLE, APPLIES_NODE, &div.to_string())])
      .unwrap_or_default();
    properties.extend(trait_properties(&traits, &name)?);
    if let Some((parent, child)) = mugration_transition(graph, context.node_key)? {
      let attribute = encode_property_token(&graph.data().traits.attribute);
      properties.push(property(
        &format!("{REF_TRAIT_TRANSITION_PREFIX}{attribute}"),
        DT_STRING,
        APPLIES_BRANCH,
        &format!("{}:{}", encode_property_token(&parent), encode_property_token(&child)),
      ));
    }
    let mut clade = empty_phyloxml_clade(
      context.node.name.clone(),
      context.edge.and_then(HasBranchLength::branch_length),
    );
    clade.confidence = input_branch_confidence(context.node.confidence, "mugration", &name)?;
    clade.property = properties;
    Ok(clade)
  })
}

pub(crate) fn timetree_to_phyloxml(
  graph: &Graph<NodeTimetree, EdgeTimetree, TimetreeGraphData>,
) -> Result<Phyloxml, Report> {
  phyloxml_from_graph(graph, "TreeTime timetree analysis", |context| {
    let name = node_name(context.node_key, context.node);
    let divergence = timetree_divergence(graph, context.node_key, context.node)?;
    ensure_finite(divergence, "timetree", &name, "divergence")?;
    ensure_optional_finite(context.node.time, "timetree", &name, "date")?;
    ensure_finite(context.edge.map_or(1.0, |edge| edge.gamma), "timetree", &name, "gamma")?;
    let mut properties = vec![
      property(REF_DIV, DT_DOUBLE, APPLIES_NODE, &divergence.to_string()),
      property(
        REF_BAD_BRANCH,
        DT_BOOLEAN,
        APPLIES_NODE,
        if context.node.bad_branch || context.node.is_outlier {
          "true"
        } else {
          "false"
        },
      ),
    ];
    if timetree_date_is_inferred(graph, context.node_key, context.node) == Some(true) {
      properties.push(property(REF_DATE_INFERRED, DT_BOOLEAN, APPLIES_NODE, "true"));
    }
    if let Some(edge) = context.edge {
      properties.push(property(REF_GAMMA, DT_DOUBLE, APPLIES_BRANCH, &edge.gamma.to_string()));
    }
    for mutation in timetree_mutations(graph, context.edge_key)? {
      properties.push(mutation_property(&mutation)?);
    }
    let date = context.node.time.map(|value| {
      let confidence = graph
        .data()
        .confidence_intervals
        .as_ref()
        .and_then(|intervals| intervals.iter().find(|interval| interval.key == context.node_key));
      PhyloxmlDate {
        desc: None,
        value: Some(value),
        minimum: confidence.map(|interval| interval.lower),
        maximum: confidence.map(|interval| interval.upper),
        unit: Some("year".to_owned()),
      }
    });
    let mut clade = empty_phyloxml_clade(
      context.node.base.name.clone(),
      context.edge.and_then(HasBranchLength::branch_length),
    );
    clade.confidence = input_branch_confidence(context.node.base.confidence, "timetree", &name)?;
    clade.date = date;
    clade.property = properties;
    clade.sequence = phyloxml_sequences(&timetree_node_sequences(graph, context.node_key));
    Ok(clade)
  })
}

fn ancestral_phyloxml_clade<D: Send + Sync>(
  context: &GraphNodeContext<NodeAncestral, EdgeAncestral, D>,
  mutations: Vec<Mutation>,
  sequences: &BTreeMap<String, String>,
) -> Result<PhyloxmlClade, Report> {
  let name = node_name(context.node_key, context.node);
  let divergence = cumulative_branch_length(context.graph, context.node_key)?;
  ensure_optional_finite(divergence, "tree output", &name, "divergence")?;
  let mut properties = divergence
    .map(|divergence| vec![property(REF_DIV, DT_DOUBLE, APPLIES_NODE, &divergence.to_string())])
    .unwrap_or_default();
  properties.extend(
    mutations
      .into_iter()
      .map(|mutation| mutation_property(&mutation))
      .collect::<Result<Vec<_>, _>>()?,
  );
  let mut clade = empty_phyloxml_clade(
    context.node.name.clone(),
    context.edge.and_then(HasBranchLength::branch_length),
  );
  clade.confidence = input_branch_confidence(context.node.confidence, "tree output", &name)?;
  clade.property = properties;
  clade.sequence = phyloxml_sequences(sequences);
  Ok(clade)
}

fn phyloxml_from_graph<N, E, D, F>(graph: &Graph<N, E, D>, title: &str, mut convert: F) -> Result<Phyloxml, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Send + Sync,
  F: FnMut(&GraphNodeContext<N, E, D>) -> Result<PhyloxmlClade, Report>,
{
  let root = graph
    .get_exactly_one_root()
    .wrap_err("When converting graph to PhyloXML")?;
  let mut node_map = btreemap! {};
  let mut queue = VecDeque::from([(Arc::clone(&root), None)]);
  while let Some((current_node, current_edge)) = queue.pop_front() {
    let node_key = current_node.read_arc().key();
    let node = current_node.read_arc().payload().read_arc();
    let edge_key = current_edge
      .as_ref()
      .map(|edge: &Arc<RwLock<Edge<E>>>| edge.read_arc().key());
    let edge = current_edge
      .as_ref()
      .map(|edge: &Arc<RwLock<Edge<E>>>| edge.read_arc().payload().read_arc());
    node_map.insert(
      node_key,
      convert(&GraphNodeContext {
        node_key,
        node: &node,
        edge_key,
        edge: edge.as_deref(),
        graph,
      })?,
    );
    for (child, edge) in graph.children_of(&current_node.read_arc()) {
      queue.push_back((child, Some(edge)));
    }
  }
  attach_phyloxml_children(graph, &root, &mut node_map)?;
  let root_key = root.read_arc().key();
  let clade = node_map
    .remove(&root_key)
    .ok_or_else(|| make_internal_report!("PhyloXML root node {root_key} was not converted"))?;
  Ok(Phyloxml {
    phylogeny: vec![PhyloxmlPhylogeny {
      rooted: true,
      rerootable: None,
      branch_length_unit: Some(BRANCH_LENGTH_UNIT.to_owned()),
      phylogeny_type: None,
      name: Some(title.to_owned()),
      id: None,
      description: None,
      date: None,
      confidence: vec![],
      clade: Some(clade),
      clade_relation: vec![],
      sequence_relation: vec![],
      property: vec![],
      other: BTreeMap::new(),
    }],
    other: BTreeMap::new(),
  })
}

fn attach_phyloxml_children<N, E, D>(
  graph: &Graph<N, E, D>,
  root: &Arc<RwLock<Node<N>>>,
  node_map: &mut BTreeMap<GraphNodeKey, PhyloxmlClade>,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Send + Sync,
{
  let mut visited = BTreeSet::new();
  let mut stack = vec![Arc::clone(root)];
  while let Some(node) = stack.pop() {
    let key = node.read_arc().key();
    if visited.contains(&key) {
      let children = graph
        .children_of(&node.read_arc())
        .into_iter()
        .map(|(child, _)| {
          let child_key = child.read_arc().key();
          node_map
            .remove(&child_key)
            .ok_or_else(|| make_internal_report!("PhyloXML child node {child_key} was not converted"))
        })
        .collect::<Result<Vec<_>, _>>()?;
      node_map
        .get_mut(&key)
        .ok_or_else(|| make_internal_report!("PhyloXML parent node {key} was not converted"))?
        .clade = children;
    } else {
      visited.insert(key);
      stack.push(Arc::clone(&node));
      stack.extend(graph.children_of(&node.read_arc()).into_iter().map(|(child, _)| child));
    }
  }
  Ok(())
}

pub(crate) fn ancestral_to_mat(graph: &GraphAncestral<AncestralGraphData>) -> Result<UsherTree, Report> {
  let reference = ancestral_root_sequences(graph)?.remove(NUC_TRACK);
  mat_from_graph(graph, reference.as_deref(), |node_key, edge_key| {
    ancestral_node_mutations(graph, node_key, Some(edge_key))
  })
}

pub(crate) fn optimize_to_mat(graph: &GraphAncestral<OptimizeGraphData>) -> Result<UsherTree, Report> {
  let reference = optimize_root_sequences(graph)?.remove(NUC_TRACK);
  mat_from_graph(graph, reference.as_deref(), |_node_key, edge_key| {
    optimize_mutations(graph, Some(edge_key))
  })
}

pub(crate) fn prune_to_mat(graph: &GraphAncestral<PruneGraphData>) -> Result<UsherTree, Report> {
  let reference = prune_root_sequences(graph)?.remove(NUC_TRACK);
  mat_from_graph(graph, reference.as_deref(), |_node_key, edge_key| {
    prune_mutations(graph, Some(edge_key))
  })
}

pub(crate) fn clock_to_mat(graph: &GraphClock<ClockGraphData>) -> Result<UsherTree, Report> {
  mutation_free_mat(graph)
}

pub(crate) fn mugration_to_mat(graph: &GraphAncestral<MugrationGraphData>) -> Result<UsherTree, Report> {
  mutation_free_mat(graph)
}

pub(crate) fn timetree_to_mat(
  graph: &Graph<NodeTimetree, EdgeTimetree, TimetreeGraphData>,
) -> Result<UsherTree, Report> {
  let reference = timetree_root_sequences(graph)?.remove(NUC_TRACK);
  mat_from_graph(graph, reference.as_deref(), |_node_key, edge_key| {
    timetree_mutations(graph, Some(edge_key))
  })
}

fn mutation_free_mat<N, E, D>(graph: &Graph<N, E, D>) -> Result<UsherTree, Report>
where
  N: GraphNode + Named + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Send + Sync,
{
  mat_from_graph(graph, None, |_node_key, _edge_key| Ok(vec![]))
}

fn mat_from_graph<N, E, D, F>(
  graph: &Graph<N, E, D>,
  reference: Option<&str>,
  mut edge_mutations: F,
) -> Result<UsherTree, Report>
where
  N: GraphNode + Named + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Send + Sync,
  F: FnMut(GraphNodeKey, GraphEdgeKey) -> Result<Vec<Mutation>, Report>,
{
  graph
    .get_exactly_one_root()
    .wrap_err("When converting graph to UShER MAT")?;
  let mut node_mutations = vec![];
  let mut condensed_nodes = vec![];
  let mut metadata = vec![];
  graph.iter_depth_first_preorder_forward(|node| {
    let name = node
      .payload
      .name()
      .map(|name| name.as_ref().to_owned())
      .unwrap_or_default();
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
    newick: nwk_write_str(graph, &NwkWriteOptions::default())?,
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

fn ancestral_root_sequences(graph: &GraphAncestral<AncestralGraphData>) -> Result<BTreeMap<String, String>, Report> {
  let mut sequences = BTreeMap::new();
  if let Some(partition) = graph.data().partition.as_ref() {
    let sequence = match partition {
      AncestralPartition::Fitch(partition) => partition.read_arc().root_sequence(graph)?,
      AncestralPartition::Sparse(partition) => partition.read_arc().root_sequence(graph)?,
      AncestralPartition::Dense(partition) => partition.read_arc().root_sequence(graph)?,
    };
    sequences.insert(NUC_TRACK.to_owned(), sequence.to_string());
  }
  if let Some(aa) = graph.data().aa_node_data.as_ref() {
    sequences.extend(aa.root_aa_sequences.clone());
  }
  Ok(sequences)
}

fn ancestral_node_sequences(
  graph: &GraphAncestral<AncestralGraphData>,
  node_key: GraphNodeKey,
) -> BTreeMap<String, String> {
  let mut sequences = graph
    .data()
    .partition
    .as_ref()
    .map(|partition| {
      let sequence = match partition {
        AncestralPartition::Fitch(partition) => partition.read_arc().node_sequence(node_key),
        AncestralPartition::Sparse(partition) => partition.read_arc().node_sequence(node_key),
        AncestralPartition::Dense(partition) => partition.read_arc().node_sequence(node_key),
      };
      btreemap! { NUC_TRACK.to_owned() => sequence.to_string() }
    })
    .unwrap_or_default();
  if graph.is_root(node_key)
    && let Some(aa) = graph.data().aa_node_data.as_ref()
  {
    sequences.extend(aa.root_aa_sequences.clone());
  }
  sequences
}

fn ancestral_genome_annotations(
  graph: &GraphAncestral<AncestralGraphData>,
  root_sequences: &BTreeMap<String, String>,
) -> Result<Option<AuspiceGenomeAnnotations>, Report> {
  let nuc = root_sequences
    .get(NUC_TRACK)
    .map(|sequence| -> Result<_, Report> {
      Ok(AuspiceGenomeAnnotationNuc {
        start: 1,
        end: isize::try_from(sequence.len()).wrap_err("Nucleotide sequence length does not fit Auspice coordinates")?,
        strand: Some("+".to_owned()),
        r#type: Some("source".to_owned()),
        other: Value::default(),
      })
    })
    .transpose()?;
  let cdses = graph
    .data()
    .aa_node_data
    .as_ref()
    .map(|data| {
      data
        .annotations
        .iter()
        .map(|(name, annotation)| Ok((name.clone(), auspice_cds_annotation(name, annotation)?)))
        .collect::<Result<BTreeMap<_, _>, Report>>()
    })
    .transpose()?
    .unwrap_or_default();
  if nuc.is_none() && cdses.is_empty() {
    return Ok(None);
  }
  Ok(Some(AuspiceGenomeAnnotations {
    nuc,
    cdses,
    other: Value::default(),
  }))
}

fn auspice_cds_annotation(
  name: &str,
  annotation: &AugurNodeDataJsonAnnotationEntry,
) -> Result<AuspiceGenomeAnnotationCds, Report> {
  let strand = annotation
    .strand
    .clone()
    .ok_or_else(|| make_report!("CDS annotation '{name}' has no strand for Auspice output"))?;
  let other = Value::Object(annotation.other.clone().into_iter().collect());
  let segments = if let Some(segments) = annotation.segments.as_ref() {
    if segments.is_empty() {
      return make_error!("CDS annotation '{name}' has no segments for Auspice output");
    }
    Segments::MultipleSegments {
      segments: segments
        .iter()
        .map(|segment| {
          Ok(StartEnd {
            start: auspice_coordinate(segment.start, name, "segment start")?,
            end: auspice_coordinate(segment.end, name, "segment end")?,
            other: Value::Object(segment.other.clone().into_iter().collect()),
          })
        })
        .collect::<Result<Vec<_>, Report>>()?,
      other,
    }
  } else {
    Segments::OneSegment(StartEnd {
      start: auspice_coordinate(
        annotation
          .start
          .ok_or_else(|| make_report!("CDS annotation '{name}' has no start for Auspice output"))?,
        name,
        "start",
      )?,
      end: auspice_coordinate(
        annotation
          .end
          .ok_or_else(|| make_report!("CDS annotation '{name}' has no end for Auspice output"))?,
        name,
        "end",
      )?,
      other,
    })
  };
  Ok(AuspiceGenomeAnnotationCds {
    r#type: annotation.entry_type.clone(),
    gene: None,
    color: None,
    display_name: None,
    description: None,
    strand: Some(strand),
    segments,
  })
}

fn auspice_coordinate(coordinate: i64, name: &str, field: &str) -> Result<isize, Report> {
  isize::try_from(coordinate)
    .wrap_err_with(|| format!("CDS annotation '{name}' {field} does not fit Auspice coordinates"))
}

fn optimize_root_sequences(graph: &GraphAncestral<OptimizeGraphData>) -> Result<BTreeMap<String, String>, Report> {
  let sequence = if let Some(partition) = graph.data().dense_partitions.first() {
    Some(partition.read_arc().root_sequence(graph)?)
  } else if let Some(partition) = graph.data().sparse_partitions.first() {
    Some(partition.read_arc().root_sequence(graph)?)
  } else {
    None
  };
  Ok(
    sequence
      .map(|sequence| btreemap! { NUC_TRACK.to_owned() => sequence.to_string() })
      .unwrap_or_default(),
  )
}

fn optimize_node_sequences(
  graph: &GraphAncestral<OptimizeGraphData>,
  node_key: GraphNodeKey,
) -> BTreeMap<String, String> {
  let sequence = if let Some(partition) = graph.data().dense_partitions.first() {
    Some(partition.read_arc().node_sequence(node_key))
  } else {
    graph
      .data()
      .sparse_partitions
      .first()
      .map(|partition| partition.read_arc().node_sequence(node_key))
  };
  sequence
    .map(|sequence| btreemap! { NUC_TRACK.to_owned() => sequence.to_string() })
    .unwrap_or_default()
}

fn prune_root_sequences(graph: &GraphAncestral<PruneGraphData>) -> Result<BTreeMap<String, String>, Report> {
  graph.data().partitions.first().map_or_else(
    || Ok(BTreeMap::new()),
    |partition| {
      Ok(btreemap! {
        NUC_TRACK.to_owned() => partition.read_arc().root_sequence(graph)?.to_string(),
      })
    },
  )
}

fn prune_node_sequences(graph: &GraphAncestral<PruneGraphData>, node_key: GraphNodeKey) -> BTreeMap<String, String> {
  graph.data().partitions.first().map_or_else(BTreeMap::new, |partition| {
    btreemap! {
      NUC_TRACK.to_owned() => partition.read_arc().node_sequence(node_key).to_string(),
    }
  })
}

fn timetree_root_sequences(
  graph: &Graph<NodeTimetree, EdgeTimetree, TimetreeGraphData>,
) -> Result<BTreeMap<String, String>, Report> {
  graph.data().partitions.first().map_or_else(
    || Ok(BTreeMap::new()),
    |partition| {
      Ok(btreemap! {
        NUC_TRACK.to_owned() => partition.read_arc().root_sequence(graph)?.to_string(),
      })
    },
  )
}

fn timetree_node_sequences(
  graph: &Graph<NodeTimetree, EdgeTimetree, TimetreeGraphData>,
  node_key: GraphNodeKey,
) -> BTreeMap<String, String> {
  graph.data().partitions.first().map_or_else(BTreeMap::new, |partition| {
    btreemap! {
      NUC_TRACK.to_owned() => partition.read_arc().node_sequence(node_key).to_string(),
    }
  })
}

fn ancestral_edge_mutations(
  graph: &GraphAncestral<AncestralGraphData>,
  edge_key: GraphEdgeKey,
) -> Result<Vec<Mutation>, Report> {
  graph.data().partition.as_ref().map_or_else(
    || Ok(vec![]),
    |partition| match partition {
      AncestralPartition::Fitch(partition) => {
        partition
          .read_arc()
          .edge_mutations(graph, edge_key, MutationTrack::Nucleotide)
      },
      AncestralPartition::Sparse(partition) => {
        partition
          .read_arc()
          .edge_mutations(graph, edge_key, MutationTrack::Nucleotide)
      },
      AncestralPartition::Dense(partition) => {
        partition
          .read_arc()
          .edge_mutations(graph, edge_key, MutationTrack::Nucleotide)
      },
    },
  )
}

fn ancestral_node_mutations(
  graph: &GraphAncestral<AncestralGraphData>,
  node_key: GraphNodeKey,
  edge_key: Option<GraphEdgeKey>,
) -> Result<Vec<Mutation>, Report> {
  let mut mutations = edge_key
    .map(|edge_key| ancestral_edge_mutations(graph, edge_key))
    .transpose()?
    .unwrap_or_default();
  if let Some(aa) = graph.data().aa_node_data.as_ref()
    && let Some(tracks) = aa.node_aa_mutations.get(&node_key)
  {
    mutations.extend(tracks.iter().flat_map(|(track, events)| {
      events.iter().cloned().map(|event| Mutation {
        track: MutationTrack::AminoAcid(track.clone()),
        event,
      })
    }));
  }
  Ok(mutations)
}

fn ancestral_all_mutations(graph: &GraphAncestral<AncestralGraphData>) -> Result<Vec<Mutation>, Report> {
  graph
    .get_nodes()
    .into_iter()
    .map(|node| {
      let node = node.read_arc();
      ancestral_node_mutations(graph, node.key(), node.inbound().first().copied())
    })
    .collect::<Result<Vec<_>, _>>()
    .map(|mutations| mutations.into_iter().flatten().collect())
}

fn optimize_has_mutations(graph: &GraphAncestral<OptimizeGraphData>) -> bool {
  !graph.data().dense_partitions.is_empty() || !graph.data().sparse_partitions.is_empty()
}

fn optimize_mutations(
  graph: &GraphAncestral<OptimizeGraphData>,
  edge_key: Option<GraphEdgeKey>,
) -> Result<Vec<Mutation>, Report> {
  let Some(edge_key) = edge_key else {
    return Ok(vec![]);
  };
  if let Some(partition) = graph.data().dense_partitions.first() {
    partition
      .read_arc()
      .edge_mutations(graph, edge_key, MutationTrack::Nucleotide)
  } else if let Some(partition) = graph.data().sparse_partitions.first() {
    partition
      .read_arc()
      .edge_mutations(graph, edge_key, MutationTrack::Nucleotide)
  } else {
    Ok(vec![])
  }
}

fn prune_mutations(
  graph: &GraphAncestral<PruneGraphData>,
  edge_key: Option<GraphEdgeKey>,
) -> Result<Vec<Mutation>, Report> {
  match (graph.data().partitions.first(), edge_key) {
    (Some(partition), Some(edge_key)) => {
      partition
        .read_arc()
        .edge_mutations(graph, edge_key, MutationTrack::Nucleotide)
    },
    _ => Ok(vec![]),
  }
}

fn timetree_mutations(
  graph: &Graph<NodeTimetree, EdgeTimetree, TimetreeGraphData>,
  edge_key: Option<GraphEdgeKey>,
) -> Result<Vec<Mutation>, Report> {
  match (graph.data().partitions.first(), edge_key) {
    (Some(partition), Some(edge_key)) => {
      partition
        .read_arc()
        .edge_mutations(graph, edge_key, MutationTrack::Nucleotide)
    },
    _ => Ok(vec![]),
  }
}

fn mugration_traits(
  graph: &GraphAncestral<MugrationGraphData>,
  node_key: GraphNodeKey,
  node_name: &str,
) -> Result<BTreeMap<String, TraitValue>, Report> {
  let Some(value) = graph.data().partition.get_reconstructed_trait(node_key) else {
    return Ok(BTreeMap::new());
  };
  let profile = graph.data().partition.get_confidence(node_key);
  if let Some(profile) = profile.as_ref() {
    for (state, probability) in graph.data().partition.states.iter().zip(profile) {
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
    .map(|profile| build_confidence_map(&graph.data().partition.states, profile))
    .unwrap_or_default();
  let entropy = profile.as_ref().map(compute_entropy);
  if let Some(entropy) = entropy {
    ensure_finite(entropy, "mugration", node_name, "trait entropy")?;
  }
  Ok(btreemap! {
    graph.data().traits.attribute.clone() => TraitValue { value, confidence, entropy },
  })
}

fn mugration_transition(
  graph: &GraphAncestral<MugrationGraphData>,
  node_key: GraphNodeKey,
) -> Result<Option<(String, String)>, Report> {
  let Some((parent_key, _edge_key)) = graph.node_parent(node_key)? else {
    return Ok(None);
  };
  let parent = graph.data().partition.get_reconstructed_trait(parent_key);
  let child = graph.data().partition.get_reconstructed_trait(node_key);
  Ok(match (parent, child) {
    (Some(parent), Some(child)) if parent != child => Some((parent, child)),
    _ => None,
  })
}

fn mugration_transition_label(
  graph: &GraphAncestral<MugrationGraphData>,
  node_key: GraphNodeKey,
) -> Result<Option<AuspiceTreeBranchAttrsLabels>, Report> {
  let Some((parent, child)) = mugration_transition(graph, node_key)? else {
    return Ok(None);
  };
  Ok(Some(AuspiceTreeBranchAttrsLabels {
    aa: None,
    clade: None,
    other: json!({ graph.data().traits.attribute.clone(): format!("{parent} → {child}") }),
  }))
}

fn timetree_divergence(
  graph: &Graph<NodeTimetree, EdgeTimetree, TimetreeGraphData>,
  node_key: GraphNodeKey,
  node: &NodeTimetree,
) -> Result<f64, Report> {
  graph.data().mutation_counts.as_ref().map_or(Ok(node.div), |counts| {
    let mut key = node_key;
    let mut count = 0;
    while let Some((parent, edge)) = graph.node_parent(key)? {
      count += counts.get(&edge).copied().unwrap_or_default();
      key = parent;
    }
    Ok(count as f64)
  })
}

fn timetree_date_confidence(
  graph: &Graph<NodeTimetree, EdgeTimetree, TimetreeGraphData>,
  node_key: GraphNodeKey,
  node_name: &str,
) -> Result<Option<[f64; 2]>, Report> {
  let confidence = graph
    .data()
    .confidence_intervals
    .as_ref()
    .and_then(|intervals| intervals.iter().find(|interval| interval.key == node_key))
    .map(|interval| [interval.lower, interval.upper]);
  if let Some([lower, upper]) = confidence {
    ensure_finite(lower, "timetree", node_name, "date confidence lower bound")?;
    ensure_finite(upper, "timetree", node_name, "date confidence upper bound")?;
    Ok(Some([format_number(lower, 3), format_number(upper, 3)]))
  } else {
    Ok(None)
  }
}

fn timetree_date_is_inferred(
  graph: &Graph<NodeTimetree, EdgeTimetree, TimetreeGraphData>,
  node_key: GraphNodeKey,
  node: &NodeTimetree,
) -> Option<bool> {
  let dates = graph.data().dates.as_ref()?;
  let name = node.base.name.as_ref();
  Some(
    name.and_then(|name| dates.get(name)).and_then(Option::as_ref).is_none()
      && node.time.is_some()
      && graph.get_node(node_key).is_some(),
  )
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
    // The nucleotide auspice mutation list mirrors the augur node-data `muts`,
    // which are substitution-only: augur export copies node-data muts verbatim
    // into `branch_attrs.mutations.nuc`, so indels (a separate track) must not
    // appear there. Amino-acid tracks keep indels to match their own node-data.
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

fn mutation_property(mutation: &Mutation) -> Result<PhyloxmlProperty, Report> {
  let track = match &mutation.track {
    MutationTrack::Nucleotide => NUC_TRACK.to_owned(),
    MutationTrack::AminoAcid(track) => format!("aa:{}", encode_property_token(track)),
  };
  let value = match &mutation.event {
    MutationEvent::Substitution(substitution) => {
      let position = substitution
        .pos()
        .checked_add(1)
        .ok_or_else(|| eyre::eyre!("Mutation coordinate overflow at {}", substitution.pos()))?;
      format!("{track}:sub:{}{position}{}", substitution.reff(), substitution.qry())
    },
    MutationEvent::Insertion(segment) => {
      let (start, end) = segment.one_based_inclusive_range()?;
      format!("{track}:ins:{start}-{end}:{}", segment.sequence)
    },
    MutationEvent::Deletion(segment) => {
      let (start, end) = segment.one_based_inclusive_range()?;
      format!("{track}:del:{start}-{end}:{}", segment.sequence)
    },
  };
  Ok(property(REF_MUTATION, DT_STRING, APPLIES_BRANCH, &value))
}

fn trait_properties(traits: &BTreeMap<String, TraitValue>, node_name: &str) -> Result<Vec<PhyloxmlProperty>, Report> {
  let mut properties = vec![];
  for (attribute, value) in traits {
    let attribute = encode_property_token(attribute);
    properties.push(property(
      &format!("{REF_TRAIT_PREFIX}{attribute}"),
      DT_STRING,
      APPLIES_NODE,
      &value.value,
    ));
    for (state, probability) in &value.confidence {
      ensure_finite(
        *probability,
        "mugration",
        node_name,
        &format!("trait confidence '{state}'"),
      )?;
      properties.push(property(
        &format!(
          "{REF_TRAIT_CONFIDENCE_PREFIX}{attribute}:{}",
          encode_property_token(state)
        ),
        DT_DOUBLE,
        APPLIES_NODE,
        &probability.to_string(),
      ));
    }
    if let Some(entropy) = value.entropy {
      ensure_finite(entropy, "mugration", node_name, "trait entropy")?;
      properties.push(property(
        &format!("{REF_TRAIT_ENTROPY_PREFIX}{attribute}"),
        DT_DOUBLE,
        APPLIES_NODE,
        &entropy.to_string(),
      ));
    }
  }
  Ok(properties)
}

fn empty_phyloxml_clade(name: Option<String>, branch_length: Option<f64>) -> PhyloxmlClade {
  PhyloxmlClade {
    name,
    branch_length_elem: branch_length,
    branch_length_attr: None,
    confidence: vec![],
    width: None,
    color: None,
    node_id: None,
    taxonomy: vec![],
    sequence: vec![],
    events: None,
    binary_characters: None,
    distribution: vec![],
    date: None,
    reference: vec![],
    property: vec![],
    clade: vec![],
    other: BTreeMap::new(),
  }
}

fn input_branch_confidence(
  confidence: Option<f64>,
  command: &str,
  node_name: &str,
) -> Result<Vec<PhyloxmlConfidence>, Report> {
  confidence.map_or_else(
    || Ok(vec![]),
    |value| {
      ensure_finite(value, command, node_name, "input branch support")?;
      Ok(vec![PhyloxmlConfidence {
        value,
        type_: TYPE_INPUT_BRANCH_SUPPORT.to_owned(),
      }])
    },
  )
}

fn phyloxml_sequences(sequences: &BTreeMap<String, String>) -> Vec<PhyloxmlSequence> {
  sequences
    .iter()
    .map(|(track, sequence)| PhyloxmlSequence {
      symbol: None,
      accession: None,
      name: Some(track.clone()),
      location: None,
      mol_seq: Some(PhyloxmlMolSeq {
        sequence: sequence.clone(),
        is_aligned: Some(true),
      }),
      uri: None,
      annotation: vec![],
      domain_architecture: None,
      other: BTreeMap::new(),
    })
    .collect()
}

fn phyloxml_date(value: f64) -> PhyloxmlDate {
  PhyloxmlDate {
    desc: None,
    value: Some(value),
    minimum: None,
    maximum: None,
    unit: Some("year".to_owned()),
  }
}

fn property(ref_: &str, datatype: &str, applies_to: &str, value: &str) -> PhyloxmlProperty {
  PhyloxmlProperty {
    value: value.to_owned(),
    ref_: ref_.to_owned(),
    unit: None,
    datatype: datatype.to_owned(),
    applies_to: applies_to.to_owned(),
    id_ref: None,
  }
}

fn encode_property_token(value: &str) -> String {
  utf8_percent_encode(value, PROPERTY_TOKEN_ENCODE_SET).to_string()
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

fn cumulative_branch_length<N, E, D>(graph: &Graph<N, E, D>, mut key: GraphNodeKey) -> Result<Option<f64>, Report>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
{
  let mut total = 0.0;
  while let Some((parent, edge_key)) = graph.node_parent(key)? {
    let edge = graph
      .get_edge(edge_key)
      .ok_or_else(|| make_internal_report!("Edge {edge_key} disappeared while summing branch lengths"))?;
    let Some(length) = edge.read_arc().payload().read_arc().branch_length() else {
      return Ok(None);
    };
    total += length;
    key = parent;
  }
  Ok(Some(total))
}

fn node_name<N: Named>(key: GraphNodeKey, node: &N) -> String {
  node
    .name()
    .map_or_else(|| format!("node_{}", key.as_usize()), |name| name.as_ref().to_owned())
}

fn ensure_optional_finite(value: Option<f64>, command: &str, node_name: &str, field: &str) -> Result<(), Report> {
  if let Some(value) = value {
    ensure_finite(value, command, node_name, field)?;
  }
  Ok(())
}

fn ensure_finite(value: f64, command: &str, node_name: &str, field: &str) -> Result<(), Report> {
  if value.is_finite() {
    Ok(())
  } else {
    make_error!("{command} node '{node_name}' has non-finite {field}={value}")
  }
}

fn finite_number(
  value: Option<f64>,
  precision: i32,
  command: &str,
  node_name: &str,
  field: &str,
) -> Result<Option<f64>, Report> {
  ensure_optional_finite(value, command, node_name, field)?;
  Ok(value.map(|value| format_number(value, precision)))
}

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

fn coloring(key: &str, title: &str, type_: &str) -> AuspiceColoring {
  AuspiceColoring {
    key: key.to_owned(),
    title: title.to_owned(),
    type_: type_.to_owned(),
    ..AuspiceColoring::default()
  }
}

fn generation_date() -> String {
  Utc::now().format("%Y-%m-%d").to_string()
}
