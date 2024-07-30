use crate::graph::edge::{Edge, GraphEdge};
use crate::graph::graph::Graph;
use crate::graph::node::GraphNode;
use crate::io::file::create_file_or_stdout;
use crate::io::file::open_file_or_stdin;
use crate::io::json::{json_read, json_write, JsonPretty};
use eyre::{Report, WrapErr};
use maplit::{btreemap, btreeset};
use parking_lot::RwLock;
use std::collections::VecDeque;
use std::io::Cursor;
use std::io::{Read, Write};
use std::path::Path;
use std::sync::Arc;

pub fn auspice_read_file<C, N, E, D>(filepath: impl AsRef<Path>) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
  C: AuspiceRead<N, E, D>,
{
  let filepath = filepath.as_ref();
  auspice_read::<C, _, _, _>(open_file_or_stdin(&Some(filepath))?)
    .wrap_err_with(|| format!("When reading Auspice v2 JSON file '{filepath:#?}'"))
}

pub fn auspice_read_str<C, N, E, D>(auspice_string: impl AsRef<str>) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
  C: AuspiceRead<N, E, D>,
{
  let auspice_string = auspice_string.as_ref();
  auspice_read::<C, _, _, _>(Cursor::new(auspice_string)).wrap_err("When reading Auspice v2 JSON string")
}

pub fn auspice_read<C, N, E, D>(reader: impl Read) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
  C: AuspiceRead<N, E, D>,
{
  let tree: AuspiceTree = json_read(reader).wrap_err("When reading Auspice v2 JSON")?;
  auspice_to_graph::<C, _, _, _>(&tree).wrap_err("When converting Auspice v2 JSON to graph")
}

pub fn auspice_write_file<C, N, E, D>(filepath: impl AsRef<Path>, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
  C: AuspiceWrite<N, E, D>,
{
  let filepath = filepath.as_ref();
  let mut f = create_file_or_stdout(filepath)?;
  auspice_write::<C, _, _, _>(&mut f, graph)
    .wrap_err_with(|| format!("When reading Auspice v2 JSON file '{filepath:#?}'"))?;
  writeln!(f)?;
  Ok(())
}

pub fn auspice_write_str<C, N, E, D>(graph: &Graph<N, E, D>) -> Result<String, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
  C: AuspiceWrite<N, E, D>,
{
  let mut buf = Vec::new();
  auspice_write::<C, _, _, _>(&mut buf, graph).wrap_err("When writing Auspice v2 JSON string")?;
  Ok(String::from_utf8(buf)?)
}

pub fn auspice_write<C, N, E, D>(writer: &mut impl Write, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
  C: AuspiceWrite<N, E, D>,
{
  let tree = auspice_from_graph::<C, _, _, _>(graph)?;
  json_write(writer, &tree, JsonPretty(true)).wrap_err("When writing Auspice v2 JSON")
}

pub struct AuspiceTreeContext<'a> {
  pub node: &'a AuspiceTreeNode,
  pub tree: &'a AuspiceTree,
}

pub trait AuspiceRead<N, E, D>: Sized
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  fn new(tree: &AuspiceTree) -> Result<Self, Report>;

  fn auspice_data_to_graph_data(&mut self, tree: &AuspiceTree) -> Result<D, Report>;

  fn auspice_node_to_graph_components(&mut self, context: &AuspiceTreeContext) -> Result<(N, E), Report>;
}

/// Convert Auspice v2 JSON to graph
pub fn auspice_to_graph<C, N, E, D>(tree: &AuspiceTree) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
  C: AuspiceRead<N, E, D>,
{
  let mut converter = C::new(tree)?;
  let mut graph = Graph::<N, E, D>::with_data(converter.auspice_data_to_graph_data(tree)?);
  let mut queue = VecDeque::from([(None, &tree.tree)]);
  while let Some((parent_key, node)) = queue.pop_front() {
    let (graph_node, graph_edge) = converter.auspice_node_to_graph_components(&AuspiceTreeContext { node, tree })?;
    let node_key = graph.add_node(graph_node);
    if let Some(parent_key) = parent_key {
      graph.add_edge(parent_key, node_key, graph_edge)?;
    }
    for child in &node.children {
      queue.push_back((Some(node_key), &child));
    }
  }
  graph.build()?;
  Ok(graph)
}

pub struct AuspiceGraphContext<'a, N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  pub node: &'a N,
  pub edge: Option<&'a E>,
  pub graph: &'a Graph<N, E, D>,
}

pub trait AuspiceWrite<N, E, D>: Sized
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  fn new(graph: &Graph<N, E, D>) -> Result<Self, Report>;

  fn auspice_data_from_graph_data(&self, graph: &Graph<N, E, D>) -> Result<AuspiceTreeData, Report>;

  fn auspice_node_from_graph_components(
    &mut self,
    context: &AuspiceGraphContext<N, E, D>,
  ) -> Result<AuspiceTreeNode, Report>;
}

/// Convert graph to Auspice v2 JSON
pub fn auspice_from_graph<C, N, E, D>(graph: &Graph<N, E, D>) -> Result<AuspiceTree, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
  C: AuspiceWrite<N, E, D>,
{
  let mut converter = C::new(graph)?;
  let root = graph.get_exactly_one_root().wrap_err("When writing Auspice v2 JSON")?;

  // Pre-order iteration to construct the nodes from graph nodes and edges
  let mut node_map = {
    let mut node_map = btreemap! {};
    let mut queue = VecDeque::from([(Arc::clone(&root), None)]);
    while let Some((current_node, current_edge)) = queue.pop_front() {
      let current_node = current_node.read_arc();

      let node = &*current_node.payload().read_arc();
      let edge = current_edge
        .as_ref()
        .map(|edge: &Arc<RwLock<Edge<E>>>| edge.read_arc().payload().read_arc());
      let edge = edge.as_deref();
      let current_tree_node =
        converter.auspice_node_from_graph_components(&AuspiceGraphContext { node, edge, graph })?;
      for (child, edge) in graph.children_of(&current_node) {
        queue.push_back((child, Some(edge)));
      }
      node_map.insert(current_node.key(), current_tree_node);
    }
    node_map
  };

  // Post-order traversal to populate .children array
  let mut visited = btreeset! {};
  let mut stack = vec![Arc::clone(&root)];
  while let Some(node) = stack.pop() {
    if visited.contains(&node.read_arc().key()) {
      let mut children_to_add = vec![];
      for (child, _) in graph.children_of(&node.read_arc()) {
        let child_key = child.read_arc().key();
        if let Some(child_node) = node_map.remove(&child_key) {
          children_to_add.push(child_node);
        }
      }
      if let Some(node) = node_map.get_mut(&node.read_arc().key()) {
        node.children.extend(children_to_add);
      }
    } else {
      visited.insert(node.read_arc().key());
      stack.push(Arc::clone(&node));
      for (child, _) in graph.children_of(&node.read_arc()) {
        stack.push(child);
      }
    }
  }

  let data = converter.auspice_data_from_graph_data(graph)?;
  let tree = node_map.remove(&root.read_arc().key()).unwrap();
  Ok(AuspiceTree { data, tree })
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::graph::graph::tests::{TestEdge, TestNode};
  use crate::io::auspice::tests::AuspiceTreeBranchAttrs;
  use crate::io::auspice::{AuspiceTreeMeta, AuspiceTreeNodeAttrs};
  use crate::make_internal_report;
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use serde_json::Value;
  use std::collections::BTreeMap;

  #[derive(Clone, Debug)]
  pub struct TestAuspiceData {
    pub version: Option<String>,
    pub meta: AuspiceTreeMeta,
    pub root_sequence: Option<BTreeMap<String, String>>,
    pub other: Value,
  }

  pub struct AuspiceWriter {}

  impl AuspiceWrite<TestNode, TestEdge, TestAuspiceData> for AuspiceWriter {
    fn new(graph: &Graph<TestNode, TestEdge, TestAuspiceData>) -> Result<Self, Report> {
      Ok(Self {})
    }

    fn auspice_data_from_graph_data(
      &self,
      graph: &Graph<TestNode, TestEdge, TestAuspiceData>,
    ) -> Result<AuspiceTreeData, Report> {
      let data = graph.data().read_arc();
      Ok(AuspiceTreeData {
        version: data.version.clone(),
        meta: data.meta.clone(),
        root_sequence: data.root_sequence.clone(),
        other: data.other.clone(),
      })
    }

    fn auspice_node_from_graph_components(
      &mut self,
      context: &AuspiceGraphContext<TestNode, TestEdge, TestAuspiceData>,
    ) -> Result<AuspiceTreeNode, Report> {
      let AuspiceGraphContext { node, edge, graph } = context;

      let name = node
        .0
        .as_ref()
        .ok_or_else(|| make_internal_report!("Encountered node with empty name"))?
        .to_owned();

      Ok(AuspiceTreeNode {
        name,
        branch_attrs: AuspiceTreeBranchAttrs::default(),
        node_attrs: AuspiceTreeNodeAttrs {
          div: edge.and_then(|edge| edge.0),
          clade_membership: None,
          region: None,
          country: None,
          division: None,
          other: Value::default(),
        },
        children: vec![],
        other: Value::default(),
      })
    }
  }

  pub struct AuspiceReader {}

  impl AuspiceRead<TestNode, TestEdge, TestAuspiceData> for AuspiceReader {
    fn new(tree: &AuspiceTree) -> Result<Self, Report> {
      Ok(Self {})
    }

    fn auspice_data_to_graph_data(&mut self, tree: &AuspiceTree) -> Result<TestAuspiceData, Report> {
      Ok(TestAuspiceData {
        version: tree.data.version.clone(),
        meta: tree.data.meta.clone(),
        root_sequence: tree.data.root_sequence.clone(),
        other: tree.data.other.clone(),
      })
    }

    fn auspice_node_to_graph_components(
      &mut self,
      context: &AuspiceTreeContext,
    ) -> Result<(TestNode, TestEdge), Report> {
      let AuspiceTreeContext { node, .. } = context;
      Ok((TestNode(Some(node.name.clone())), TestEdge(node.node_attrs.div)))
    }
  }

  #[test]
  fn test_auspice_roundtrip() -> Result<(), Report> {
    let input = indoc!(
      // language=json
      r#"{
        "version": "v2",
        "meta": {
          "description": "This is a test!",
          "name": "Test"
        },
        "root_sequence": {
          "nuc": "ACGTACGTACGTACGTACGTACGT"
        },
        "unknown_field": 42,
        "tree": {
          "name": "root",
          "node_attrs": {},
          "children": [
            {
              "name": "AB",
              "node_attrs": {
                "div": 3.0
              },
              "children": [
                {
                  "name": "A",
                  "node_attrs": {
                    "div": 8.0
                  }
                },
                {
                  "name": "B",
                  "node_attrs": {
                    "div": 5.0
                  }
                }
              ]
            },
            {
              "name": "C",
              "node_attrs": {
                "div": 7.0
              }
            }
          ]
        }
      }"#
    );

    let graph = auspice_read_str::<AuspiceReader, _, _, _>(input)?;
    let output = auspice_write_str::<AuspiceWriter, _, _, _>(&graph)?;
    assert_eq!(input, output);
    Ok(())
  }
}

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::slice::Iter;
use traversal::{Bft, DftPost, DftPre};

#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct AuspiceGraphMeta {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub auspice_tree_version: Option<String>,

  pub meta: AuspiceTreeMeta,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct AuspiceTreeNodeAttr {
  pub value: String,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

impl AuspiceTreeNodeAttr {
  pub fn new(value: &str) -> Self {
    Self {
      value: value.to_owned(),
      other: serde_json::Value::default(),
    }
  }
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct AuspiceTreeNodeAttrF64 {
  pub value: f64,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

impl AuspiceTreeNodeAttrF64 {
  pub fn new(value: f64) -> Self {
    Self {
      value,
      other: serde_json::Value::default(),
    }
  }
}

#[derive(Clone, Default, Eq, PartialEq, Serialize, Deserialize, Debug)]
pub struct AuspiceTreeBranchAttrsLabels {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub aa: Option<String>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub clade: Option<String>,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

#[derive(Clone, Default, Eq, PartialEq, Serialize, Deserialize, Debug)]
pub struct AuspiceTreeBranchAttrs {
  pub mutations: BTreeMap<String, Vec<String>>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub labels: Option<AuspiceTreeBranchAttrsLabels>,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

impl AuspiceTreeBranchAttrs {
  #[inline]
  pub fn is_default(&self) -> bool {
    self == &Self::default()
  }
}

#[derive(Clone, Default, Serialize, Deserialize, Debug)]
pub struct AuspiceTreeNodeAttrs {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub div: Option<f64>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub clade_membership: Option<AuspiceTreeNodeAttr>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub region: Option<AuspiceTreeNodeAttr>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub country: Option<AuspiceTreeNodeAttr>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub division: Option<AuspiceTreeNodeAttr>,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct AuspiceTreeNode {
  pub name: String,

  #[serde(default, skip_serializing_if = "AuspiceTreeBranchAttrs::is_default")]
  pub branch_attrs: AuspiceTreeBranchAttrs,

  pub node_attrs: AuspiceTreeNodeAttrs,

  #[serde(skip_serializing_if = "Vec::is_empty")]
  #[serde(default)]
  pub children: Vec<AuspiceTreeNode>,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct AuspiceColoring {
  #[serde(rename = "type")]
  pub type_: String,

  pub key: String,

  pub title: String,

  #[serde(skip_serializing_if = "Vec::is_empty")]
  #[serde(default)]
  pub scale: Vec<[String; 2]>,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

#[derive(Clone, Default, Serialize, Deserialize, Eq, PartialEq, Debug)]
pub struct AuspiceDisplayDefaults {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub branch_label: Option<String>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub color_by: Option<String>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub distance_measure: Option<String>,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

impl AuspiceDisplayDefaults {
  pub fn is_empty(&self) -> bool {
    self == &Self::default()
  }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AuspiceGenomeAnnotationNuc {
  pub start: isize,

  pub end: isize,

  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub strand: Option<String>,

  #[serde(rename = "type", default, skip_serializing_if = "Option::is_none")]
  pub r#type: Option<String>,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct StartEnd {
  pub start: isize,
  pub end: isize,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(untagged)]
pub enum Segments {
  OneSegment(StartEnd),
  MultipleSegments {
    segments: Vec<StartEnd>,

    #[serde(flatten)]
    other: serde_json::Value,
  },
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AuspiceGenomeAnnotationCds {
  #[serde(rename = "type", default, skip_serializing_if = "Option::is_none")]
  pub r#type: Option<String>,

  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub gene: Option<String>,

  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub color: Option<String>,

  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub display_name: Option<String>,

  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub description: Option<String>,

  #[serde(rename = "type", default, skip_serializing_if = "Option::is_none")]
  pub strand: Option<String>,

  #[serde(flatten)]
  pub segments: Segments,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AuspiceGenomeAnnotations {
  pub nuc: AuspiceGenomeAnnotationNuc,

  #[serde(flatten)]
  pub cdses: BTreeMap<String, AuspiceGenomeAnnotationCds>,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

#[derive(Clone, Default, Serialize, Deserialize, Debug)]
pub struct AuspiceTreeMeta {
  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub title: Option<String>,

  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub description: Option<String>,

  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub updated: Option<String>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub genome_annotations: Option<AuspiceGenomeAnnotations>,

  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub colorings: Vec<AuspiceColoring>,

  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub panels: Vec<String>,

  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub filters: Vec<String>,

  #[serde(default, skip_serializing_if = "AuspiceDisplayDefaults::is_empty")]
  pub display_defaults: AuspiceDisplayDefaults,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub geo_resolutions: Option<serde_json::Value>,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

#[repr(u8)]
#[derive(Debug, Clone, Copy, Eq, PartialEq, Default)]
pub enum DivergenceUnits {
  NumSubstitutionsPerYearPerSite,
  #[default]
  NumSubstitutionsPerYear,
}

impl DivergenceUnits {
  ///
  /// Guesses the unit of measurement of divergence, based on the greatest value of divergence on the tree
  ///
  pub fn guess_from_max_divergence(max_divergence: f64) -> DivergenceUnits {
    // FIXME: This should be fixed upstream in augur & auspice, but it is hard to do without breaking Auspice JSON v2 format.
    // Taken from: https://github.com/nextstrain/auspice/blob/6a2d0f276fccf05bfc7084608bb0010a79086c83/src/components/tree/phyloTree/renderers.js#L376
    // A quote from there:
    //  > Prior to Jan 2020, the divergence measure was always "subs per site per year"
    //  > however certain datasets changed this to "subs per year" across entire sequence.
    //  > This distinction is not set in the JSON, so in order to correctly display the rate
    //  > we will "guess" this here. A future augur update will export this in a JSON key,
    //  > removing the need to guess
    //
    // HACK: Arbitrary threshold to make a guess
    const HACK_MAX_DIVERGENCE_THRESHOLD: f64 = 5.0;
    if max_divergence <= HACK_MAX_DIVERGENCE_THRESHOLD {
      DivergenceUnits::NumSubstitutionsPerYearPerSite
    } else {
      DivergenceUnits::NumSubstitutionsPerYear
    }
  }
}

#[derive(Clone, Default, Serialize, Deserialize, Debug)]
pub struct AuspiceTreeData {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub version: Option<String>,

  pub meta: AuspiceTreeMeta,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub root_sequence: Option<BTreeMap<String, String>>,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct AuspiceTree {
  #[serde(flatten)]
  pub data: AuspiceTreeData,

  pub tree: AuspiceTreeNode,
}

pub type AuspiceTreeNodeIter<'a> = Iter<'a, AuspiceTreeNode>;

pub type AuspiceTreeNodeIterFn<'a> = fn(&'a AuspiceTreeNode) -> AuspiceTreeNodeIter<'_>;

impl AuspiceTree {
  /// Returns iterator for breadth-first tree traversal
  pub fn iter_breadth_first<'a>(
    &'a self,
  ) -> Bft<'a, AuspiceTreeNode, AuspiceTreeNodeIterFn<'a>, AuspiceTreeNodeIter<'a>> {
    Bft::new(&self.tree, |node: &'a AuspiceTreeNode| node.children.iter())
  }

  /// Returns iterator for depth-first pre-order tree traversal
  pub fn iter_depth_first_preorder<'a>(
    &'a self,
  ) -> DftPre<'a, AuspiceTreeNode, AuspiceTreeNodeIterFn<'a>, AuspiceTreeNodeIter<'a>> {
    DftPre::new(&self.tree, |node: &'a AuspiceTreeNode| node.children.iter())
  }

  /// Returns iterator for depth-first post-order tree traversal
  pub fn iter_depth_first_postorder<'a>(
    &'a self,
  ) -> DftPost<'a, AuspiceTreeNode, AuspiceTreeNodeIterFn<'a>, AuspiceTreeNodeIter<'a>> {
    DftPost::new(&self.tree, |node: &'a AuspiceTreeNode| node.children.iter())
  }

  fn map_nodes_rec(index: usize, node: &AuspiceTreeNode, action: fn((usize, &AuspiceTreeNode))) {
    action((index, node));
    for child in &node.children {
      Self::map_nodes_rec(index + 1, child, action);
    }
  }

  /// Iterates over nodes and applies a function to each. Immutable version.
  pub fn map_nodes(&self, action: fn((usize, &AuspiceTreeNode))) {
    Self::map_nodes_rec(0, &self.tree, action);
  }

  fn map_nodes_mut_rec(index: usize, node: &mut AuspiceTreeNode, action: fn((usize, &mut AuspiceTreeNode))) {
    action((index, node));
    for child in &mut node.children {
      Self::map_nodes_mut_rec(index + 1, child, action);
    }
  }

  /// Iterates over nodes and applies a function to each. Mutable version.
  pub fn map_nodes_mut(&mut self, action: fn((usize, &mut AuspiceTreeNode))) {
    Self::map_nodes_mut_rec(0, &mut self.tree, action);
  }

  pub fn root_sequence(&self) -> Option<&str> {
    self
      .data
      .root_sequence
      .as_ref()
      .and_then(|root_sequence| root_sequence.get("nuc"))
      .map(String::as_str)
  }
}
