use crate::graph::assign_node_names::assign_node_names;
use crate::graph::edge::{Edge, GraphEdge};
use crate::graph::graph::{Graph, GraphElement, GraphElementRef};
use crate::graph::node::{GraphNode, Named};
use crate::io::auspice::auspice::{AuspiceTree, AuspiceTreeMeta, AuspiceTreeNode};
use crate::io::file::create_file_or_stdout;
use crate::io::file::open_file_or_stdin;
use crate::io::json::{json_read, json_write, JsonPretty};
use crate::o;
use eyre::{Report, WrapErr};
use indexmap::{indexmap, indexset};
use parking_lot::RwLock;
use serde_json::Value;
use std::collections::VecDeque;
use std::io::Cursor;
use std::io::{Read, Write};
use std::path::Path;
use std::sync::Arc;

pub fn auspice_read_file<N, E>(filepath: impl AsRef<Path>) -> Result<Graph<N, E>, Report>
where
  N: GraphNode + Named,
  E: GraphEdge,
  GraphElement<N, E>: FromAuspiceNode,
{
  let filepath = filepath.as_ref();
  auspice_read(open_file_or_stdin(&Some(filepath))?)
    .wrap_err_with(|| format!("When reading Auspice v2 JSON file '{filepath:#?}'"))
}

pub fn auspice_read_str<N, E>(auspice_string: impl AsRef<str>) -> Result<Graph<N, E>, Report>
where
  N: GraphNode + Named,
  E: GraphEdge,
  GraphElement<N, E>: FromAuspiceNode,
{
  let auspice_string = auspice_string.as_ref();
  auspice_read(Cursor::new(auspice_string)).wrap_err("When reading Auspice v2 JSON string")
}

pub fn auspice_read<N, E>(reader: impl Read) -> Result<Graph<N, E>, Report>
where
  N: GraphNode + Named,
  E: GraphEdge,
  GraphElement<N, E>: FromAuspiceNode,
{
  let tree: AuspiceTree = json_read(reader).wrap_err("When reading Auspice v2 JSON")?;
  let tree_root = tree.tree;

  let mut queue = VecDeque::new();
  queue.push_back((None, tree_root));

  let mut graph = Graph::<N, E>::new();
  while let Some((parent_key, node)) = queue.pop_front() {
    let ge = GraphElement::from_auspice_node(&node);
    let node_key = graph.add_node(ge.node);
    if let Some(parent_key) = parent_key {
      graph.add_edge(parent_key, node_key, ge.edge)?;
    }
    for child in node.children {
      queue.push_back((Some(node_key), child));
    }
  }

  graph.build()?;
  assign_node_names(&graph);
  Ok(graph)
}

pub fn auspice_write_file<N, E>(filepath: impl AsRef<Path>, graph: &Graph<N, E>) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  for<'a> GraphElementRef<'a, N, E>: ToAuspiceNode,
{
  let filepath = filepath.as_ref();
  let mut f = create_file_or_stdout(filepath)?;
  auspice_write(&mut f, graph).wrap_err_with(|| format!("When reading Auspice v2 JSON file '{filepath:#?}'"))?;
  writeln!(f)?;
  Ok(())
}

pub fn auspice_write_str<N, E>(graph: &Graph<N, E>) -> Result<String, Report>
where
  N: GraphNode,
  E: GraphEdge,
  for<'a> GraphElementRef<'a, N, E>: ToAuspiceNode,
{
  let mut buf = Vec::new();
  auspice_write(&mut buf, graph).wrap_err("When writing Auspice v2 JSON string")?;
  Ok(String::from_utf8(buf)?)
}

pub fn auspice_write<N, E>(writer: &mut impl Write, graph: &Graph<N, E>) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  for<'a> GraphElementRef<'a, N, E>: ToAuspiceNode,
{
  let root = graph.get_exactly_one_root().wrap_err("When writing Auspice v2 JSON")?;

  // Pre-order iteration to construct the nodes from graph nodes and edges
  let mut node_map = {
    let mut node_map = indexmap! {};
    let mut queue = VecDeque::from([(Arc::clone(&root), None)]);
    while let Some((current_node, current_edge)) = queue.pop_front() {
      let current_node = current_node.read_arc();

      let current_node_payload = &*current_node.payload().read_arc();
      let current_edge_payload = current_edge
        .as_ref()
        .map(|edge: &Arc<RwLock<Edge<E>>>| edge.read_arc().payload().read_arc());

      let current_tree_node =
        GraphElementRef::new(current_node_payload, current_edge_payload.as_deref()).to_auspice_node();
      for (child, edge) in graph.children_of(&current_node) {
        queue.push_back((child, Some(edge)));
      }
      node_map.insert(current_node.key(), current_tree_node);
    }
    node_map
  };

  // Post-order traversal to populate .children array
  let mut visited = indexset! {};
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

  let root = node_map.remove(&root.read_arc().key()).unwrap();

  // TODO: fill missing fields
  let tree = AuspiceTree {
    version: Some(o!("v2")),
    meta: AuspiceTreeMeta::default(),
    tree: root,
    root_sequence: None,
    other: Value::default(),
  };

  json_write(writer, &tree, JsonPretty(true)).wrap_err("When writing Auspice v2 JSON")
}

/// Defines conversion to tree node when writing to Auspice v2 JSON
pub trait ToAuspiceNode: Sized {
  fn to_auspice_node(&self) -> AuspiceTreeNode;
}

/// Defines conversion from tree node when reading from Auspice v2 JSON
pub trait FromAuspiceNode: Sized {
  fn from_auspice_node(tree_node: &AuspiceTreeNode) -> Self;
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::graph::graph::tests::{TestEdge, TestNode};
  use crate::io::auspice::auspice::AuspiceTreeNodeAttrs;
  use crate::io::auspice::tests::auspice::AuspiceTreeBranchAttrs;
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;

  impl ToAuspiceNode for GraphElementRef<'_, TestNode, TestEdge> {
    fn to_auspice_node(&self) -> AuspiceTreeNode {
      AuspiceTreeNode {
        name: self.node.0.clone(),
        branch_attrs: AuspiceTreeBranchAttrs::default(),
        node_attrs: AuspiceTreeNodeAttrs {
          div: self.edge.map(|edge| edge.0),
          clade_membership: None,
          region: None,
          country: None,
          division: None,
          other: Value::default(),
        },
        children: vec![],
        other: Value::default(),
      }
    }
  }

  impl FromAuspiceNode for GraphElement<TestNode, TestEdge> {
    fn from_auspice_node(node: &AuspiceTreeNode) -> Self {
      Self {
        node: TestNode(node.name.clone()),
        edge: TestEdge(node.node_attrs.div.unwrap_or_default()),
      }
    }
  }

  #[test]
  fn test_auspice_roundtrip() -> Result<(), Report> {
    let input = indoc!(
      // language=json
      r#"{
        "version": "v2",
        "meta": {},
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

    let graph = auspice_read_str::<TestNode, TestEdge>(input)?;
    let output = auspice_write_str(&graph)?;
    assert_eq!(input, output);
    Ok(())
  }
}

mod auspice {
  #![allow(dead_code)]

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

  #[derive(Clone, Serialize, Deserialize, Debug)]
  pub struct AuspiceTree {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub version: Option<String>,

    pub meta: AuspiceTreeMeta,

    pub tree: AuspiceTreeNode,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub root_sequence: Option<BTreeMap<String, String>>,

    #[serde(flatten)]
    pub other: serde_json::Value,
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
        .root_sequence
        .as_ref()
        .and_then(|root_sequence| root_sequence.get("nuc"))
        .map(String::as_str)
    }
  }
}
