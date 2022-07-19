use crate::graph::graph::{Graph as GenericGraph, Weighted};
use crate::io::nwk::read_nwk;
use color_eyre::Section;
use eyre::{eyre, Report, WrapErr};
use indexmap::IndexMap;
use ndarray::{array, Array1, Array2};
use std::fmt::{Display, Formatter};
use std::path::Path;

pub type AncestralGraph = GenericGraph<Node, Edge>;

#[derive(Clone, Debug, PartialEq)]
pub enum NodeType {
  Leaf(String),
  Internal(f64),
}

#[derive(Clone, Debug, PartialEq)]
pub struct Node {
  pub name: String,
  pub node_type: NodeType,
  pub joint_Lx: Array2<f64>,   // likelihood
  pub joint_Cx: Array2<usize>, // max likelihood indices
  pub mask: Option<usize>,
  pub seq: Array1<char>,
  pub seq_ii: Array1<usize>,
}

impl Node {
  pub fn leaf(name: &str) -> Self {
    Self {
      name: name.to_owned(),
      node_type: NodeType::Leaf(name.to_owned()),
      joint_Lx: Array2::<f64>::default((0, 0)),
      joint_Cx: Array2::<usize>::default((0, 0)),
      mask: None,
      seq: Array1::<char>::default(0),
      seq_ii: Array1::<usize>::default(0),
    }
  }

  pub fn internal(name: &str, weight: f64) -> Self {
    Self {
      name: name.to_owned(),
      node_type: NodeType::Internal(weight),
      joint_Lx: array![[]],
      joint_Cx: array![[]],
      mask: None,
      seq: Array1::<char>::default(0),
      seq_ii: Array1::<usize>::default(0),
    }
  }
}

impl Display for Node {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    match &self.node_type {
      NodeType::Leaf(name) => write!(f, "{name}"),
      NodeType::Internal(weight) => write!(f, "{weight:1.4}"),
    }
  }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Edge {
  weight: f64,
}

impl Weighted for Edge {}

impl Display for Edge {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    write!(f, "{:1.4}", &self.weight)
  }
}

pub fn infer_graph() -> Result<AncestralGraph, Report> {
  unimplemented!("Not implemented: Graph inference is not yet implemented. Please provide the `--tree` argument.")
}

pub fn create_graph(tree_path: impl AsRef<Path>) -> Result<AncestralGraph, Report> {
  let nwk_tree = read_nwk(tree_path)
    .wrap_err("When parsing input tree")
    .with_section(|| "Note: only Newick format is currently supported")?;

  let mut graph = AncestralGraph::new();

  // Insert nodes
  let mut index_map = IndexMap::<usize, usize>::new(); // Map of internal `nwk` node indices to `Graph` node indices
  let mut node_counter: usize = 0;
  for (nwk_idx, nwk_node) in nwk_tree.g.raw_nodes().iter().enumerate() {
    // Attempt to parse weight as float. If not a float, then it's a named leaf node, otherwise - internal node.
    let inserted_node_idx = if let Ok(weight) = nwk_node.weight.parse::<f64>() {
      let node = Node::internal(&format!("NODE_{node_counter:07}"), weight);
      node_counter += 1;
      graph.add_node(node)
    } else {
      let name = nwk_node.weight.as_str();
      if name == "N/A" {
        let node = Node::internal(&format!("NODE_{node_counter:07}"), 0.0);
        node_counter += 1;
        graph.add_node(node)
      } else {
        graph.add_node(Node::leaf(name))
      }
    };

    index_map.insert(nwk_idx, inserted_node_idx);
  }

  // Insert edges
  for (nwk_idx, nwk_edge) in nwk_tree.g.raw_edges().iter().enumerate() {
    let weight: f64 = nwk_edge.weight as f64;
    let source: usize = nwk_edge.source().index();
    let target: usize = nwk_edge.target().index();

    let source = index_map
      .get(&source)
      .ok_or_else(|| eyre!("When inserting edge {nwk_idx}: Node with index {source} not found."))?;

    let target = index_map
      .get(&target)
      .ok_or_else(|| eyre!("When inserting edge {nwk_idx}: Node with index {target} not found."))?;

    graph.add_edge(*source, *target, Edge { weight });
  }

  graph.build()?;

  Ok(graph)
}
