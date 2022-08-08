use crate::graph::graph::{Graph as GenericGraph, Weighted};
use crate::io::nwk::read_nwk;
use color_eyre::Section;
use eyre::{eyre, Report, WrapErr};
use indexmap::IndexMap;
use ndarray::{array, Array1};
use std::fmt::{Display, Formatter};
use std::path::Path;

pub type ClockGraph = GenericGraph<Node, Edge>;

// {'O': array([    36153.374401089997263625,         0.401970000000000049,  72615145.592690706253051758,       806.093643935642603537,         0.015141327100000003,        18.000000000000000000]),
//  '_color': None,
//  '_ii': array([18]),
//  '_v': 0.04481212509893939,
//  'bad_branch': False,
//  'branch_length': 0.00142,
//  'clades': [],
//  'clock_length': 0.00142,
//  'comment': None,
//  'confidence': None,
//  'dist2root': 0.04481212509893939,
//  'mask': None,
//  'mutation_length': 0.00142,
//  'name': 'A/Boston/DOA2_107/2012|CY148382|11/01/2012|USA|12_13|H3N2/1-1409',
//  'original_length': 0.00142,
//  'raw_date_constraint': 2012.83778234,
//  'tt': <treetime.treetime.TreeTime object at 0x7fd6fcd78760>,
//  'up': Clade(bad_branch=False, branch_length=0.00055, clock_length=0.00055, confidence=0.45, mutation_length=0.00055, name='NODE_0000002', original_length=0.00055),
//  'width': None}

#[derive(Clone, Debug, PartialEq)]
pub enum NodeType {
  Root,
  Internal,
  Leaf(String),
}

#[derive(Clone, Debug, PartialEq)]
pub struct Node {
  pub name: String,
  pub node_type: NodeType,
  pub bad_branch: bool,
  pub original_length: f64,
  pub branch_length: f64,
  pub mutation_length: f64,
  pub dist2root: f64,
  pub raw_date_constraint: f64,
  pub Q: Array1<f64>,
  pub Qtot: Array1<f64>,
  pub O: Array1<f64>,
}

impl Node {
  pub fn leaf(name: &str) -> Self {
    Self {
      name: name.to_owned(),
      node_type: NodeType::Leaf(name.to_owned()),
      bad_branch: false,
      original_length: 0.0,
      branch_length: 0.0,
      mutation_length: 0.0,
      dist2root: 0.0,
      raw_date_constraint: 0.0,
      Q: Array1::zeros(6),
      Qtot: Array1::zeros(6),
      O: Array1::zeros(6),
    }
  }

  pub fn internal(name: &str) -> Self {
    Self {
      name: name.to_owned(),
      node_type: NodeType::Internal,
      bad_branch: false,
      original_length: 0.0,
      branch_length: 0.0,
      mutation_length: 0.0,
      dist2root: 0.0,
      raw_date_constraint: 0.0,
      Q: Array1::zeros(6),
      Qtot: Array1::zeros(6),
      O: Array1::zeros(6),
    }
  }
}

impl Display for Node {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    write!(f, "{}", self.name)
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

pub fn infer_graph() -> Result<ClockGraph, Report> {
  unimplemented!("Not implemented: Graph inference is not yet implemented. Please provide the `--tree` argument.")
}

pub fn create_graph(tree_path: impl AsRef<Path>) -> Result<ClockGraph, Report> {
  let nwk_tree = read_nwk(tree_path)
    .wrap_err("When parsing input tree")
    .with_section(|| "Note: only Newick format is currently supported")?;

  let mut graph = ClockGraph::new();

  // Insert nodes
  let mut index_map = IndexMap::<usize, usize>::new(); // Map of internal `nwk` node indices to `Graph` node indices
  let mut node_counter: usize = 0;
  for (nwk_idx, nwk_node) in nwk_tree.g.raw_nodes().iter().enumerate() {
    // Attempt to parse weight as float. If not a float, then it's a named leaf node, otherwise - internal node.
    let inserted_node_idx = if let Ok(weight) = nwk_node.weight.parse::<f64>() {
      let node = Node::internal(&format!("NODE_{node_counter:07}"));
      node_counter += 1;
      graph.add_node(node)
    } else {
      let name = nwk_node.weight.as_str();
      if name == "N/A" {
        let node = Node::internal(&format!("NODE_{node_counter:07}"));
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

  for root in graph.get_roots() {
    root.write().payload().write().node_type = NodeType::Root;
  }

  Ok(graph)
}
