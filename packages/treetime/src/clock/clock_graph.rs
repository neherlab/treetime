use crate::graph::graph::{Graph as GenericGraph, GraphNodeBackward, Weighted};
use crate::io::dates::{DateOrRange, DatesMap};
use crate::io::nwk::read_nwk;
use crate::make_error;
use color_eyre::Section;
use eyre::{eyre, Report, WrapErr};
use indexmap::IndexMap;
use ndarray::Array1;
use std::fmt::{Display, Formatter};
use std::path::Path;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering::Relaxed;

pub type ClockGraph = GenericGraph<Node, Edge>;

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
  pub raw_date_constraint: Option<f64>,
  pub Q: Array1<f64>,
  pub Qtot: Array1<f64>,
  pub O: Array1<f64>,
  pub v: f64,
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
      raw_date_constraint: None,
      Q: Array1::zeros(6),
      Qtot: Array1::zeros(6),
      O: Array1::zeros(6),
      v: 0.0,
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
      raw_date_constraint: None,
      Q: Array1::zeros(6),
      Qtot: Array1::zeros(6),
      O: Array1::zeros(6),
      v: 0.0,
    }
  }

  #[inline]
  pub const fn is_root(&self) -> bool {
    matches!(self.node_type, NodeType::Root)
  }

  #[inline]
  pub const fn is_leaf(&self) -> bool {
    matches!(self.node_type, NodeType::Leaf(_))
  }
}

impl Display for Node {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    write!(f, "{}", self.name)
  }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Edge {
  pub weight: f64,
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

pub fn create_graph(tree_path: impl AsRef<Path>, dates: &DatesMap) -> Result<ClockGraph, Report> {
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

  // Mark roots
  for root in graph.get_roots() {
    root.write().payload().write().node_type = NodeType::Root;
  }

  assign_dates(&mut graph, dates)?;

  Ok(graph)
}

pub fn assign_dates(graph: &mut ClockGraph, dates: &DatesMap) -> Result<(), Report> {
  // Add dates and mark bad branches
  let num_bad_branches = AtomicUsize::new(0);
  let num_dates = AtomicUsize::new(0);
  graph.par_iter_breadth_first_backward(
    |GraphNodeBackward {
       is_root,
       is_leaf,
       key,
       payload: node,
       children,
     }| {
      node.raw_date_constraint = dates.get(&node.name).map(DateOrRange::mean);
      if node.raw_date_constraint.is_none() {
        if node.is_leaf() {
          // Terminal nodes without date constraints marked as 'bad'
          node.bad_branch = true;
          num_bad_branches.fetch_add(1, Relaxed);
        } else {
          // If all branches downstream are 'bad', and there is no date constraint for this node, the branch is marked as 'bad'
          node.bad_branch = children.iter().all(|child| child.node.read().bad_branch);
        }
      } else {
        num_dates.fetch_add(1, Relaxed);
      }
    },
  );

  let num_bad_branches = num_bad_branches.into_inner();
  let num_dates = num_dates.into_inner();
  if num_bad_branches > graph.num_leaves() - 3 {
    return make_error!(
      "Not enough date constraints: total leaves {}, total dates: {}, bad branches: {}",
      graph.num_leaves(),
      num_dates,
      num_bad_branches
    );
  };

  Ok(())
}
