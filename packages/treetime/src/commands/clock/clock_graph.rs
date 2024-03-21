use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::create_graph_from_nwk::create_graph_from_nwk_file;
use crate::graph::edge::{GraphEdge, Weighted};
use crate::graph::graph::{Graph, GraphNodeBackward, GraphNodeForward};
use crate::graph::node::{GraphNode, Named, NodeType, WithNwkComments};
use crate::io::dates_csv::{DateOrRange, DatesMap};
use crate::make_error;
use eyre::Report;
use ndarray::Array1;
use std::fmt::{Display, Formatter};
use std::path::Path;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering::Relaxed;

pub type ClockGraph = Graph<Node, Edge>;

#[derive(Clone, Debug, PartialEq)]
pub struct Node {
  pub name: String,
  pub node_type: NodeType,
  pub bad_branch: bool,
  pub dist2root: f64,
  pub raw_date_constraint: Option<f64>,
  pub Q: Array1<f64>,
  pub Qtot: Array1<f64>,
  pub O: Array1<f64>,
  pub v: f64,
}

impl GraphNode for Node {
  fn root(name: &str) -> Self {
    Self {
      name: name.to_owned(),
      node_type: NodeType::Root(name.to_owned()),
      bad_branch: false,
      dist2root: 0.0,
      raw_date_constraint: None,
      Q: Array1::zeros(6),
      Qtot: Array1::zeros(6),
      O: Array1::zeros(6),
      v: 0.0,
    }
  }

  fn internal(name: &str) -> Self {
    Self {
      name: name.to_owned(),
      node_type: NodeType::Internal(name.to_owned()),
      bad_branch: false,
      dist2root: 0.0,
      raw_date_constraint: None,
      Q: Array1::zeros(6),
      Qtot: Array1::zeros(6),
      O: Array1::zeros(6),
      v: 0.0,
    }
  }

  fn leaf(name: &str) -> Self {
    Self {
      name: name.to_owned(),
      node_type: NodeType::Leaf(name.to_owned()),
      bad_branch: false,
      dist2root: 0.0,
      raw_date_constraint: None,
      Q: Array1::zeros(6),
      Qtot: Array1::zeros(6),
      O: Array1::zeros(6),
      v: 0.0,
    }
  }

  fn set_node_type(&mut self, node_type: NodeType) {
    self.node_type = node_type;
  }
}

impl WithNwkComments for Node {}

impl Named for Node {
  fn name(&self) -> &str {
    &self.name
  }

  fn set_name(&mut self, name: &str) {
    self.name = name.to_owned();
  }
}

impl Node {
  #[inline]
  pub const fn is_root(&self) -> bool {
    matches!(self.node_type, NodeType::Root(_))
  }

  #[inline]
  pub const fn is_internal(&self) -> bool {
    matches!(self.node_type, NodeType::Internal(_))
  }

  #[inline]
  pub const fn is_leaf(&self) -> bool {
    matches!(self.node_type, NodeType::Leaf(_))
  }
}

impl Display for Node {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    match &self.node_type {
      NodeType::Root(weight) => write!(f, "{weight:1.4}"),
      NodeType::Internal(weight) => write!(f, "{weight:1.4}"),
      NodeType::Leaf(name) => write!(f, "{name}"),
    }
  }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Edge {
  pub weight: f64,
}

impl GraphEdge for Edge {
  fn new(weight: f64) -> Self {
    Self { weight }
  }
}

impl Weighted for Edge {
  fn weight(&self) -> f64 {
    self.weight
  }
}

impl Display for Edge {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    write!(f, "{:1.4}", &self.weight)
  }
}

pub fn infer_graph() -> Result<ClockGraph, Report> {
  unimplemented!("Not implemented: Graph inference is not yet implemented. Please provide the `--tree` argument.")
}

pub fn create_graph<P: AsRef<Path>>(tree_path: P, dates: &DatesMap) -> Result<ClockGraph, Report> {
  let mut graph = create_graph_from_nwk_file::<Node, Edge>(tree_path)?;
  assign_dates(&mut graph, dates)?;
  calculate_distances_to_root(&mut graph);
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
       payload: mut node,
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
          node.bad_branch = children.iter().all(|(child, _)| child.read().bad_branch);
        }
      } else {
        num_dates.fetch_add(1, Relaxed);
      }
      GraphTraversalContinuation::Continue
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

/// For every node, calculate distance of that node to the root
pub fn calculate_distances_to_root(graph: &mut Graph<Node, Edge>) {
  graph.par_iter_breadth_first_forward(
    |GraphNodeForward {
       payload: mut node,
       parents,
       ..
     }| {
      if node.is_root() {
        node.dist2root = 0.0;
      } else {
        assert!(parents.len() <= 1, "Multiple parents are not supported yet");

        node.dist2root = parents
          .iter()
          .map(|(parent, edge)| {
            let parent = parent.read();
            let edge = edge.read();
            parent.dist2root + edge.weight
          })
          .next().unwrap_or(0.0) // Returns first item or 0.0
          ;
      }
      GraphTraversalContinuation::Continue
    },
  );
}
