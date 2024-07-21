use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::edge::{GraphEdge, Weighted};
use crate::graph::graph::{Graph, GraphNodeBackward, GraphNodeForward};
use crate::graph::node::{GraphNode, Named};
use crate::io::dates_csv::{DateOrRange, DatesMap};
use crate::io::graphviz::{EdgeToGraphViz, NodeToGraphviz};
use crate::io::nwk::{format_weight, nwk_read_file, EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk, NwkWriteOptions};
use crate::make_error;
use eyre::Report;
use maplit::btreemap;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::Path;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering::Relaxed;

pub type ClockGraph = Graph<Node, Edge>;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Node {
  pub name: String,
  pub bad_branch: bool,
  pub dist2root: f64,
  pub raw_date_constraint: Option<f64>,
  pub Q: Array1<f64>,
  pub Qtot: Array1<f64>,
  pub O: Array1<f64>,
  pub v: f64,
}

impl GraphNode for Node {}

impl Named for Node {
  fn name(&self) -> &str {
    &self.name
  }

  fn set_name(&mut self, name: impl AsRef<str>) {
    self.name = name.as_ref().to_owned();
  }
}

impl NodeFromNwk for Node {
  fn from_nwk(name: impl AsRef<str>, _: &BTreeMap<String, String>) -> Result<Self, Report> {
    Ok(Self {
      name: name.as_ref().to_owned(),
      bad_branch: false,
      dist2root: 0.0,
      raw_date_constraint: None,
      Q: Array1::zeros(6),
      Qtot: Array1::zeros(6),
      O: Array1::zeros(6),
      v: 0.0,
    })
  }
}

impl NodeToNwk for Node {
  fn nwk_name(&self) -> Option<impl AsRef<str>> {
    Some(&self.name)
  }

  fn nwk_comments(&self) -> BTreeMap<String, String> {
    btreemap! {}
  }
}

impl NodeToGraphviz for Node {
  fn to_graphviz_label(&self) -> String {
    self.name.clone()
  }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Edge {
  pub weight: f64,
}

impl GraphEdge for Edge {}

impl Weighted for Edge {
  fn weight(&self) -> f64 {
    self.weight
  }
}

impl EdgeFromNwk for Edge {
  fn from_nwk(weight: Option<f64>) -> Result<Self, Report> {
    Ok(Self {
      weight: weight.unwrap_or_default(),
    })
  }
}

impl EdgeToNwk for Edge {
  fn nwk_weight(&self) -> Option<f64> {
    Some(self.weight)
  }
}

impl EdgeToGraphViz for Edge {
  fn to_graphviz_label(&self) -> String {
    format_weight(self.weight, &NwkWriteOptions::default())
  }

  fn to_graphviz_weight(&self) -> f64 {
    self.weight
  }
}

pub fn infer_graph() -> Result<ClockGraph, Report> {
  unimplemented!("Not implemented: Graph inference is not yet implemented. Please provide the `--tree` argument.")
}

pub fn create_graph<P: AsRef<Path>>(tree_path: P, dates: &DatesMap) -> Result<ClockGraph, Report> {
  let mut graph = nwk_read_file(tree_path)?;
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
       ..
     }| {
      node.raw_date_constraint = dates.get(&node.name).map(DateOrRange::mean);
      if node.raw_date_constraint.is_none() {
        if is_leaf {
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
       is_root,
       ..
     }| {
      if is_root {
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
