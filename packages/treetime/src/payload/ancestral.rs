#[cfg(test)]
mod __tests__;

use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use treetime_graph::edge::{GraphEdge, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_io::graphviz::{EdgeToGraphviz, NodeToGraphviz};
use treetime_io::nwk::{EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk};

pub type GraphAncestral<D = ()> = Graph<NodeAncestral, EdgeAncestral, D>;

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct NodeAncestral {}

impl NodeFromNwk for NodeAncestral {
  fn from_nwk(
    _name: Option<impl AsRef<str>>,
    _confidence: Option<f64>,
    _: &BTreeMap<String, String>,
  ) -> Result<Self, Report> {
    Ok(Self {})
  }
}

impl NodeToNwk for NodeAncestral {
  fn nwk_comments(&self) -> BTreeMap<String, String> {
    BTreeMap::new()
  }
}

impl GraphNode for NodeAncestral {}

impl NodeToGraphviz for NodeAncestral {}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct EdgeAncestral {
  pub branch_length: Option<f64>,
}

impl GraphEdge for EdgeAncestral {}

impl HasBranchLength for EdgeAncestral {
  fn branch_length(&self) -> Option<f64> {
    self.branch_length
  }

  fn set_branch_length(&mut self, weight: Option<f64>) {
    self.branch_length = weight;
  }
}

impl EdgeFromNwk for EdgeAncestral {
  fn from_nwk(branch_length: Option<f64>) -> Result<Self, Report> {
    Ok(Self { branch_length })
  }
}

impl EdgeToNwk for EdgeAncestral {
  fn nwk_weight(&self) -> Option<f64> {
    self.branch_length()
  }
}

impl EdgeToGraphviz for EdgeAncestral {}
