#[cfg(test)]
mod __tests__;

use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdge;
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
pub struct EdgeAncestral {}

impl GraphEdge for EdgeAncestral {}

impl EdgeFromNwk for EdgeAncestral {
  fn from_nwk(_branch_length: Option<f64>) -> Result<Self, Report> {
    Ok(Self {})
  }
}

impl EdgeToNwk for EdgeAncestral {
  fn nwk_weight(&self) -> Option<f64> {
    None
  }
}

impl EdgeToGraphviz for EdgeAncestral {}
