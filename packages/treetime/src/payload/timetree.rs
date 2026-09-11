#[cfg(test)]
mod __tests__;

use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdge;
use treetime_graph::node::GraphNode;
use treetime_io::graphviz::{EdgeToGraphviz, NodeToGraphviz};
use treetime_io::nwk::{EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk};

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct NodeTimetree {}

impl GraphNode for NodeTimetree {}

impl NodeFromNwk for NodeTimetree {
  fn from_nwk(
    _name: Option<impl AsRef<str>>,
    _confidence: Option<f64>,
    _: &BTreeMap<String, String>,
  ) -> Result<Self, Report> {
    Ok(Self {})
  }
}

impl NodeToNwk for NodeTimetree {
  fn nwk_comments(&self) -> BTreeMap<String, String> {
    BTreeMap::new()
  }
}

impl NodeToGraphviz for NodeTimetree {}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct EdgeTimetree {}

impl GraphEdge for EdgeTimetree {}

impl EdgeFromNwk for EdgeTimetree {
  fn from_nwk(_branch_length: Option<f64>) -> Result<Self, Report> {
    Ok(Self {})
  }
}

impl EdgeToNwk for EdgeTimetree {
  fn nwk_weight(&self) -> Option<f64> {
    None
  }
}

impl EdgeToGraphviz for EdgeTimetree {}
