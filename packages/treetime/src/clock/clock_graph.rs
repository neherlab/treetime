use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_io::graphviz::{EdgeToGraphviz, NodeToGraphviz};
use treetime_io::nwk::{EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk};

pub type GraphClock<D = ()> = Graph<NodeClock, EdgeClock, D>;

#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct NodeClock {}

impl GraphNode for NodeClock {}

impl NodeFromNwk for NodeClock {
  fn from_nwk(
    _name: Option<impl AsRef<str>>,
    _confidence: Option<f64>,
    _: &BTreeMap<String, String>,
  ) -> Result<Self, Report> {
    Ok(Self {})
  }
}

impl NodeToNwk for NodeClock {
  fn nwk_comments(&self) -> BTreeMap<String, String> {
    BTreeMap::new()
  }
}

impl NodeToGraphviz for NodeClock {}

#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct EdgeClock {}

impl GraphEdge for EdgeClock {}

impl EdgeFromNwk for EdgeClock {
  fn from_nwk(_branch_length: Option<f64>) -> Result<Self, Report> {
    Ok(Self {})
  }
}

impl EdgeToNwk for EdgeClock {
  fn nwk_weight(&self) -> Option<f64> {
    None
  }
}

impl EdgeToGraphviz for EdgeClock {}
