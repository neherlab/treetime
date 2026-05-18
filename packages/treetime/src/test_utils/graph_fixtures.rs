use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use treetime_graph::edge::{GraphEdge, HasBranchLength};
use treetime_graph::node::{GraphNode, Named};
use treetime_io::graphviz::{EdgeToGraphviz, NodeToGraphviz};
use treetime_io::nwk::{EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk, NwkWriteOptions, format_weight};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct TestNode(pub Option<String>);

impl GraphNode for TestNode {}

impl Named for TestNode {
  fn name(&self) -> Option<impl AsRef<str>> {
    self.0.as_deref()
  }
  fn set_name(&mut self, name: Option<impl AsRef<str>>) {
    self.0 = name.map(|n| n.as_ref().to_owned());
  }
}

impl NodeFromNwk for TestNode {
  fn from_nwk(name: Option<impl AsRef<str>>, _: &BTreeMap<String, String>) -> Result<Self, Report> {
    Ok(Self(name.map(|n| n.as_ref().to_owned())))
  }
}

impl NodeToNwk for TestNode {
  fn nwk_name(&self) -> Option<impl AsRef<str>> {
    self.0.as_deref()
  }
}

impl NodeToGraphviz for TestNode {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self.0.as_deref()
  }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TestEdge(pub Option<f64>);

impl GraphEdge for TestEdge {}

impl HasBranchLength for TestEdge {
  fn branch_length(&self) -> Option<f64> {
    self.0
  }
  fn set_branch_length(&mut self, weight: Option<f64>) {
    self.0 = weight;
  }
}

impl EdgeFromNwk for TestEdge {
  fn from_nwk(weight: Option<f64>) -> Result<Self, Report> {
    Ok(Self(weight))
  }
}

impl EdgeToNwk for TestEdge {
  fn nwk_weight(&self) -> Option<f64> {
    self.0
  }
}

impl EdgeToGraphviz for TestEdge {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self.0.map(|weight| format_weight(weight, &NwkWriteOptions::default()))
  }

  fn to_graphviz_weight(&self) -> Option<f64> {
    self.0
  }
}
