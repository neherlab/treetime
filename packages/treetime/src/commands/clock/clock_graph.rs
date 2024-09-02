use crate::commands::clock::clock_set::ClockSet;
use crate::graph::edge::{GraphEdge, Weighted};
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, Named};
use crate::io::graphviz::{EdgeToGraphViz, NodeToGraphviz};
use crate::io::nwk::{format_weight, EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk, NwkWriteOptions};
use crate::o;
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub type ClockGraph = Graph<ClockNodePayload, ClockEdgePayload, ()>;

#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct ClockNodePayload {
  pub name: Option<String>,
  pub date: Option<f64>,
  pub bad_branch: bool,
  pub div: f64,
  pub is_outlier: bool,
  pub total: ClockSet,
  pub to_parent: ClockSet,
  pub to_children: BTreeMap<String, ClockSet>,
  pub from_children: BTreeMap<String, ClockSet>,
}

impl GraphNode for ClockNodePayload {}

impl NodeFromNwk for ClockNodePayload {
  fn from_nwk(name: Option<impl AsRef<str>>, _: &BTreeMap<String, String>) -> Result<Self, Report> {
    Ok(Self {
      name: name.map(|s| s.as_ref().to_owned()),
      ..ClockNodePayload::default()
    })
  }
}

impl NodeToNwk for ClockNodePayload {
  fn nwk_name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn nwk_comments(&self) -> BTreeMap<String, String> {
    let mutations: String = "".to_owned(); // TODO: fill mutations
    BTreeMap::from([(o!("mutations"), mutations)])
  }
}

impl Named for ClockNodePayload {
  fn name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn set_name(&mut self, name: Option<impl AsRef<str>>) {
    self.name = name.map(|n| n.as_ref().to_owned());
  }
}

impl NodeToGraphviz for ClockNodePayload {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }
}

#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct ClockEdgePayload {
  pub branch_length: Option<f64>,
}

impl GraphEdge for ClockEdgePayload {}

impl Weighted for ClockEdgePayload {
  fn weight(&self) -> Option<f64> {
    self.branch_length
  }

  fn set_weight(&mut self, weight: Option<f64>) {
    self.branch_length = weight;
  }
}

impl EdgeFromNwk for ClockEdgePayload {
  fn from_nwk(branch_length: Option<f64>) -> Result<Self, Report> {
    Ok(Self { branch_length })
  }
}

impl EdgeToNwk for ClockEdgePayload {
  fn nwk_weight(&self) -> Option<f64> {
    self.weight()
  }
}

impl EdgeToGraphViz for ClockEdgePayload {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self
      .weight()
      .map(|weight| format_weight(weight, &NwkWriteOptions::default()))
  }

  fn to_graphviz_weight(&self) -> Option<f64> {
    self.weight()
  }
}
