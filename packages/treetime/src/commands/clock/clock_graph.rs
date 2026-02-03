use crate::commands::clock::clock_set::ClockSet;
use crate::commands::clock::clock_traits::{ClockEdge, ClockNode};
use crate::graph::edge::{GraphEdge, Weighted};
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, Named, Outlier};
use crate::io::graphviz::{EdgeToGraphViz, NodeToGraphviz};
use crate::io::nwk::{EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk, NwkWriteOptions, format_weight};
use crate::o;
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub type GraphClock = Graph<NodeClock, EdgeClock, ()>;

#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct NodeClock {
  pub name: Option<String>,
  pub date: Option<f64>,
  pub bad_branch: bool,
  pub div: f64,
  pub is_outlier: bool,
  pub clock_set: ClockSet,
}

impl GraphNode for NodeClock {}

impl NodeFromNwk for NodeClock {
  fn from_nwk(name: Option<impl AsRef<str>>, _: &BTreeMap<String, String>) -> Result<Self, Report> {
    Ok(Self {
      name: name.map(|s| s.as_ref().to_owned()),
      ..NodeClock::default()
    })
  }
}

impl NodeToNwk for NodeClock {
  fn nwk_name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn nwk_comments(&self) -> BTreeMap<String, String> {
    let mutations: String = "".to_owned(); // TODO: fill mutations
    BTreeMap::from([(o!("mutations"), mutations)])
  }
}

impl Named for NodeClock {
  fn name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn set_name(&mut self, name: Option<impl AsRef<str>>) {
    self.name = name.map(|n| n.as_ref().to_owned());
  }
}

impl Outlier for NodeClock {
  fn is_outlier(&self) -> bool {
    self.is_outlier
  }

  fn set_is_outlier(&mut self, is_outlier: bool) {
    self.is_outlier = is_outlier;
  }
}

impl NodeToGraphviz for NodeClock {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }
}

#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct EdgeClock {
  pub branch_length: Option<f64>,
  pub clock_to_parent: ClockSet,
  pub clock_to_child: ClockSet,
  pub clock_from_child: ClockSet, // this is the propagated 'to_parent' msg. only need to avoid recalculation of propagated message
}

impl GraphEdge for EdgeClock {}

impl Weighted for EdgeClock {
  fn weight(&self) -> Option<f64> {
    self.branch_length
  }

  fn set_weight(&mut self, weight: Option<f64>) {
    self.branch_length = weight;
  }
}

impl EdgeFromNwk for EdgeClock {
  fn from_nwk(branch_length: Option<f64>) -> Result<Self, Report> {
    Ok(Self {
      branch_length,
      ..EdgeClock::default()
    })
  }
}

impl EdgeToNwk for EdgeClock {
  fn nwk_weight(&self) -> Option<f64> {
    self.weight()
  }
}

impl EdgeToGraphViz for EdgeClock {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self
      .weight()
      .map(|weight| format_weight(weight, &NwkWriteOptions::default()))
  }

  fn to_graphviz_weight(&self) -> Option<f64> {
    self.weight()
  }
}

impl ClockNode for NodeClock {
  fn likely_time(&self) -> Option<f64> {
    self.date
  }

  fn div(&self) -> f64 {
    self.div
  }

  fn set_div(&mut self, div: f64) {
    self.div = div;
  }

  fn clock_set(&self) -> &ClockSet {
    &self.clock_set
  }

  fn clock_set_mut(&mut self) -> &mut ClockSet {
    &mut self.clock_set
  }
}

impl ClockEdge for EdgeClock {
  fn branch_length(&self) -> Option<f64> {
    self.branch_length
  }

  fn set_branch_length(&mut self, length: Option<f64>) {
    self.branch_length = length;
  }

  fn to_parent(&self) -> &ClockSet {
    &self.clock_to_parent
  }

  fn to_parent_mut(&mut self) -> &mut ClockSet {
    &mut self.clock_to_parent
  }

  fn to_child(&self) -> &ClockSet {
    &self.clock_to_child
  }

  fn to_child_mut(&mut self) -> &mut ClockSet {
    &mut self.clock_to_child
  }

  fn from_child(&self) -> &ClockSet {
    &self.clock_from_child
  }

  fn from_child_mut(&mut self) -> &mut ClockSet {
    &mut self.clock_from_child
  }
}
