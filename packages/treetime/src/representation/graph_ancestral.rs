use crate::commands::clock::clock_set::ClockSet;
use crate::commands::clock::clock_traits::{ClockEdge, ClockNode};
use crate::commands::clock::date_constraints::DateConstraintNode;
use crate::distribution::distribution::Distribution;
use crate::graph::edge::{GraphEdge, Weighted};
use crate::graph::graph::Graph;
use crate::graph::node::{Described, GraphNode, Named};
use crate::io::graphviz::{EdgeToGraphViz, NodeToGraphviz};
use crate::io::nwk::{EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk, NwkWriteOptions, format_weight};
use crate::o;
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::sync::Arc;

pub type GraphAncestral = Graph<NodeAncestral, EdgeAncestral, ()>;

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct NodeAncestral {
  pub name: Option<String>,
  pub desc: Option<String>,
  pub time: Option<f64>,
  pub time_before_present: Option<f64>,
  pub time_distribution: Option<Arc<Distribution>>,
  pub bad_branch: bool,
  pub div: f64,
  pub is_outlier: bool,
  pub clock_set: ClockSet,
}

impl NodeFromNwk for NodeAncestral {
  fn from_nwk(name: Option<impl AsRef<str>>, _: &BTreeMap<String, String>) -> Result<Self, Report> {
    Ok(Self {
      name: name.map(|s| s.as_ref().to_owned()),
      ..NodeAncestral::default()
    })
  }
}

impl NodeToNwk for NodeAncestral {
  fn nwk_name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn nwk_comments(&self) -> BTreeMap<String, String> {
    let mut comments = BTreeMap::new();
    let mutations: String = "".to_owned(); // TODO: fill mutations
    comments.insert(o!("mutations"), mutations);
    if let Some(time) = self.time {
      comments.insert("date".to_owned(), format!("{time:.2}"));
    }
    comments
  }
}

impl GraphNode for NodeAncestral {}

impl Named for NodeAncestral {
  fn name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn set_name(&mut self, name: Option<impl AsRef<str>>) {
    self.name = name.map(|n| n.as_ref().to_owned());
  }
}

impl Described for NodeAncestral {
  fn desc(&self) -> &Option<String> {
    &self.desc
  }

  fn set_desc(&mut self, desc: Option<String>) {
    self.desc = desc;
  }
}

impl DateConstraintNode for NodeAncestral {
  fn get_time_distribution(&self) -> &Option<Arc<Distribution>> {
    &self.time_distribution
  }

  fn set_time_distribution(&mut self, dist: Option<Arc<Distribution>>) {
    self.time_distribution = dist;
  }

  fn get_bad_branch(&self) -> bool {
    self.bad_branch
  }

  fn set_bad_branch(&mut self, bad: bool) {
    self.bad_branch = bad;
  }
}

impl NodeToGraphviz for NodeAncestral {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct EdgeAncestral {
  pub branch_length: Option<f64>,
  pub branch_length_distribution: Option<Arc<Distribution>>,
  pub msg_to_parent: Option<Arc<Distribution>>,
  #[serde(skip)]
  pub clock_to_parent: ClockSet,
  #[serde(skip)]
  pub clock_to_child: ClockSet,
  #[serde(skip)]
  pub clock_from_child: ClockSet,
}

impl GraphEdge for EdgeAncestral {}

impl Weighted for EdgeAncestral {
  fn weight(&self) -> Option<f64> {
    self.branch_length
  }

  fn set_weight(&mut self, weight: Option<f64>) {
    self.branch_length = weight;
  }
}

impl EdgeFromNwk for EdgeAncestral {
  fn from_nwk(branch_length: Option<f64>) -> Result<Self, Report> {
    Ok(Self {
      branch_length,
      ..Self::default()
    })
  }
}

impl EdgeToNwk for EdgeAncestral {
  fn nwk_weight(&self) -> Option<f64> {
    self.weight()
  }
}

impl EdgeToGraphViz for EdgeAncestral {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self
      .weight()
      .map(|weight| format_weight(weight, &NwkWriteOptions::default()))
  }

  fn to_graphviz_weight(&self) -> Option<f64> {
    self.weight()
  }
}

impl ClockNode for NodeAncestral {
  fn likely_time(&self) -> Option<f64> {
    self.time_distribution.as_ref().and_then(|dist| dist.likely_time())
  }

  fn div(&self) -> f64 {
    self.div
  }

  fn set_div(&mut self, div: f64) {
    self.div = div;
  }

  fn is_outlier(&self) -> bool {
    self.is_outlier
  }

  fn clock_set(&self) -> &ClockSet {
    &self.clock_set
  }

  fn clock_set_mut(&mut self) -> &mut ClockSet {
    &mut self.clock_set
  }
}

impl ClockEdge for EdgeAncestral {
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
