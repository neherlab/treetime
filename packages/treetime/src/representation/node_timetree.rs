use crate::commands::clock::clock_set::ClockSet;
use crate::commands::clock::clock_traits::ClockNode;
use crate::commands::timetree::data::date_constraints::DateConstraintNode;
use crate::commands::timetree::timetree_traits::TimetreeNode;
use crate::distribution::distribution::Distribution;
use crate::graph::node::{Described, GraphNode, Named};
use crate::io::graphviz::NodeToGraphviz;
use crate::io::nwk::{NodeFromNwk, NodeToNwk};
use crate::representation::graph_ancestral::NodeAncestral;
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::sync::Arc;

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct NodeTimetree {
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

impl GraphNode for NodeTimetree {}

impl Named for NodeTimetree {
  fn name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn set_name(&mut self, name: Option<impl AsRef<str>>) {
    self.name = name.map(|n| n.as_ref().to_owned());
  }
}

impl Described for NodeTimetree {
  fn desc(&self) -> &Option<String> {
    &self.desc
  }

  fn set_desc(&mut self, desc: Option<String>) {
    self.desc = desc;
  }
}

impl ClockNode for NodeTimetree {
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

impl DateConstraintNode for NodeTimetree {
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

impl NodeFromNwk for NodeTimetree {
  fn from_nwk(name: Option<impl AsRef<str>>, _: &BTreeMap<String, String>) -> Result<Self, Report> {
    Ok(Self {
      name: name.map(|s| s.as_ref().to_owned()),
      ..NodeTimetree::default()
    })
  }
}

impl NodeToNwk for NodeTimetree {
  fn nwk_name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn nwk_comments(&self) -> BTreeMap<String, String> {
    let mut comments = BTreeMap::new();
    if let Some(time) = self.time {
      comments.insert("date".to_owned(), format!("{time:.2}"));
    }
    comments
  }
}

impl NodeToGraphviz for NodeTimetree {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }
}

impl TimetreeNode for NodeTimetree {
  fn time_distribution(&self) -> &Option<Arc<Distribution>> {
    &self.time_distribution
  }

  fn set_time_distribution(&mut self, dist: Option<Arc<Distribution>>) {
    self.time_distribution = dist;
  }

  fn time(&self) -> Option<f64> {
    self.time
  }

  fn set_time(&mut self, time: Option<f64>) {
    self.time = time;
  }
}

impl From<&NodeAncestral> for NodeTimetree {
  fn from(node: &NodeAncestral) -> Self {
    Self {
      name: node.name.clone(),
      desc: node.desc.clone(),
      time: node.time,
      time_before_present: node.time_before_present,
      time_distribution: node.time_distribution.clone(),
      bad_branch: node.bad_branch,
      div: node.div,
      is_outlier: node.is_outlier,
      clock_set: node.clock_set.clone(),
    }
  }
}
