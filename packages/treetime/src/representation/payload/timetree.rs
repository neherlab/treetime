use crate::commands::clock::clock_set::ClockSet;
use crate::commands::clock::clock_traits::{ClockEdge, ClockNode};
use crate::commands::clock::date_constraints::DateConstraintNode;
use crate::commands::timetree::timetree_traits::{TimetreeEdge, TimetreeNode};
use crate::representation::payload::ancestral::{EdgeAncestral, NodeAncestral};
use eyre::Report;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_graph::edge::{BranchDistribution, ClockMessages, GraphEdge, HasBranchLength, TimeLength};
use treetime_graph::node::{Described, Divergence, GraphNode, Named, Outlier, TimeConstraint};
use treetime_io::graphviz::{EdgeToGraphviz, NodeToGraphviz};
use treetime_io::nwk::{EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk, NwkWriteOptions, format_weight};

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct NodeTimetree {
  pub base: NodeAncestral,
  pub time: Option<f64>,
  pub time_distribution: Option<Arc<Distribution>>,
  pub bad_branch: bool,
  pub div: f64,
  pub is_outlier: bool,
  pub clock_set: ClockSet,
}

impl GraphNode for NodeTimetree {}

impl Named for NodeTimetree {
  fn name(&self) -> Option<impl AsRef<str>> {
    self.base.name()
  }

  fn set_name(&mut self, name: Option<impl AsRef<str>>) {
    self.base.set_name(name);
  }
}

impl Described for NodeTimetree {
  fn desc(&self) -> &Option<String> {
    self.base.desc()
  }

  fn set_desc(&mut self, desc: Option<String>) {
    self.base.set_desc(desc);
  }
}

impl Divergence for NodeTimetree {
  fn div(&self) -> Option<f64> {
    Some(self.div)
  }

  fn set_div(&mut self, div: Option<f64>) {
    if let Some(div) = div {
      self.div = div;
    }
  }
}

impl Outlier for NodeTimetree {
  fn is_outlier(&self) -> bool {
    self.is_outlier
  }

  fn set_is_outlier(&mut self, is_outlier: bool) {
    self.is_outlier = is_outlier;
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

  fn clock_set(&self) -> &ClockSet {
    &self.clock_set
  }

  fn clock_set_mut(&mut self) -> &mut ClockSet {
    &mut self.clock_set
  }
}

impl TimeConstraint<Arc<Distribution>> for NodeTimetree {
  fn time_distribution(&self) -> &Option<Arc<Distribution>> {
    &self.time_distribution
  }

  fn set_time_distribution(&mut self, dist: Option<Arc<Distribution>>) {
    self.time_distribution = dist;
  }

  fn bad_branch(&self) -> bool {
    self.bad_branch
  }

  fn set_bad_branch(&mut self, bad: bool) {
    self.bad_branch = bad;
  }
}

impl DateConstraintNode for NodeTimetree {}

impl NodeFromNwk for NodeTimetree {
  fn from_nwk(name: Option<impl AsRef<str>>, comments: &BTreeMap<String, String>) -> Result<Self, Report> {
    Ok(Self {
      base: NodeAncestral::from_nwk(name, comments)?,
      ..NodeTimetree::default()
    })
  }
}

impl NodeToNwk for NodeTimetree {
  fn nwk_name(&self) -> Option<impl AsRef<str>> {
    self.base.nwk_name()
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
    self.base.name()
  }
}

impl TimetreeNode for NodeTimetree {
  fn time(&self) -> Option<f64> {
    self.time
  }

  fn set_time(&mut self, time: Option<f64>) {
    self.time = time;
  }
}

#[derive(Clone, SmartDefault, Debug, Serialize, Deserialize)]
pub struct EdgeTimetree {
  pub base: EdgeAncestral,
  pub time_length: Option<f64>,
  pub branch_length_distribution: Option<Arc<Distribution>>,
  pub msg_to_parent: Option<Arc<Distribution>>,
  /// Branch-specific rate multiplier for relaxed molecular clock.
  /// Default 1.0 means branch evolves at the average clock rate.
  /// Values > 1.0 indicate faster evolution, < 1.0 slower.
  #[default = 1.0]
  pub gamma: f64,
  #[serde(skip)]
  pub clock_to_parent: ClockSet,
  #[serde(skip)]
  pub clock_to_child: ClockSet,
  #[serde(skip)]
  pub clock_from_child: ClockSet,
}

impl GraphEdge for EdgeTimetree {}

impl HasBranchLength for EdgeTimetree {
  fn branch_length(&self) -> Option<f64> {
    self.base.branch_length()
  }

  fn set_branch_length(&mut self, weight: Option<f64>) {
    self.base.set_branch_length(weight);
  }
}

impl ClockEdge for EdgeTimetree {
  fn gamma(&self) -> f64 {
    self.gamma
  }
}

impl ClockMessages<ClockSet> for EdgeTimetree {
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

impl BranchDistribution<Arc<Distribution>> for EdgeTimetree {
  fn branch_length_distribution(&self) -> &Option<Arc<Distribution>> {
    &self.branch_length_distribution
  }

  fn set_branch_length_distribution(&mut self, dist: Option<Arc<Distribution>>) {
    self.branch_length_distribution = dist;
  }

  fn msg_to_parent(&self) -> &Option<Arc<Distribution>> {
    &self.msg_to_parent
  }

  fn set_msg_to_parent(&mut self, msg: Option<Arc<Distribution>>) {
    self.msg_to_parent = msg;
  }
}

impl TimeLength for EdgeTimetree {
  fn time_length(&self) -> Option<f64> {
    self.time_length
  }

  fn set_time_length(&mut self, length: Option<f64>) {
    self.time_length = length;
  }
}

impl EdgeFromNwk for EdgeTimetree {
  fn from_nwk(branch_length: Option<f64>) -> Result<Self, Report> {
    Ok(Self {
      base: EdgeAncestral::from_nwk(branch_length)?,
      ..Self::default()
    })
  }
}

impl EdgeToNwk for EdgeTimetree {
  fn nwk_weight(&self) -> Option<f64> {
    self.time_length
  }
}

impl EdgeToGraphviz for EdgeTimetree {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self
      .branch_length()
      .map(|weight| format_weight(weight, &NwkWriteOptions::default()))
  }

  fn to_graphviz_weight(&self) -> Option<f64> {
    self.branch_length()
  }
}

impl TimetreeEdge for EdgeTimetree {
  fn set_gamma(&mut self, gamma: f64) {
    self.gamma = gamma;
  }
}
