#[cfg(test)]
mod __tests__;

use crate::payload::ancestral::{EdgeAncestral, NodeAncestral};
use crate::payload::clock_set::ClockSet;
use crate::payload::traits::{ClockEdge, ClockNode, DateConstraintNode, TimetreeEdge, TimetreeNode};
use eyre::Report;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::{Distribution, NegLog};
use treetime_graph::edge::{BranchDistribution, ClockMessages, GraphEdge, HasBranchLength, TimeLength};
use treetime_graph::node::{Described, Divergence, GraphNode, Named, Outlier, TimeConstraint};
use treetime_io::graphviz::{EdgeToGraphviz, NodeToGraphviz};
use treetime_io::nwk::{EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk};

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct NodeTimetree {
  pub base: NodeAncestral,
  pub time: Option<f64>,
  /// Distribution of the date given as input, `None` for a node no date was given for. Written
  /// once by [`load_date_constraints`](crate::clock::date_constraints::load_date_constraints) and
  /// then held fixed: `time_distribution` is refined in place by the inference passes, and every
  /// backward pass lifts this back into it so the input is never inferred away.
  pub date_constraint: Option<Arc<Distribution<NegLog>>>,
  pub time_distribution: Option<Arc<Distribution<NegLog>>>,
  pub bad_branch: bool,
  pub div: f64,
  pub is_outlier: bool,
  pub clock_set: ClockSet,
  /// Node dates inferred at three clock rates [lower, central, upper], sorted by date.
  /// Populated by rate susceptibility analysis (--confidence with --covariation or --clock-std-dev).
  /// See Sagulenko, Puller & Neher 2018, Section 2.2 (marginal date inference and confidence intervals).
  pub rate_susceptibility_dates: Option<[f64; 3]>,
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
  /// The observed date where there is one, the current estimate otherwise.
  ///
  /// The clock regression and the clock filter read this for leaves, and must see the date as
  /// given: the forward pass refines the time distribution of a leaf whose date is uncertain, and
  /// regressing on that refined date would feed the tree's own inference back into the clock.
  fn likely_time(&self) -> Option<f64> {
    self
      .date_constraint
      .as_ref()
      .or(self.time_distribution.as_ref())
      .and_then(|dist| dist.likely_time())
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

impl TimeConstraint<Arc<Distribution<NegLog>>> for NodeTimetree {
  fn date_constraint(&self) -> &Option<Arc<Distribution<NegLog>>> {
    &self.date_constraint
  }

  fn set_date_constraint(&mut self, dist: Option<Arc<Distribution<NegLog>>>) {
    self.date_constraint = dist;
  }

  fn time_distribution(&self) -> &Option<Arc<Distribution<NegLog>>> {
    &self.time_distribution
  }

  fn set_time_distribution(&mut self, dist: Option<Arc<Distribution<NegLog>>>) {
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
  fn from_nwk(
    name: Option<impl AsRef<str>>,
    confidence: Option<f64>,
    comments: &BTreeMap<String, String>,
  ) -> Result<Self, Report> {
    Ok(Self {
      base: NodeAncestral::from_nwk(name, confidence, comments)?,
      ..NodeTimetree::default()
    })
  }
}

impl NodeToNwk for NodeTimetree {
  fn nwk_name(&self) -> Option<impl AsRef<str>> {
    self.base.nwk_name()
  }

  fn nwk_comments(&self) -> BTreeMap<String, String> {
    let mut comments = self.base.nwk_comments();
    if let Some(time) = self.time {
      // 2-decimal precision matches v0 Newick format. See kb/decisions/timetree-nwk-date-two-decimal-precision.md
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
  /// Clock-constrained branch length in substitutions per site, `clock_rate * gamma * dt` over
  /// the inferred node times. Distinct from `base.branch_length`, which stays the ML or input
  /// length the branch-length grid is centred on and the root-to-tip regression reads; this one
  /// is what sequence profiles propagate along. v0 keeps the same two-length split as
  /// `branch_length` and `mutation_length`.
  pub clock_branch_length: Option<f64>,
  pub branch_length_distribution: Option<Arc<Distribution<NegLog>>>,
  pub msg_to_parent: Option<Arc<Distribution<NegLog>>>,
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

  /// The clock-constrained length once it has been committed, the ML or input length before
  /// that and after a topology change has invalidated it.
  fn profile_branch_length(&self) -> Option<f64> {
    self.clock_branch_length.or_else(|| self.branch_length())
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

impl BranchDistribution<Arc<Distribution<NegLog>>> for EdgeTimetree {
  fn branch_length_distribution(&self) -> &Option<Arc<Distribution<NegLog>>> {
    &self.branch_length_distribution
  }

  fn set_branch_length_distribution(&mut self, dist: Option<Arc<Distribution<NegLog>>>) {
    self.branch_length_distribution = dist;
  }

  fn msg_to_parent(&self) -> &Option<Arc<Distribution<NegLog>>> {
    &self.msg_to_parent
  }

  fn set_msg_to_parent(&mut self, msg: Option<Arc<Distribution<NegLog>>>) {
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

impl EdgeToGraphviz for EdgeTimetree {}

impl TimetreeEdge for EdgeTimetree {
  fn set_gamma(&mut self, gamma: f64) {
    self.gamma = gamma;
  }

  fn clock_branch_length(&self) -> Option<f64> {
    self.clock_branch_length
  }

  fn set_clock_branch_length(&mut self, length: Option<f64>) {
    self.clock_branch_length = length;
  }
}
