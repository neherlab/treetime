use crate::commands::clock::clock_set::ClockSet;
use crate::commands::clock::clock_traits::ClockEdge;
use crate::commands::timetree::timetree_traits::TimetreeEdge;
use crate::distribution::distribution::Distribution;
use crate::graph::edge::{GraphEdge, NumMuts, Weighted};
use crate::io::graphviz::EdgeToGraphViz;
use crate::io::nwk::{EdgeFromNwk, EdgeToNwk, NwkWriteOptions, format_weight};
use crate::representation::graph_ancestral::{EdgeAncestral, SparseSeqEdge};
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct EdgeTimetree {
  pub sparse_partitions: Vec<SparseSeqEdge>,
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

impl GraphEdge for EdgeTimetree {}

impl Weighted for EdgeTimetree {
  fn weight(&self) -> Option<f64> {
    self.branch_length
  }

  fn set_weight(&mut self, weight: Option<f64>) {
    self.branch_length = weight;
  }
}

impl NumMuts for EdgeTimetree {
  fn num_muts(&self) -> Option<usize> {
    if self.sparse_partitions.is_empty() {
      None
    } else {
      Some(self.sparse_partitions.iter().map(|p| p.subs.len()).sum())
    }
  }
}

impl ClockEdge for EdgeTimetree {
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

impl EdgeFromNwk for EdgeTimetree {
  fn from_nwk(branch_length: Option<f64>) -> Result<Self, Report> {
    Ok(Self {
      branch_length,
      ..Self::default()
    })
  }
}

impl EdgeToNwk for EdgeTimetree {
  fn nwk_weight(&self) -> Option<f64> {
    self.weight()
  }
}

impl EdgeToGraphViz for EdgeTimetree {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self
      .weight()
      .map(|weight| format_weight(weight, &NwkWriteOptions::default()))
  }

  fn to_graphviz_weight(&self) -> Option<f64> {
    self.weight()
  }
}

impl TimetreeEdge for EdgeTimetree {
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

impl From<&EdgeAncestral> for EdgeTimetree {
  fn from(edge: &EdgeAncestral) -> Self {
    Self {
      sparse_partitions: edge.sparse_partitions.clone(),
      branch_length: edge.branch_length,
      branch_length_distribution: edge.branch_length_distribution.clone(),
      msg_to_parent: edge.msg_to_parent.clone(),
      clock_to_parent: edge.clock_to_parent.clone(),
      clock_to_child: edge.clock_to_child.clone(),
      clock_from_child: edge.clock_from_child.clone(),
    }
  }
}
