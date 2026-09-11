use crate::payload::clock_set::ClockSet;
use crate::payload::traits::ClockEdge;
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use treetime_graph::edge::{ClockMessages, GraphEdge, HasBranchLength, TimeLength};
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
pub struct EdgeClock {
  pub branch_length: Option<f64>,
  pub clock_to_parent: ClockSet,
  pub clock_to_child: ClockSet,
  pub clock_from_child: ClockSet, // this is the propagated 'to_parent' msg. only need to avoid recalculation of propagated message
}

impl GraphEdge for EdgeClock {}

impl HasBranchLength for EdgeClock {
  fn branch_length(&self) -> Option<f64> {
    self.branch_length
  }

  fn set_branch_length(&mut self, weight: Option<f64>) {
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
    self.branch_length()
  }
}

impl EdgeToGraphviz for EdgeClock {}

impl ClockMessages<ClockSet> for EdgeClock {
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

impl TimeLength for EdgeClock {
  fn time_length(&self) -> Option<f64> {
    None
  }

  fn set_time_length(&mut self, _length: Option<f64>) {}
}

impl ClockEdge for EdgeClock {}
