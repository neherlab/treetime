use crate::alphabet::alphabet::Alphabet;
use crate::commands::timetree::date_constraints::DateConstraint;
use crate::distribution::distribution::Distribution;
use crate::graph::edge::GraphEdgeKey;
use crate::graph::graph::{GraphNodeBackward, GraphNodeForward};
use crate::graph::node::GraphNodeKey;
use crate::gtr::gtr::GTR;
use crate::representation::graph_ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use crate::representation::partition_timetree::{PartitionTimetree, PartitionTimetreeOps};
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::sync::Arc;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NodeTimetreeSparse {
  pub constraint: Option<DateConstraint>,
  pub bad_branch: bool,
  pub time_before_present: Option<f64>,
  pub posterior: Option<Arc<Distribution>>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EdgeTimetreeSparse {
  pub clock_length: Option<f64>,
  pub msg_to_parent: Option<Arc<Distribution>>,
  pub msg_to_child: Option<Arc<Distribution>>,
  pub msg_from_child: Option<Arc<Distribution>>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PartitionTimetreeSparse {
  pub index: usize,
  pub gtr: GTR,
  pub alphabet: Alphabet,
  pub length: usize,
  pub sequence_length: Option<usize>,
  pub nodes: BTreeMap<GraphNodeKey, NodeTimetreeSparse>,
  pub edges: BTreeMap<GraphEdgeKey, EdgeTimetreeSparse>,
}

impl PartitionTimetree for PartitionTimetreeSparse {}

impl PartitionTimetreeOps for PartitionTimetreeSparse {
  fn attach_date_constraints(
    &mut self,
    _graph: &GraphAncestral,
    constraints: &BTreeMap<GraphNodeKey, DateConstraint>,
  ) -> Result<(), Report> {
    for (node_key, constraint) in constraints {
      self
        .nodes
        .entry(*node_key)
        .or_insert_with(|| NodeTimetreeSparse {
          constraint: None,
          bad_branch: false,
          time_before_present: None,
          posterior: None,
        })
        .constraint = Some(constraint.clone());
    }
    Ok(())
  }

  fn process_node_backward(
    &mut self,
    _node: &GraphNodeBackward<NodeAncestral, EdgeAncestral, ()>,
  ) -> Result<(), Report> {
    todo!("Implement backward pass for sparse timetree partition")
  }

  fn process_node_forward(
    &mut self,
    _graph: &GraphAncestral,
    _node: &GraphNodeForward<NodeAncestral, EdgeAncestral, ()>,
  ) -> Result<(), Report> {
    todo!("Implement forward pass for sparse timetree partition")
  }

  fn extract_node_time(&self, node_key: GraphNodeKey) -> Option<f64> {
    self.nodes.get(&node_key)?.time_before_present
  }

  fn extract_node_posterior(&self, node_key: GraphNodeKey) -> Option<Arc<Distribution>> {
    self.nodes.get(&node_key)?.posterior.clone()
  }

  fn get_sequence_length(&self) -> Option<usize> {
    self.sequence_length
  }
}
