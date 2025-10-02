use crate::alphabet::alphabet::Alphabet;
use crate::commands::timetree::date_constraints::DateConstraint;
use crate::distribution::distribution::Distribution;
use crate::graph::edge::GraphEdgeKey;
use crate::graph::graph::{GraphNodeBackward, GraphNodeForward};
use crate::graph::node::GraphNodeKey;
use crate::gtr::gtr::GTR;
use crate::io::fasta::FastaRecord;
use crate::representation::graph_ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use crate::representation::log_lh::HasLogLh;
use crate::representation::partition_marginal::{PartitionMarginal, PartitionMarginalOps};
use crate::representation::partition_timetree::{PartitionTimetree, PartitionTimetreeOps};
use crate::representation::seq::Seq;
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::sync::Arc;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NodeTimetreeDense {
  pub constraint: Option<DateConstraint>,
  pub bad_branch: bool,
  pub time_before_present: Option<f64>,
  pub posterior: Option<Arc<Distribution>>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EdgeTimetreeDense {
  pub clock_length: Option<f64>,
  pub msg_to_parent: Option<Arc<Distribution>>,
  pub msg_to_child: Option<Arc<Distribution>>,
  pub msg_from_child: Option<Arc<Distribution>>,
}

#[derive(Clone, Debug)]
pub struct PartitionTimetreeDense {
  pub index: usize,
  pub gtr: GTR,
  pub alphabet: Alphabet,
  pub length: usize,
  pub sequence_length: Option<usize>,
  pub nodes: BTreeMap<GraphNodeKey, NodeTimetreeDense>,
  pub edges: BTreeMap<GraphEdgeKey, EdgeTimetreeDense>,
}

impl PartitionTimetree for PartitionTimetreeDense {}

impl PartitionTimetreeOps for PartitionTimetreeDense {
  fn attach_date_constraints(
    &mut self,
    _graph: &GraphAncestral,
    constraints: &BTreeMap<GraphNodeKey, DateConstraint>,
  ) -> Result<(), Report> {
    for (node_key, constraint) in constraints {
      self
        .nodes
        .entry(*node_key)
        .or_insert_with(|| NodeTimetreeDense {
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
    todo!("Implement backward pass for dense timetree partition")
  }

  fn process_node_forward(
    &mut self,
    _graph: &GraphAncestral,
    _node: &GraphNodeForward<NodeAncestral, EdgeAncestral, ()>,
  ) -> Result<(), Report> {
    todo!("Implement forward pass for dense timetree partition")
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

impl HasLogLh for PartitionTimetreeDense {
  fn get_log_lh(&self, _node_key: GraphNodeKey) -> f64 {
    0.0
  }
}

impl PartitionMarginal for PartitionTimetreeDense {}

impl PartitionMarginalOps for PartitionTimetreeDense {
  fn attach_sequences(&mut self, _graph: &GraphAncestral, _aln: &[FastaRecord]) -> Result<(), Report> {
    Ok(())
  }

  fn process_node_backward(
    &mut self,
    node: &GraphNodeBackward<NodeAncestral, EdgeAncestral, ()>,
  ) -> Result<(), Report> {
    PartitionTimetreeOps::process_node_backward(self, node)
  }

  fn process_node_forward(
    &mut self,
    graph: &GraphAncestral,
    node: &GraphNodeForward<NodeAncestral, EdgeAncestral, ()>,
  ) -> Result<(), Report> {
    PartitionTimetreeOps::process_node_forward(self, graph, node)
  }

  fn extract_ancestral_sequence(&mut self, _node_key: GraphNodeKey) -> Seq {
    crate::seq![]
  }

  fn reconstruct_node_sequence(
    &mut self,
    _node: &GraphNodeForward<NodeAncestral, EdgeAncestral, ()>,
    _include_leaves: bool,
  ) -> Option<Seq> {
    None
  }

  fn get_length(&self) -> usize {
    self.length
  }
}
