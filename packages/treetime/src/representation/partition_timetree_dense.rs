use crate::alphabet::alphabet::Alphabet;
use crate::distribution::distribution::Distribution;
use crate::graph::edge::GraphEdgeKey;
use crate::graph::graph::{GraphNodeBackward, GraphNodeForward};
use crate::graph::node::GraphNodeKey;
use crate::gtr::gtr::GTR;
use crate::io::fasta::FastaRecord;
use crate::representation::edge_timetree::EdgeTimetree;
use crate::representation::log_lh::HasLogLh;
use crate::representation::node_timetree::NodeTimetree;
use crate::representation::partition_marginal::{PartitionMarginal, PartitionMarginalOps};
use crate::representation::partition_timetree::{GraphTimetree, PartitionTimetree, PartitionTimetreeOps};
use crate::representation::seq::Seq;
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::sync::Arc;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NodeTimetreeDense {
  pub seq: super::graph_dense::DenseSeqInfo,
  pub profile: super::graph_dense::DenseSeqDis,
  pub msg_from_parent: Option<Arc<Distribution>>,
  pub msg_from_children: Option<Arc<Distribution>>,
}

#[derive(Clone, Debug, Serialize, Deserialize, Default)]
pub struct EdgeTimetreeDense {
  pub msg_to_parent: super::graph_dense::DenseSeqDis,
  pub msg_to_child: super::graph_dense::DenseSeqDis,
  pub msg_from_child: super::graph_dense::DenseSeqDis,
  pub time_msg_to_parent: Option<Arc<Distribution>>,
  pub time_msg_to_child: Option<Arc<Distribution>>,
  pub time_msg_from_child: Option<Arc<Distribution>>,
}

#[derive(Clone, Debug)]
pub struct PartitionTimetreeDense {
  pub index: usize,
  pub gtr: GTR,
  pub alphabet: Alphabet,
  pub sequence_length: Option<usize>,
  pub nodes: BTreeMap<GraphNodeKey, NodeTimetreeDense>,
  pub edges: BTreeMap<GraphEdgeKey, EdgeTimetreeDense>,
}

impl PartitionTimetree for PartitionTimetreeDense {}

impl PartitionTimetreeOps<NodeTimetree, EdgeTimetree> for PartitionTimetreeDense {
  fn initialize_nodes(&mut self, graph: &GraphTimetree) -> Result<(), Report> {
    for node_ref in graph.get_nodes() {
      let node_key = node_ref.read_arc().key();
      self.nodes.entry(node_key).or_insert_with(|| NodeTimetreeDense {
        msg_from_parent: None,
        msg_from_children: None,
      });
    }
    Ok(())
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

impl PartitionMarginalOps<NodeTimetree, EdgeTimetree> for PartitionTimetreeDense {
  fn attach_sequences(&mut self, _graph: &GraphTimetree, _aln: &[FastaRecord]) -> Result<(), Report> {
    Ok(())
  }

  fn process_node_backward(
    &mut self,
    _node: &GraphNodeBackward<NodeTimetree, EdgeTimetree, ()>,
  ) -> Result<(), Report> {
    Ok(())
  }

  fn process_node_forward(
    &mut self,
    _graph: &GraphTimetree,
    _node: &GraphNodeForward<NodeTimetree, EdgeTimetree, ()>,
  ) -> Result<(), Report> {
    Ok(())
  }

  fn extract_ancestral_sequence(&mut self, _node_key: GraphNodeKey) -> Seq {
    crate::seq![]
  }

  fn reconstruct_node_sequence(
    &mut self,
    _node: &GraphNodeForward<NodeTimetree, EdgeTimetree, ()>,
    _include_leaves: bool,
  ) -> Option<Seq> {
    None
  }

  fn get_sequence_length(&self) -> Option<usize> {
    self.sequence_length
  }
}
