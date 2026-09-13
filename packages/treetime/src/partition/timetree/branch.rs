use crate::gtr::gtr::GTR;
use crate::partition::timetree::partition::PartitionTimetree;
use crate::partition::traits::{HasGtr, PartitionBranchOps};
use crate::seq::indel::InDel;
use crate::seq::mutation::Sub;
use eyre::Report;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;

impl HasGtr for PartitionTimetree {
  fn gtr(&self) -> &GTR {
    match self {
      Self::Dense(family) => family.partition.gtr(),
      Self::Sparse(family) => family.partition.gtr(),
    }
  }

  fn gtr_mut(&mut self) -> &mut GTR {
    match self {
      Self::Dense(family) => family.partition.gtr_mut(),
      Self::Sparse(family) => family.partition.gtr_mut(),
    }
  }

  fn sequence_length(&self) -> usize {
    match self {
      Self::Dense(family) => family.partition.length,
      Self::Sparse(family) => family.partition.length,
    }
  }
}

impl PartitionBranchOps for PartitionTimetree {
  fn sequence_length(&self) -> usize {
    match self {
      Self::Dense(family) => family.partition.length,
      Self::Sparse(family) => family.partition.length,
    }
  }

  fn edge_subs(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    match self {
      Self::Dense(family) => family.readout().edge_subs(graph, edge_key),
      Self::Sparse(family) => family.readout().edge_subs(graph, edge_key),
    }
  }

  fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<InDel> {
    match self {
      Self::Dense(family) => family.readout().edge_indels(edge_key),
      Self::Sparse(family) => family.readout().edge_indels(edge_key),
    }
  }

  fn root_sequence(&self, graph: &Graph) -> Result<Seq, Report> {
    match self {
      Self::Dense(family) => family.readout().root_sequence(graph),
      Self::Sparse(family) => family.readout().root_sequence(graph),
    }
  }

  fn node_sequence(&self, node_key: GraphNodeKey) -> Seq {
    match self {
      Self::Dense(family) => family.readout().node_sequence(node_key),
      Self::Sparse(family) => family.readout().node_sequence(node_key),
    }
  }

  fn edge_effective_length(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<usize, Report> {
    match self {
      Self::Dense(family) => family.readout().edge_effective_length(graph, edge_key),
      Self::Sparse(family) => family.readout().edge_effective_length(graph, edge_key),
    }
  }
}
