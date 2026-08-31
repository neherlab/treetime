use crate::gtr::gtr::GTR;
use crate::partition::timetree::partition::PartitionTimetree;
use crate::partition::traits::{BranchTopology, HasGtr, PartitionBranchOps};
use crate::seq::indel::InDel;
use crate::seq::mutation::Sub;
use eyre::Report;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;

impl HasGtr for PartitionTimetree {
  fn gtr(&self) -> &GTR {
    match self {
      Self::Dense(partition) => partition.gtr(),
      Self::Sparse(partition) => partition.gtr(),
    }
  }

  fn gtr_mut(&mut self) -> &mut GTR {
    match self {
      Self::Dense(partition) => partition.gtr_mut(),
      Self::Sparse(partition) => partition.gtr_mut(),
    }
  }

  fn sequence_length(&self) -> usize {
    match self {
      Self::Dense(partition) => HasGtr::sequence_length(partition),
      Self::Sparse(partition) => HasGtr::sequence_length(partition),
    }
  }
}

impl PartitionBranchOps for PartitionTimetree {
  fn sequence_length(&self) -> usize {
    match self {
      Self::Dense(partition) => PartitionBranchOps::sequence_length(partition),
      Self::Sparse(partition) => PartitionBranchOps::sequence_length(partition),
    }
  }

  fn edge_subs(&self, graph: &dyn BranchTopology, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    match self {
      Self::Dense(partition) => partition.edge_subs(graph, edge_key),
      Self::Sparse(partition) => partition.edge_subs(graph, edge_key),
    }
  }

  fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<InDel> {
    match self {
      Self::Dense(partition) => partition.edge_indels(edge_key),
      Self::Sparse(partition) => partition.edge_indels(edge_key),
    }
  }

  fn root_sequence(&self, graph: &dyn BranchTopology) -> Result<Seq, Report> {
    match self {
      Self::Dense(partition) => partition.root_sequence(graph),
      Self::Sparse(partition) => partition.root_sequence(graph),
    }
  }

  fn node_sequence(&self, node_key: GraphNodeKey) -> Seq {
    match self {
      Self::Dense(partition) => partition.node_sequence(node_key),
      Self::Sparse(partition) => partition.node_sequence(node_key),
    }
  }

  fn edge_effective_length(&self, graph: &dyn BranchTopology, edge_key: GraphEdgeKey) -> Result<usize, Report> {
    match self {
      Self::Dense(partition) => partition.edge_effective_length(graph, edge_key),
      Self::Sparse(partition) => partition.edge_effective_length(graph, edge_key),
    }
  }
}
