use crate::ancestral::sample::SampleMode;
use crate::gtr::gtr::GTR;
use crate::partition::marginal_dense::PartitionMarginalDense;
use crate::partition::marginal_sparse::PartitionMarginalSparse;
use crate::partition::optimization_contribution::OptimizationContribution;
use crate::partition::traits::{
  BranchTopology, HasGtr, HasLogLh, PartitionBranchOps, PartitionMarginalOps, PartitionMarginalPasses,
  PartitionOptimizeOps, PartitionRerootOps, PartitionTimetreeOps,
};
use crate::payload::timetree::{EdgeTimetree, NodeTimetree};
use crate::seq::indel::InDel;
use crate::seq::mutation::Sub;
use eyre::Report;
use parking_lot::RwLock;
use serde::Serialize;
use std::sync::Arc;
use treetime_graph::edge::{EdgeOptimizeOps, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::GraphNodeForward;
use treetime_graph::node::{GraphNode, GraphNodeKey, Named};
use treetime_graph::reroot::RerootChanges;
use treetime_io::fasta::FastaRecord;
use treetime_primitives::Seq;

pub type GraphTimetree<D = ()> = Graph<NodeTimetree, EdgeTimetree, D>;
pub type PartitionTimetreeRef = Arc<RwLock<PartitionTimetree>>;
pub type PartitionTimetreeAllVec = Vec<PartitionTimetreeRef>;

#[derive(Debug, Serialize)]
pub enum PartitionTimetree {
  Dense(PartitionMarginalDense),
  Sparse(PartitionMarginalSparse),
}

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

impl HasLogLh for PartitionTimetree {
  fn get_log_lh(&self, node_key: GraphNodeKey) -> f64 {
    match self {
      Self::Dense(partition) => partition.get_log_lh(node_key),
      Self::Sparse(partition) => partition.get_log_lh(node_key),
    }
  }

  fn reset_node_log_likelihoods(&mut self) {
    match self {
      Self::Dense(partition) => partition.reset_node_log_likelihoods(),
      Self::Sparse(partition) => partition.reset_node_log_likelihoods(),
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

impl PartitionOptimizeOps for PartitionTimetree {
  fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report> {
    match self {
      Self::Dense(partition) => partition.create_edge_contribution(edge_key),
      Self::Sparse(partition) => partition.create_edge_contribution(edge_key),
    }
  }

  fn edge_indel_count(&self, edge_key: GraphEdgeKey) -> usize {
    match self {
      Self::Dense(partition) => partition.edge_indel_count(edge_key),
      Self::Sparse(partition) => partition.edge_indel_count(edge_key),
    }
  }
}

impl PartitionRerootOps for PartitionTimetree {
  fn apply_reroot(&mut self, changes: &RerootChanges) -> Result<(), Report> {
    match self {
      Self::Dense(partition) => partition.apply_reroot(changes),
      Self::Sparse(partition) => partition.apply_reroot(changes),
    }
  }
}

impl<N, E> PartitionTimetreeOps<N, E> for PartitionTimetree
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  fn reconcile_topology(&mut self, graph: &Graph<N, E, ()>) {
    match self {
      Self::Dense(partition) => partition.reconcile_topology(graph),
      Self::Sparse(partition) => partition.reconcile_topology(graph),
    }
  }
}

impl PartitionMarginalPasses<NodeTimetree, EdgeTimetree> for PartitionTimetree {
  fn process_backward_pass(&mut self, graph: &GraphTimetree) -> Result<(), Report> {
    match self {
      Self::Dense(partition) => partition.process_backward_pass(graph),
      Self::Sparse(partition) => partition.process_backward_pass(graph),
    }
  }

  fn process_forward_pass(&mut self, graph: &GraphTimetree) -> Result<(), Report> {
    match self {
      Self::Dense(partition) => partition.process_forward_pass(graph),
      Self::Sparse(partition) => partition.process_forward_pass(graph),
    }
  }

  fn get_sequence_length(&self) -> usize {
    match self {
      Self::Dense(partition) => partition.get_sequence_length(),
      Self::Sparse(partition) => partition.get_sequence_length(),
    }
  }
}

impl PartitionMarginalOps<NodeTimetree, EdgeTimetree> for PartitionTimetree {
  fn attach_sequences(&mut self, graph: &GraphTimetree, aln: &[FastaRecord]) -> Result<(), Report> {
    match self {
      Self::Dense(partition) => partition.attach_sequences(graph, aln),
      Self::Sparse(partition) => partition.attach_sequences(graph, aln),
    }
  }

  fn extract_ancestral_sequence(&self, node_key: GraphNodeKey) -> Seq {
    match self {
      Self::Dense(partition) => {
        <PartitionMarginalDense as PartitionMarginalOps<NodeTimetree, EdgeTimetree>>::extract_ancestral_sequence(
          partition, node_key,
        )
      },
      Self::Sparse(partition) => {
        <PartitionMarginalSparse as PartitionMarginalOps<NodeTimetree, EdgeTimetree>>::extract_ancestral_sequence(
          partition, node_key,
        )
      },
    }
  }

  fn reconstruct_node_sequence(
    &mut self,
    node: &GraphNodeForward<NodeTimetree, EdgeTimetree>,
    include_leaves: bool,
    sample_mode: SampleMode,
    rng: &mut dyn rand::RngCore,
  ) -> Option<Seq> {
    match self {
      Self::Dense(partition) => partition.reconstruct_node_sequence(node, include_leaves, sample_mode, rng),
      Self::Sparse(partition) => partition.reconstruct_node_sequence(node, include_leaves, sample_mode, rng),
    }
  }
}
