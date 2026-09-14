use crate::gtr::gtr::GTR;
use crate::partition::timetree::partition::PartitionTimetree;
use crate::partition::traits::HasGtr;
use crate::seq::indel::InDel;
use crate::seq::mutation::{Mutation, MutationTrack, Sub, combine_edge_mutations};
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

/// Output-side read access over a completed timetree reconstruction, dispatching each operation to the
/// concrete representation. The timetree tree writers read `node_sequence`/`root_sequence`/
/// `edge_mutations`, and the divergence reader reads `edge_subs`.
impl PartitionTimetree {
  /// MAP-derived nucleotide substitutions on one edge (parent -> child).
  pub fn edge_subs(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    match self {
      Self::Dense(family) => family.edge_subs(graph, edge_key),
      Self::Sparse(family) => family.edge_subs(edge_key),
    }
  }

  /// Grouped aligned insertions and deletions for one edge.
  pub fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<InDel> {
    match self {
      Self::Dense(family) => family.edge_indels(edge_key),
      Self::Sparse(family) => family.edge_indels(edge_key),
    }
  }

  /// The reconstructed root sequence.
  pub fn root_sequence(&self, graph: &Graph) -> Result<Seq, Report> {
    match self {
      Self::Dense(family) => family.root_sequence(graph),
      Self::Sparse(family) => family.root_sequence(graph),
    }
  }

  /// The reconstructed sequence for one node.
  pub fn node_sequence(&self, node_key: GraphNodeKey) -> Seq {
    match self {
      Self::Dense(family) => family.node_sequence(node_key),
      Self::Sparse(family) => family.node_sequence(node_key),
    }
  }

  /// The substitutions and indels on one edge as one mutation list on the given track.
  pub fn edge_mutations(
    &self,
    graph: &Graph,
    edge_key: GraphEdgeKey,
    track: MutationTrack,
  ) -> Result<Vec<Mutation>, Report> {
    combine_edge_mutations(self.edge_subs(graph, edge_key)?, &self.edge_indels(edge_key), &track)
  }
}
