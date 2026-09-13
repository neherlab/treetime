use crate::gtr::gtr::GTR;
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::seq::indel::InDel;
use crate::seq::mutation::{Mutation, MutationTrack, Sub};
use eyre::Report;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;

/// Access to a partition's substitution model and sequence length, shared by the concrete
/// representations so rate normalization can stay generic over them.
pub trait HasGtr {
  fn gtr(&self) -> &GTR;
  fn gtr_mut(&mut self) -> &mut GTR;
  fn sequence_length(&self) -> usize;

  fn weighted_rate(&self) -> f64 {
    self.sequence_length() as f64 * self.gtr().mu
  }

  fn normalize_rate(&mut self, scale: f64) {
    self.gtr_mut().mu /= scale;
  }
}

/// Read accessors over a completed marginal reconstruction, shared by dense and sparse representations.
///
/// Implemented by short-lived per-representation read views that borrow the durable partition inputs
/// together with the node states and edge messages/estimates the passes returned. The view is a
/// transient read projection assembled at a consumer boundary, never a stored stage-filled object.
///
/// Requires marginal inference to have run. Pre-marginal consumers (GTR inference, prune/merge) read
/// Fitch data directly.
pub trait PartitionBranchOps: Send + Sync {
  /// Return the sequence length represented by this partition.
  fn sequence_length(&self) -> usize;

  /// Return MAP-derived nucleotide substitutions for one edge.
  fn edge_subs(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report>;

  /// Return grouped aligned insertions and deletions for one edge.
  fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<InDel>;

  /// Return the reconstructed root sequence represented by this partition.
  fn root_sequence(&self, graph: &Graph) -> Result<Seq, Report>;

  /// Return the reconstructed sequence for one node.
  fn node_sequence(&self, node_key: GraphNodeKey) -> Seq;

  fn edge_mutations(
    &self,
    graph: &Graph,
    edge_key: GraphEdgeKey,
    track: MutationTrack,
  ) -> Result<Vec<Mutation>, Report> {
    self
      .edge_subs(graph, edge_key)?
      .into_iter()
      .map(|substitution| Ok(Mutation::substitution(track.clone(), substitution)))
      .chain(
        self
          .edge_indels(edge_key)
          .iter()
          .map(|indel| Mutation::indel(track.clone(), indel)),
      )
      .collect()
  }

  /// Return the number of alignment positions where both parent and child
  /// have canonical (non-gap, non-ambiguous) states for one edge.
  fn edge_effective_length(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<usize, Report>;
}

/// Optimize-specific read accessors, extending [`PartitionBranchOps`] with the per-edge likelihood
/// contribution and indel count the branch-length optimizer needs. Implemented by the same short-lived
/// per-representation read views.
pub trait PartitionOptimizeOps: PartitionBranchOps {
  /// Return the precomputed likelihood contribution for one edge.
  fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report>;

  /// Return the number of indel events on one edge.
  fn edge_indel_count(&self, edge_key: GraphEdgeKey) -> usize;
}
