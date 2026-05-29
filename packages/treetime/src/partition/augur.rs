use crate::partition::fitch::PartitionFitch;
use crate::partition::marginal_dense::PartitionMarginalDense;
use crate::partition::marginal_sparse::PartitionMarginalSparse;
use crate::partition::traits::{BranchTopology, PartitionBranchOps, PartitionMarginalOps};
use crate::payload::ancestral::{EdgeAncestral, NodeAncestral};
use crate::seq::mutation::Sub;
use eyre::Report;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{AsciiChar, Seq};

/// Partition-side data access for the `ancestral` augur node data JSON writer.
///
/// The writer needs, per node, a reconstructed sequence and the substitutions on
/// the node's parent edge, plus the alignment length. These are available from
/// every partition representation that carries a reconstruction (Fitch parsimony,
/// marginal sparse, marginal dense), but via different concrete data paths. This
/// trait erases those differences so a single writer serves all of them.
///
/// `root_sequence` resolves the root via the graph and returns its reconstructed
/// sequence, guaranteeing the JSON `reference.nuc` equals the root node's
/// `sequence` field (augur derives both from the reconstructed root).
pub trait AugurNodeDataJsonAncestralPartition {
  /// Alignment length (number of sites).
  fn sequence_length(&self) -> usize;

  /// Reconstructed sequence for one node.
  fn node_sequence(&self, node_key: GraphNodeKey) -> Seq;

  /// Substitutions on the parent edge of one node (parent -> child).
  fn edge_subs(&self, graph: &dyn BranchTopology, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report>;

  /// Ambiguous (unknown) character of the partition alphabet, used to fill
  /// masked positions in per-node output sequences.
  fn ambiguous_char(&self) -> AsciiChar;

  /// Reconstructed root sequence, used as the JSON reference.
  fn root_sequence(&self, graph: &dyn BranchTopology) -> Result<Seq, Report> {
    Ok(self.node_sequence(graph.root_key()?))
  }
}

impl AugurNodeDataJsonAncestralPartition for PartitionFitch {
  fn sequence_length(&self) -> usize {
    self.length
  }

  fn node_sequence(&self, node_key: GraphNodeKey) -> Seq {
    self.nodes[&node_key].seq.sequence.clone()
  }

  fn edge_subs(&self, _graph: &dyn BranchTopology, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    Ok(self.edges[&edge_key].fitch_subs().to_vec())
  }

  fn ambiguous_char(&self) -> AsciiChar {
    self.alphabet.unknown()
  }
}

impl AugurNodeDataJsonAncestralPartition for PartitionMarginalSparse {
  fn sequence_length(&self) -> usize {
    self.length
  }

  fn node_sequence(&self, node_key: GraphNodeKey) -> Seq {
    self.nodes[&node_key].seq.sequence.clone()
  }

  fn edge_subs(&self, graph: &dyn BranchTopology, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    PartitionBranchOps::edge_subs(self, graph, edge_key)
  }

  fn ambiguous_char(&self) -> AsciiChar {
    self.alphabet.unknown()
  }
}

impl AugurNodeDataJsonAncestralPartition for PartitionMarginalDense {
  fn sequence_length(&self) -> usize {
    self.length
  }

  fn node_sequence(&self, node_key: GraphNodeKey) -> Seq {
    // Dense partitions do not store full sequences; they are derived from the
    // per-site marginal posteriors (MAP states) via the marginal ops trait.
    <Self as PartitionMarginalOps<NodeAncestral, EdgeAncestral>>::extract_ancestral_sequence(self, node_key)
  }

  fn edge_subs(&self, graph: &dyn BranchTopology, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    PartitionBranchOps::edge_subs(self, graph, edge_key)
  }

  fn ambiguous_char(&self) -> AsciiChar {
    self.alphabet.unknown()
  }
}
