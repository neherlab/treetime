use crate::alphabet::alphabet::Alphabet;
use crate::error::OperationError;
use crate::gtr::gtr::GTR;
use crate::partition::fitch::partition::PartitionFitch;
use crate::partition::marginal::reconstruction::MarginalReconstruction;
use crate::seq::indel::InDel;
use crate::seq::mutation::{MutationTrack, SequenceMutations, stream_sequence_mutations};
use crate::seq::sink::SeqSink;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;

#[expect(
  clippy::large_enum_variant,
  reason = "one value per run moves between steps by value; it is never stored in a collection"
)]
pub(crate) enum AncestralPartition {
  Fitch(PartitionFitch),
  Marginal {
    reconstruction: MarginalReconstruction,
    sampled: BTreeMap<GraphNodeKey, Seq>,
    impute: bool,
  },
}

impl AncestralPartition {
  pub(crate) fn gtr(&self) -> Option<&GTR> {
    match self {
      Self::Fitch(_) => None,
      Self::Marginal { reconstruction, .. } => Some(reconstruction.gtr()),
    }
  }

  pub(crate) fn alphabet(&self) -> &Alphabet {
    match self {
      Self::Fitch(partition) => &partition.alphabet,
      Self::Marginal { reconstruction, .. } => reconstruction.alphabet(),
    }
  }

  pub(crate) fn stream_sequences(
    &self,
    graph: &Graph,
    track: &MutationTrack,
    include_leaves: bool,
    sink: Option<&mut dyn SeqSink>,
  ) -> Result<SequenceMutations, OperationError> {
    stream_sequence_mutations(
      graph,
      self.alphabet(),
      track,
      include_leaves,
      |node_key| self.node_sequence(graph, node_key),
      |edge_key| self.edge_indels(edge_key),
      sink,
    )
  }

  fn node_sequence(&self, graph: &Graph, node_key: GraphNodeKey) -> Result<Seq, Report> {
    match self {
      Self::Fitch(partition) => Ok(partition.node_sequence(node_key)),
      Self::Marginal {
        reconstruction,
        sampled,
        impute,
      } => sampled
        .get(&node_key)
        .cloned()
        .map_or_else(|| reconstruction.node_sequence(graph, *impute, node_key), Ok),
    }
  }

  fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<InDel> {
    match self {
      Self::Fitch(partition) => partition.edge_indels(edge_key),
      Self::Marginal { reconstruction, .. } => reconstruction.edge_indels(edge_key),
    }
  }
}
