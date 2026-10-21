use crate::alphabet::alphabet::Alphabet;
use crate::error::OperationError;
use crate::gtr::gtr::GTR;
use crate::partition::fitch::partition::PartitionFitch;
use crate::partition::marginal::reconstruction::MarginalReconstruction;
use crate::seq::indel::InDel;
use crate::seq::mutation::{SequenceMutations, stream_sequence_mutations};
use crate::seq::sink::{SeqStream, SeqTrack};
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{LogLh, Seq};

#[expect(
  clippy::large_enum_variant,
  reason = "a run holds one value per alignment file, a handful at most, so boxing the large variant saves no meaningful memory"
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

  pub(crate) fn log_lh(&self, graph: &Graph) -> Result<Option<LogLh>, Report> {
    match self {
      Self::Fitch(_) => Ok(None),
      Self::Marginal { reconstruction, .. } => reconstruction.graph_log_lh(graph).map(Some),
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
    track: SeqTrack<'_>,
    include_leaves: bool,
    report_unknown: bool,
    sink: Option<SeqStream<'_>>,
  ) -> Result<SequenceMutations, OperationError> {
    if let Self::Marginal {
      reconstruction,
      sampled,
      impute,
    } = self
      && sampled.is_empty()
    {
      return reconstruction.stream_sequences(graph, *impute, track, include_leaves, report_unknown, sink);
    }
    stream_sequence_mutations(
      graph,
      self.alphabet(),
      track,
      include_leaves,
      report_unknown,
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
