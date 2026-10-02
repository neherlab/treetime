use crate::gtr::gtr::GTR;
use crate::partition::fitch::partition::PartitionFitch;
use crate::partition::marginal::reconstruction::MarginalReconstruction;
use crate::seq::indel::InDel;
use crate::seq::mutation::{Mutation, MutationTrack, Sub, combine_edge_mutations};
use eyre::Report;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{AsciiChar, Seq};

#[expect(
  clippy::large_enum_variant,
  reason = "one value per run moves between steps by value; it is never stored in a collection"
)]
#[derive(Clone, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum AncestralPartition {
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

  pub fn sequence_length(&self) -> usize {
    match self {
      Self::Fitch(partition) => partition.sequence_length(),
      Self::Marginal { reconstruction, .. } => reconstruction.sequence_length(),
    }
  }

  pub fn ambiguous_char(&self) -> AsciiChar {
    match self {
      Self::Fitch(partition) => partition.ambiguous_char(),
      Self::Marginal { reconstruction, .. } => reconstruction.ambiguous_char(),
    }
  }

  pub fn augur_node_sequence(&self, graph: &Graph, node_key: GraphNodeKey) -> Result<Seq, Report> {
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

  pub fn root_sequence(&self, graph: &Graph) -> Result<Seq, Report> {
    match self {
      Self::Fitch(partition) => partition.root_sequence(graph),
      Self::Marginal { reconstruction, .. } => reconstruction.root_sequence(graph),
    }
  }

  pub fn augur_root_sequence(&self, graph: &Graph) -> Result<Seq, Report> {
    self.augur_node_sequence(graph, graph.root_key()?)
  }

  pub fn edge_subs(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    match self {
      Self::Fitch(partition) => partition.edge_subs(edge_key),
      Self::Marginal { reconstruction, .. } => reconstruction.edge_subs(graph, edge_key),
    }
  }

  pub(crate) fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<InDel> {
    match self {
      Self::Fitch(partition) => partition.edge_indels(edge_key),
      Self::Marginal { reconstruction, .. } => reconstruction.edge_indels(edge_key),
    }
  }

  pub fn edge_mutations(
    &self,
    graph: &Graph,
    edge_key: GraphEdgeKey,
    track: &MutationTrack,
  ) -> Result<Vec<Mutation>, Report> {
    combine_edge_mutations(self.edge_subs(graph, edge_key)?, &self.edge_indels(edge_key), track)
  }
}
