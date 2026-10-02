use crate::partition::marginal::shared::update::MarginalPasses;
use crate::partition::timetree::partition::PartitionTimetree;
use crate::seq::indel::InDel;
use crate::seq::mutation::{Mutation, MutationTrack, Sub, combine_edge_mutations};
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{LogLh, Seq};

impl PartitionTimetree {
  pub(crate) fn graph_log_lh(&self, graph: &Graph) -> Result<LogLh, Report> {
    let root_key = graph.get_exactly_one_root()?.key();
    Ok(match self {
      Self::Dense(family) => family.partition.get_log_lh(&family.node_states, root_key),
      Self::Sparse(family) => family.partition.get_log_lh(&family.node_states, root_key),
    })
  }

  pub(crate) fn marginal_update(
    self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  ) -> Result<Self, Report> {
    Ok(match self {
      Self::Dense(family) => Self::Dense(family.marginal_update(graph, branch_lengths)?.0),
      Self::Sparse(family) => Self::Sparse(family.marginal_update(graph, branch_lengths)?.0),
    })
  }

  pub(crate) fn extract_ancestral_sequence(&self, node_key: GraphNodeKey) -> Seq {
    match self {
      Self::Dense(family) => family
        .partition
        .extract_ancestral_sequence(&family.node_states, node_key),
      Self::Sparse(family) => family
        .partition
        .extract_ancestral_sequence(&family.node_states, node_key),
    }
  }

  pub(crate) fn node_sequence(&self, graph: &Graph, impute: bool, node_key: GraphNodeKey) -> Result<Seq, Report> {
    match self {
      Self::Dense(family) => family.node_sequence(graph, impute, node_key),
      Self::Sparse(family) => family.node_sequence(graph, impute, node_key),
    }
  }

  pub fn edge_subs(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    match self {
      Self::Dense(family) => family.edge_subs(graph, edge_key),
      Self::Sparse(family) => family.edge_subs(edge_key),
    }
  }

  pub(crate) fn edge_sub_count(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<Option<usize>, Report> {
    match self {
      Self::Dense(family) => Ok(Some(family.edge_subs(graph, edge_key)?.len())),
      Self::Sparse(family) => Ok(family.edges.estimates.get(&edge_key).map(Vec::len)),
    }
  }

  fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<InDel> {
    match self {
      Self::Dense(family) => family.edge_indels(edge_key),
      Self::Sparse(family) => family.edge_indels(edge_key),
    }
  }

  pub fn root_sequence(&self, graph: &Graph) -> Result<Seq, Report> {
    match self {
      Self::Dense(family) => family.root_sequence(graph),
      Self::Sparse(family) => family.root_sequence(graph),
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
