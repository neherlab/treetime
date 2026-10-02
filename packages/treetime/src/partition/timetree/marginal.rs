use crate::ancestral::reconstruction::Reconstruction;
use crate::ancestral::sample::SampleMode;
use crate::ancestral::tip_states::TipStates;
use crate::partition::marginal::shared::update::MarginalPasses;
use crate::partition::timetree::partition::PartitionTimetree;
use crate::seq::indel::InDel;
use crate::seq::mutation::{Mutation, MutationTrack, Sub, combine_edge_mutations};
use eyre::Report;
use rand::RngCore;
use rayon::prelude::*;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{LogLh, Seq};

pub(crate) fn graph_log_lh(graph: &Graph, partitions: &[PartitionTimetree]) -> Result<LogLh, Report> {
  let root_key = graph.get_exactly_one_root()?.key();
  let log_lh = partitions
    .par_iter()
    .map(|partition| partition.get_log_lh(root_key))
    .collect::<Vec<_>>()
    .into_iter()
    .sum();
  Ok(log_lh)
}

pub(crate) fn marginal_update_timetree(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  partitions: Vec<PartitionTimetree>,
) -> Result<(Vec<PartitionTimetree>, LogLh), Report> {
  partitions
    .into_iter()
    .try_fold((Vec::new(), LogLh::ZERO), |(mut updated, total), partition| {
      let (partition, log_lh) = partition.marginal_update(graph, branch_lengths)?;
      updated.push(partition);
      Ok((updated, total + log_lh))
    })
}

impl PartitionTimetree {
  pub(crate) fn get_sequence_length(&self) -> usize {
    match self {
      Self::Dense(family) => family.partition.length,
      Self::Sparse(family) => family.partition.length,
    }
  }

  fn get_log_lh(&self, node_key: GraphNodeKey) -> LogLh {
    match self {
      Self::Dense(family) => family.partition.get_log_lh(&family.node_states, node_key),
      Self::Sparse(family) => family.partition.get_log_lh(&family.node_states, node_key),
    }
  }

  fn marginal_update(
    self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  ) -> Result<(Self, LogLh), Report> {
    match self {
      Self::Dense(family) => {
        let (family, log_lh) = family.marginal_update(graph, branch_lengths)?;
        Ok((Self::Dense(family), log_lh))
      },
      Self::Sparse(family) => {
        let (family, log_lh) = family.marginal_update(graph, branch_lengths)?;
        Ok((Self::Sparse(family), log_lh))
      },
    }
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

  pub(crate) fn reconstruct_sequences(
    &self,
    graph: &Graph,
    tips: TipStates,
    sample_mode: SampleMode,
    rng: &mut dyn RngCore,
  ) -> Result<Reconstruction, Report> {
    match self {
      Self::Dense(family) => family.reconstruct_sequences(graph, tips, sample_mode, rng),
      Self::Sparse(family) => family.reconstruct_sequences(graph, tips, sample_mode, rng),
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
