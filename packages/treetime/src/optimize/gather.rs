use crate::partition::marginal::reconstruction::ReconstructionRef;
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::partition::set::PartitionSet;
use eyre::Report;
use itertools::Itertools;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use smallvec::SmallVec;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;

pub(crate) trait EdgeContributions: Sync {
  fn with_edge<R>(
    &self,
    edge_key: GraphEdgeKey,
    f: impl FnOnce(&[OptimizationContribution]) -> Result<R, Report>,
  ) -> Result<R, Report>;
}

impl EdgeContributions for PartitionSet {
  fn with_edge<R>(
    &self,
    edge_key: GraphEdgeKey,
    f: impl FnOnce(&[OptimizationContribution]) -> Result<R, Report>,
  ) -> Result<R, Report> {
    let contributions: SmallVec<[OptimizationContribution; 1]> = self
      .iter()
      .map(|partition| partition.create_edge_contribution(edge_key))
      .try_collect()?;
    f(&contributions)
  }
}

#[cfg(test)]
impl EdgeContributions for BTreeMap<GraphEdgeKey, Vec<OptimizationContribution>> {
  fn with_edge<R>(
    &self,
    edge_key: GraphEdgeKey,
    f: impl FnOnce(&[OptimizationContribution]) -> Result<R, Report>,
  ) -> Result<R, Report> {
    f(&self[&edge_key])
  }
}

#[cfg(test)]
pub(crate) fn gather_edge_contributions(
  graph: &Graph,
  partitions: &PartitionSet,
) -> Result<BTreeMap<GraphEdgeKey, Vec<OptimizationContribution>>, Report> {
  graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(
      |edge_ref| -> Result<(GraphEdgeKey, Vec<OptimizationContribution>), Report> {
        let edge_key = edge_ref.key();
        let contributions = partitions
          .iter()
          .map(|partition| partition.create_edge_contribution(edge_key))
          .try_collect()?;
        Ok((edge_key, contributions))
      },
    )
    .collect()
}

pub(crate) struct InitialGuessEvidence {
  pub total_length: usize,
  pub sub_counts: BTreeMap<GraphEdgeKey, usize>,
  pub effective_lengths: BTreeMap<GraphEdgeKey, usize>,
}

impl InitialGuessEvidence {
  pub(crate) fn gather(graph: &Graph, partitions: &PartitionSet) -> Result<Self, Report> {
    Ok(Self {
      total_length: partitions.total_length(),
      sub_counts: gather_edge_sums(graph, partitions, |partition, edge_key| {
        Ok(partition.edge_subs(graph, edge_key)?.len())
      })?,
      effective_lengths: gather_edge_sums(graph, partitions, |partition, edge_key| {
        partition.edge_effective_length(graph, edge_key)
      })?,
    })
  }
}

pub(crate) fn gather_edge_indel_counts(graph: &Graph, partitions: &PartitionSet) -> BTreeMap<GraphEdgeKey, usize> {
  graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|edge_ref| {
      let edge_key = edge_ref.key();
      let count = partitions
        .iter()
        .map(|partition| partition.edge_indel_count(edge_key))
        .sum();
      (edge_key, count)
    })
    .collect()
}

pub(crate) fn gather_partition_edge_indel_counts(
  graph: &Graph,
  partitions: &PartitionSet,
) -> Vec<BTreeMap<GraphEdgeKey, usize>> {
  partitions
    .iter()
    .map(|partition| {
      graph
        .get_edges()
        .map(|edge_ref| {
          let edge_key = edge_ref.key();
          (edge_key, partition.edge_indel_count(edge_key))
        })
        .collect()
    })
    .collect()
}

fn gather_edge_sums(
  graph: &Graph,
  partitions: &PartitionSet,
  value: impl Fn(ReconstructionRef<'_>, GraphEdgeKey) -> Result<usize, Report> + Sync,
) -> Result<BTreeMap<GraphEdgeKey, usize>, Report> {
  graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|edge_ref| -> Result<(GraphEdgeKey, usize), Report> {
      let edge_key = edge_ref.key();
      let sum = partitions
        .iter()
        .map(|partition| value(partition, edge_key))
        .sum::<Result<usize, Report>>()?;
      Ok((edge_key, sum))
    })
    .collect()
}
