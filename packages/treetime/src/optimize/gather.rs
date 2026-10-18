use crate::partition::marginal::reconstruction::{MarginalReconstruction, ReconstructionRef};
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::partition::set::PartitionSet;
use eyre::Report;
use itertools::Itertools;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use smallvec::SmallVec;
use std::collections::BTreeMap;
use std::iter;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;

pub(crate) trait EdgeContributions: Sync {
  fn with_edge<R>(
    &self,
    edge_key: GraphEdgeKey,
    f: impl FnOnce(&[OptimizationContribution]) -> Result<R, Report>,
  ) -> Result<R, Report>;
}

pub(crate) trait SequencePartitions: Sync {
  fn partitions(&self) -> impl Iterator<Item = ReconstructionRef<'_>>;
}

impl SequencePartitions for PartitionSet {
  fn partitions(&self) -> impl Iterator<Item = ReconstructionRef<'_>> {
    self.iter()
  }
}

impl SequencePartitions for MarginalReconstruction {
  fn partitions(&self) -> impl Iterator<Item = ReconstructionRef<'_>> {
    iter::once(self.view())
  }
}

impl<P: SequencePartitions> EdgeContributions for P {
  fn with_edge<R>(
    &self,
    edge_key: GraphEdgeKey,
    f: impl FnOnce(&[OptimizationContribution]) -> Result<R, Report>,
  ) -> Result<R, Report> {
    let contributions: SmallVec<[OptimizationContribution; 1]> = self
      .partitions()
      .map(|partition| partition.create_edge_contribution(edge_key))
      .try_collect()?;
    f(&contributions)
  }
}

pub(crate) fn total_length(partitions: &impl SequencePartitions) -> usize {
  partitions.partitions().map(ReconstructionRef::length).sum()
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
  partitions: &impl SequencePartitions,
) -> Result<BTreeMap<GraphEdgeKey, Vec<OptimizationContribution>>, Report> {
  graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(
      |edge_ref| -> Result<(GraphEdgeKey, Vec<OptimizationContribution>), Report> {
        let edge_key = edge_ref.key();
        let contributions = partitions
          .partitions()
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
  pub(crate) fn gather(graph: &Graph, partitions: &impl SequencePartitions) -> Result<Self, Report> {
    Ok(Self {
      total_length: total_length(partitions),
      sub_counts: gather_edge_sums(graph, partitions, |partition, edge_key| {
        Ok(partition.edge_subs(graph, edge_key)?.len())
      })?,
      effective_lengths: gather_edge_sums(graph, partitions, |partition, edge_key| {
        partition.edge_effective_length(graph, edge_key)
      })?,
    })
  }
}

pub(crate) fn gather_edge_indel_counts(
  graph: &Graph,
  partitions: &impl SequencePartitions,
) -> BTreeMap<GraphEdgeKey, usize> {
  graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|edge_ref| {
      let edge_key = edge_ref.key();
      let count = partitions
        .partitions()
        .map(|partition| partition.edge_indel_count(edge_key))
        .sum();
      (edge_key, count)
    })
    .collect()
}

pub(crate) fn gather_partition_edge_indel_counts(
  graph: &Graph,
  partitions: &impl SequencePartitions,
) -> Vec<BTreeMap<GraphEdgeKey, usize>> {
  partitions
    .partitions()
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

fn gather_edge_sums<P: SequencePartitions>(
  graph: &Graph,
  partitions: &P,
  value: impl Fn(ReconstructionRef<'_>, GraphEdgeKey) -> Result<usize, Report> + Sync,
) -> Result<BTreeMap<GraphEdgeKey, usize>, Report> {
  graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|edge_ref| -> Result<(GraphEdgeKey, usize), Report> {
      let edge_key = edge_ref.key();
      let sum = partitions
        .partitions()
        .map(|partition| value(partition, edge_key))
        .sum::<Result<usize, Report>>()?;
      Ok((edge_key, sum))
    })
    .collect()
}
