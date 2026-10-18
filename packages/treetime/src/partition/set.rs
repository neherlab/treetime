use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::partition::create::{Representation, build_marginal_partition};
use crate::partition::marginal::reconstruction::{
  DenseReconstruction, MarginalReconstruction, ReconstructionRef, SparseReconstruction,
};
use crate::partition::result::{PartitionSource, SplitSource, split_sources};
use crate::progress::LogSink;
use crate::seq::alignment::SequenceSource;
use eyre::Report;
use itertools::{Either, Itertools, izip};
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::reroot::RerootResult;
use treetime_primitives::LogLh;
use treetime_utils::{make_internal_error, make_internal_report};

pub(crate) fn build_partition_set(
  representation: Representation,
  model: GtrModelName,
  graph: &Graph,
  sources: Vec<SequenceSource>,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  log: &dyn LogSink,
) -> Result<(Vec<PartitionSource>, PartitionSet), Report> {
  let (sources, reconstructions): (Vec<_>, Vec<_>) = split_sources(sources)
    .into_iter()
    .map(
      |SplitSource {
         source,
         alphabet,
         nodes,
       }| {
        let reconstruction =
          build_marginal_partition(representation, model, graph, alphabet, nodes, branch_lengths, log)?;
        Ok((source, reconstruction))
      },
    )
    .collect::<Result<Vec<_>, Report>>()?
    .into_iter()
    .unzip();
  let set = PartitionSet::from_reconstructions(reconstructions)?;
  let total_length = set.total_length();
  Ok((sources, set.with_total_length(total_length)))
}

#[derive(Clone, Debug)]
pub enum PartitionSet {
  Dense(Vec<DenseReconstruction>),
  Sparse(Vec<SparseReconstruction>),
}

impl PartitionSet {
  pub(crate) fn from_reconstructions(reconstructions: Vec<MarginalReconstruction>) -> Result<Self, Report> {
    let (dense, sparse): (Vec<_>, Vec<_>) =
      reconstructions
        .into_iter()
        .partition_map(|reconstruction| match reconstruction {
          MarginalReconstruction::Dense(dense) => Either::Left(dense),
          MarginalReconstruction::Sparse(sparse) => Either::Right(sparse),
        });
    match (dense.is_empty(), sparse.is_empty()) {
      (_, true) => Ok(Self::Dense(dense)),
      (true, false) => Ok(Self::Sparse(sparse)),
      (false, false) => Err(make_internal_report!(
        "A partition set mixes dense and sparse partitions; every partition of a run uses one representation"
      )),
    }
  }

  pub(crate) fn len(&self) -> usize {
    match self {
      Self::Dense(partitions) => partitions.len(),
      Self::Sparse(partitions) => partitions.len(),
    }
  }

  pub(crate) fn iter(&self) -> impl Iterator<Item = ReconstructionRef<'_>> {
    match self {
      Self::Dense(partitions) => Either::Left(partitions.iter().map(ReconstructionRef::Dense)),
      Self::Sparse(partitions) => Either::Right(partitions.iter().map(ReconstructionRef::Sparse)),
    }
  }

  pub(crate) fn total_length(&self) -> usize {
    self.iter().map(ReconstructionRef::length).sum()
  }

  pub(crate) fn graph_log_lh(&self, graph: &Graph) -> Result<LogLh, Report> {
    self.iter().map(|partition| partition.graph_log_lh(graph)).sum()
  }

  pub(crate) fn edge_sub_count(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<Option<usize>, Report> {
    self
      .iter()
      .map(|partition| partition.edge_sub_count(graph, edge_key))
      .sum()
  }

  pub(crate) fn marginal_update(
    self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  ) -> Result<(Self, LogLh), Report> {
    Ok(match self {
      Self::Dense(partitions) => {
        let (partitions, log_lh) =
          update_each(partitions, |partition| partition.marginal_update(graph, branch_lengths))?;
        (Self::Dense(partitions), log_lh)
      },
      Self::Sparse(partitions) => {
        let (partitions, log_lh) =
          update_each(partitions, |partition| partition.marginal_update(graph, branch_lengths))?;
        (Self::Sparse(partitions), log_lh)
      },
    })
  }

  pub(crate) fn apply_reroot(self, changes: &RerootResult) -> Result<Self, Report> {
    Ok(match self {
      Self::Dense(partitions) => Self::Dense(partitions.into_iter().map(DenseReconstruction::reseeded).collect()),
      Self::Sparse(partitions) => Self::Sparse(
        partitions
          .into_iter()
          .map(|partition| partition.apply_reroot(changes))
          .try_collect()?,
      ),
    })
  }

  #[must_use]
  pub(crate) fn reconcile_topology(self, graph: &Graph) -> Self {
    match self {
      Self::Dense(partitions) => Self::Dense(partitions.into_iter().map(DenseReconstruction::reseeded).collect()),
      Self::Sparse(partitions) => Self::Sparse(
        partitions
          .into_iter()
          .map(|partition| partition.reconcile_topology(graph))
          .collect(),
      ),
    }
  }

  #[must_use]
  pub(crate) fn reconcile_node_states(self, graph: &Graph) -> Self {
    match self {
      Self::Dense(partitions) => Self::Dense(partitions.into_iter().map(DenseReconstruction::reseeded).collect()),
      Self::Sparse(partitions) => Self::Sparse(
        partitions
          .into_iter()
          .map(|partition| partition.reconcile_node_states(graph))
          .collect(),
      ),
    }
  }

  pub(crate) fn estimate_rates(
    self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    log: &dyn LogSink,
  ) -> Result<Self, Report> {
    Ok(match self {
      Self::Dense(partitions) => Self::Dense(
        partitions
          .into_iter()
          .map(|partition| partition.with_estimated_rate(graph, branch_lengths, log))
          .try_collect()?,
      ),
      Self::Sparse(partitions) => Self::Sparse(
        partitions
          .into_iter()
          .map(|partition| partition.with_estimated_rate(graph, branch_lengths, log))
          .try_collect()?,
      ),
    })
  }

  pub(crate) fn rates(&self) -> Vec<f64> {
    self.iter().map(|partition| partition.gtr().mu()).collect()
  }

  pub(crate) fn set_rates(&mut self, rates: &[f64]) -> Result<(), Report> {
    let mut gtrs = self.lengths_and_gtrs_mut();
    if gtrs.len() != rates.len() {
      return make_internal_error!(
        "Cannot set {} relative rates on a partition set with {} partitions",
        rates.len(),
        gtrs.len()
      );
    }
    izip!(&mut gtrs, rates).try_for_each(|((_, gtr), &rate)| gtr.set_mu(rate))
  }

  pub(crate) fn sparse_mut(&mut self) -> &mut [SparseReconstruction] {
    match self {
      Self::Dense(_) => &mut [],
      Self::Sparse(partitions) => partitions,
    }
  }

  pub(crate) fn lengths_and_gtrs_mut(&mut self) -> Vec<(usize, &mut GTR)> {
    match self {
      Self::Dense(partitions) => partitions
        .iter_mut()
        .map(|partition| (partition.partition.length, &mut partition.gtr))
        .collect(),
      Self::Sparse(partitions) => partitions
        .iter_mut()
        .map(|partition| (partition.partition.length, &mut partition.gtr))
        .collect(),
    }
  }

  #[must_use]
  fn with_total_length(self, total_length: usize) -> Self {
    match self {
      Self::Dense(partitions) => Self::Dense(
        partitions
          .into_iter()
          .map(|partition| partition.with_total_length(total_length))
          .collect(),
      ),
      Self::Sparse(partitions) => Self::Sparse(
        partitions
          .into_iter()
          .map(|partition| partition.with_total_length(total_length))
          .collect(),
      ),
    }
  }
}

fn update_each<R>(
  partitions: Vec<R>,
  update: impl Fn(R) -> Result<(R, LogLh), Report>,
) -> Result<(Vec<R>, LogLh), Report> {
  let updated: Vec<(R, LogLh)> = partitions.into_iter().map(update).try_collect()?;
  let log_lh = updated.iter().map(|(_, log_lh)| *log_lh).sum();
  Ok((updated.into_iter().map(|(partition, _)| partition).collect(), log_lh))
}
