use crate::gtr::gamma::{GammaRates, INITIAL_GAMMA_SHAPE};
use crate::gtr::gamma_shape::{ShapeEstimate, optimize_gamma_shape};
use crate::gtr::gtr::GTR;
use crate::gtr::partition_model::PartitionModel;
use crate::partition::assembly::{AlignmentGroup, alignment_groups};
use crate::partition::create::{PartitionInput, Representation, build_marginal_partitions};
use crate::partition::inferred_rate::InferredRate;
use crate::partition::marginal::reconstruction::{
  DenseReconstruction, MarginalReconstruction, ReconstructionRef, SparseReconstruction,
};
use crate::partition::marginal::sparse::alignment_indels::merge_sparse_edges;
use crate::partition::marginal::sparse::reroot::reroot_sparse;
use crate::partition::resolve::PartitionLayout;
use crate::partition::split::{ModelSettings, PartitionData, SplitPartition};
use crate::progress::LogSink;
use crate::progress_info;
use eyre::Report;
use itertools::{Either, Itertools, izip};
use log::debug;
use std::borrow::Cow;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::EdgeMap;
use treetime_graph::reroot::{EdgeMergeInfo, RerootResult};
use treetime_graph::scope::Inheritance;
use treetime_primitives::LogLh;
use treetime_utils::{make_internal_error, make_internal_report};

pub(crate) fn build_partition_set(
  representation: Representation,
  graph: &Graph,
  partitions: Vec<SplitPartition>,
  models: ModelSettings,
  branch_lengths: &EdgeMap<f64>,
  log: &dyn LogSink,
) -> Result<PartitionSet, Report> {
  let inheritances = vec![Inheritance::ALL; partitions.len()];
  build_inherited_partition_set(
    representation,
    graph,
    partitions,
    inheritances,
    models,
    branch_lengths,
    log,
  )
}

pub(crate) fn build_inherited_partition_set(
  representation: Representation,
  graph: &Graph,
  partitions: Vec<SplitPartition>,
  inheritances: Vec<Inheritance>,
  models: ModelSettings,
  branch_lengths: &EdgeMap<f64>,
  log: &dyn LogSink,
) -> Result<PartitionSet, Report> {
  let rate = InferredRate::for_partition_count(partitions.len());
  let (heads, inputs): (Vec<_>, Vec<_>) = izip!(partitions, inheritances)
    .map(|(SplitPartition { data, model }, inheritance)| {
      let PartitionData {
        layout,
        alphabet,
        nodes,
      } = data;
      ((layout, model), (inheritance, alphabet, nodes))
    })
    .unzip();
  let inputs = izip!(&heads, inputs)
    .map(|((_, model), (inheritance, alphabet, nodes))| PartitionInput {
      model,
      rate,
      inheritance,
      alphabet,
      nodes,
    })
    .collect_vec();
  let reconstructions = build_marginal_partitions(representation, inputs, models, graph, branch_lengths, log)?;
  let reconstructions = izip!(heads, reconstructions)
    .map(|((layout, _), reconstruction)| reconstruction.with_layout(layout))
    .collect_vec();
  let set = PartitionSet::from_reconstructions(reconstructions)?;
  let total_length = set.total_length();
  Ok(set.with_total_length(total_length))
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

  pub(crate) fn into_reconstructions(self) -> Vec<MarginalReconstruction> {
    match self {
      Self::Dense(partitions) => partitions.into_iter().map(MarginalReconstruction::Dense).collect(),
      Self::Sparse(partitions) => partitions.into_iter().map(MarginalReconstruction::Sparse).collect(),
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

  pub(crate) fn layouts(&self) -> Vec<Cow<'_, PartitionLayout>> {
    layouts_or_whole(self.iter().map(|partition| (partition.layout(), partition.length())))
  }

  pub(crate) fn alignment_groups(&self) -> Result<Vec<AlignmentGroup>, Report> {
    alignment_groups(self.iter().map(|partition| (partition.layout(), partition.length())))
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

  pub(crate) fn marginal_update(self, graph: &Graph, branch_lengths: &EdgeMap<f64>) -> Result<(Self, LogLh), Report> {
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
      Self::Sparse(partitions) => Self::Sparse(reroot_sparse(partitions, changes)?),
    })
  }

  pub(crate) fn reconcile_topology(self, graph: &Graph, merges: &[EdgeMergeInfo]) -> Result<Self, Report> {
    Ok(match self {
      Self::Dense(partitions) => Self::Dense(partitions.into_iter().map(DenseReconstruction::reseeded).collect()),
      Self::Sparse(mut partitions) => {
        merge_sparse_edges(&mut partitions, merges)?;
        Self::Sparse(
          partitions
            .into_iter()
            .map(|partition| partition.reconcile_topology(graph, merges))
            .try_collect()?,
        )
      },
    })
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
    branch_lengths: &EdgeMap<f64>,
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

  pub(crate) fn estimate_gamma_shapes(
    self,
    graph: &Graph,
    branch_lengths: &EdgeMap<f64>,
    log: &dyn LogSink,
  ) -> Result<Self, Report> {
    let groups = self.gamma_shape_groups();
    if groups.is_empty() {
      return Ok(self);
    }
    let partitions = self.iter().collect_vec();
    let mut shapes = vec![None; partitions.len()];
    for group in &groups {
      let start = partitions[group[0]]
        .model()
        .gamma()
        .shape()
        .unwrap_or(INITIAL_GAMMA_SHAPE);
      let ShapeEstimate { shape, evaluations } = optimize_gamma_shape(start, |shape| {
        group
          .iter()
          .map(|&index| {
            let partition = partitions[index];
            let gamma = GammaRates::new(partition.model().categories(), shape)?;
            partition.mixture_log_lh(&partition.model().clone().with_gamma(gamma), graph, branch_lengths)
          })
          .sum()
      })?;
      let labels = group
        .iter()
        .map(|&index| partition_label(partitions[index], index))
        .join(", ");
      progress_info!(log, "Gamma shape of {labels}: {shape:.4}");
      debug!("Gamma shape of {labels}: {evaluations} likelihood evaluations, starting from {start:.4}");
      for &index in group {
        shapes[index] = Some(shape);
      }
    }
    let with_shape = |shape: Option<f64>, model: &PartitionModel| -> Result<Option<PartitionModel>, Report> {
      shape
        .map(|shape| Ok(model.clone().with_gamma(GammaRates::new(model.categories(), shape)?)))
        .transpose()
    };
    Ok(match self {
      Self::Dense(partitions) => Self::Dense(
        izip!(partitions, shapes)
          .map(|(partition, shape)| {
            Ok(match with_shape(shape, &partition.model)? {
              Some(model) => partition.with_model(model),
              None => partition,
            })
          })
          .collect::<Result<_, Report>>()?,
      ),
      Self::Sparse(partitions) => Self::Sparse(
        izip!(partitions, shapes)
          .map(|(partition, shape)| {
            Ok(match with_shape(shape, &partition.model)? {
              Some(model) => partition.with_model(model),
              None => partition,
            })
          })
          .collect::<Result<_, Report>>()?,
      ),
    })
  }

  fn gamma_shape_groups(&self) -> Vec<Vec<usize>> {
    let (shared, own): (Vec<_>, Vec<_>) = self
      .iter()
      .enumerate()
      .filter(|(_, partition)| partition.model().categories() > 1)
      .partition(|(_, partition)| partition.model().shares_gamma_shape());
    let own = own.into_iter().map(|(index, _)| vec![index]);
    let shared = (!shared.is_empty()).then(|| shared.into_iter().map(|(index, _)| index).collect_vec());
    shared.into_iter().chain(own).sorted_by_key(|group| group[0]).collect()
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
        .map(|partition| (partition.partition.length, partition.model.gtr_mut()))
        .collect(),
      Self::Sparse(partitions) => partitions
        .iter_mut()
        .map(|partition| (partition.partition.length, partition.model.gtr_mut()))
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

pub(crate) fn layouts_or_whole<'a>(
  partitions: impl IntoIterator<Item = (Option<&'a PartitionLayout>, usize)>,
) -> Vec<Cow<'a, PartitionLayout>> {
  partitions
    .into_iter()
    .zip(1..)
    .map(|((layout, length), number)| {
      layout.map_or_else(|| Cow::Owned(PartitionLayout::whole(number, length)), Cow::Borrowed)
    })
    .collect()
}

fn update_each<R>(
  partitions: Vec<R>,
  update: impl Fn(R) -> Result<(R, LogLh), Report>,
) -> Result<(Vec<R>, LogLh), Report> {
  let updated: Vec<(R, LogLh)> = partitions.into_iter().map(update).try_collect()?;
  let log_lh = updated.iter().map(|(_, log_lh)| *log_lh).sum();
  Ok((updated.into_iter().map(|(partition, _)| partition).collect(), log_lh))
}

fn partition_label(partition: ReconstructionRef<'_>, index: usize) -> String {
  partition
    .layout()
    .map_or_else(|| format!("partition {}", index + 1), PartitionLayout::label)
}
