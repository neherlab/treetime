use crate::gtr::gtr::GTR;
use crate::gtr::partition_model::PartitionModel;
use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::partition::marginal::mixture::site_log_lh::{DenseMixture, dense_mixture};
use crate::partition::marginal::shared::pass::{
  IndexedForwardMessages, IndexedPartition, indexed_backward_with_site_log_lh, indexed_forward_messages,
  indexed_forward_sequences,
};
use crate::partition::marginal::shared::update::{MarginalBackward, MarginalEdges, MarginalUpdate};
use crate::partition::storage::dense::{DenseEdgeBackward, DenseEdgeEstimate, DenseEdgeForward, DenseNodeState};
use eyre::Report;
use itertools::Itertools;
use ndarray::{Axis, Zip};
use treetime_graph::graph::Graph;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_primitives::LogLh;
use treetime_utils::{make_internal_error, make_internal_report};

pub(crate) fn dense_mixture_update(
  partition: &PartitionMarginalDense,
  model: &PartitionModel,
  graph: &Graph,
  branch_lengths: &EdgeMap<f64>,
) -> Result<(DenseMixtureUpdate, DenseMixtureState), Report> {
  let mut profiles = vec![];
  let mut categories = vec![];
  let mut site_log_lh = vec![];
  for gtr in model.category_gtrs()? {
    let (MarginalBackward { node_states, backward }, category_site_log_lh) =
      indexed_backward_with_site_log_lh(partition, &gtr, graph, branch_lengths, partition.length)?;
    site_log_lh.push(category_site_log_lh);
    let IndexedForwardMessages { node_states, forward } =
      indexed_forward_messages(partition.inputs(), &gtr, graph, branch_lengths, node_states, &backward)?;
    profiles.push(node_states);
    categories.push(DenseCategory { gtr, backward, forward });
  }
  let mixture = dense_mixture(&site_log_lh, model.gamma().weights())?;
  let node_states = mix_dense_node_states(profiles, &mixture)?;
  let (node_states, estimates) = indexed_forward_sequences(partition.indexed_kind(), graph, node_states)?;
  let update = MarginalUpdate {
    node_states,
    edges: MarginalEdges {
      backward: EdgeMap::new(),
      forward: EdgeMap::new(),
      estimates,
    },
    log_lh: mixture.log_lh,
  };
  Ok((
    update,
    DenseMixtureState {
      categories,
      mixture,
      branch_lengths: branch_lengths.clone(),
    },
  ))
}

pub(crate) fn dense_mixture_log_lh(
  partition: &PartitionMarginalDense,
  model: &PartitionModel,
  graph: &Graph,
  branch_lengths: &EdgeMap<f64>,
) -> Result<LogLh, Report> {
  let site_log_lh = model
    .category_gtrs()?
    .iter()
    .map(|gtr| {
      let (_, site_log_lh) =
        indexed_backward_with_site_log_lh(partition, gtr, graph, branch_lengths, partition.length)?;
      Ok(site_log_lh)
    })
    .try_collect::<_, Vec<_>, Report>()?;
  Ok(dense_mixture(&site_log_lh, model.gamma().weights())?.log_lh)
}

pub(crate) type DenseMixtureUpdate =
  MarginalUpdate<DenseNodeState, DenseEdgeBackward, DenseEdgeForward, DenseEdgeEstimate>;

#[derive(Clone, Debug)]
pub(crate) struct DenseMixtureState {
  pub(crate) categories: Vec<DenseCategory>,
  pub(crate) mixture: DenseMixture,
  pub(crate) branch_lengths: EdgeMap<f64>,
}

#[derive(Clone, Debug)]
pub(crate) struct DenseCategory {
  pub(crate) gtr: GTR,
  pub(crate) backward: EdgeMap<DenseEdgeBackward>,
  pub(crate) forward: EdgeMap<DenseEdgeForward>,
}

fn mix_dense_node_states(
  categories: Vec<NodeMap<DenseNodeState>>,
  mixture: &DenseMixture,
) -> Result<NodeMap<DenseNodeState>, Report> {
  let mut categories = categories.into_iter().enumerate();
  let Some((_, mut mixed)) = categories.next() else {
    return make_internal_error!("A rate mixture needs at least one category");
  };
  let weight = |k: usize| mixture.posterior.column(k).insert_axis(Axis(1));
  for node in mixed.values_mut() {
    node.profile.dis *= &weight(0);
    node.profile.log_lh = mixture.log_lh;
  }
  for (k, category) in categories {
    for (key, node) in category {
      let target = mixed
        .get_mut(key)
        .ok_or_else(|| make_internal_report!("Rate category {k} has node {key}, which the first category lacks"))?;
      Zip::from(&mut target.profile.dis)
        .and(&node.profile.dis)
        .and_broadcast(&weight(k))
        .for_each(|total, &value, &weight| *total += value * weight);
    }
  }
  Ok(mixed)
}
