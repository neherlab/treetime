use crate::alphabet::state_alphabet::StateAlphabet;
use crate::partition::marginal::shared::data::DenseInputs;
use crate::partition::marginal::shared::pass::{IndexedKind, IndexedPartition};
use crate::partition::storage::dense::DenseNodeState;
use crate::{make_internal_error, make_internal_report};
use eyre::Report;
use ndarray::{Array1, Array2};
use std::collections::BTreeMap;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::array::ndarray::argmax_first;

#[derive(Clone, Debug)]
pub struct PartitionMarginalDiscrete {
  pub(crate) inputs: DenseInputs,
  pub(crate) states: StateAlphabet,
  pub(crate) obs_leaves: BTreeMap<GraphNodeKey, Array2<f64>>,
}

impl PartitionMarginalDiscrete {
  pub(crate) fn new(
    states: StateAlphabet,
    graph: &Graph,
    observations: &BTreeMap<GraphNodeKey, StateObservation>,
    min_branch_length: f64,
    filter_uninformative_root: bool,
  ) -> Result<Self, Report> {
    let obs_leaves = graph
      .get_leaves()
      .map(|leaf| {
        let leaf_key = leaf.key();
        let observation = observations
          .get(&leaf_key)
          .ok_or_else(|| make_internal_report!("Leaf {leaf_key} has no state observation"))?;
        Ok((leaf_key, leaf_profile(*observation, states.len())?))
      })
      .collect::<Result<_, Report>>()?;
    Ok(Self {
      inputs: DenseInputs {
        min_branch_length,
        filter_uninformative_root,
      },
      states,
      obs_leaves,
    })
  }

  pub fn n_states(&self) -> usize {
    self.states.len()
  }

  pub fn get_reconstructed_trait(
    &self,
    node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
    node_key: GraphNodeKey,
  ) -> Option<String> {
    let node = node_states.get(&node_key)?;
    let row = node.profile.dis.row(0);
    let argmax = argmax_first(&row)?;
    self.states.name(argmax).map(ToOwned::to_owned)
  }

  pub fn get_confidence(
    &self,
    node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
    node_key: GraphNodeKey,
  ) -> Option<Array1<f64>> {
    let node = node_states.get(&node_key)?;
    Some(node.profile.dis.row(0).to_owned())
  }
}

impl IndexedPartition for PartitionMarginalDiscrete {
  fn inputs(&self) -> &DenseInputs {
    &self.inputs
  }

  fn indexed_kind(&self) -> IndexedKind<'_> {
    IndexedKind::Discrete {
      leaves: &self.obs_leaves,
    }
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StateObservation {
  Known(usize),
  Unknown,
}

fn leaf_profile(observation: StateObservation, n_states: usize) -> Result<Array2<f64>, Report> {
  match observation {
    StateObservation::Known(index) if index >= n_states => {
      make_internal_error!("Observed state index {index} is outside the {n_states} states of the partition")
    },
    StateObservation::Known(index) => Ok(Array2::from_shape_fn((1, n_states), |(_, state)| {
      if state == index { 1.0 } else { 0.0 }
    })),
    StateObservation::Unknown => Ok(Array2::ones((1, n_states))),
  }
}
