use crate::alphabet::state_alphabet::StateAlphabet;
use crate::partition::marginal::discrete::state_posterior::StatePosterior;
use crate::partition::marginal::shared::data::DenseInputs;
use crate::partition::marginal::shared::pass::{IndexedKind, IndexedPartition};
use crate::partition::storage::dense::DenseNodeState;
use crate::{make_internal_error, make_internal_report};
use eyre::Report;
use ndarray::Array2;
use std::collections::BTreeMap;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

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

  pub(crate) fn state_posteriors(
    &self,
    graph: &Graph,
    node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
  ) -> Result<BTreeMap<GraphNodeKey, StatePosterior>, Report> {
    graph
      .get_nodes()
      .map(|node| {
        let key = node.key();
        let profile = node_states[&key].profile.dis.row(0).to_owned();
        if profile.len() != self.states.len() {
          return make_internal_error!(
            "Posterior profile of node {key} has {} states, the partition has {}",
            profile.len(),
            self.states.len()
          );
        }
        Ok((key, StatePosterior::new(profile)?))
      })
      .collect()
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
