use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::MutationCounts;
use crate::partition::marginal::discrete::input::{one_hot_profile, uniform_profile, validate_trait_names};
use crate::partition::marginal::shared::data::{DenseInputs, count_transitions_dense};
use crate::partition::marginal::shared::pass::{IndexedKind, indexed_backward, indexed_forward};
use crate::partition::marginal::shared::update::{MarginalBackward, MarginalForward, MarginalPasses};
use crate::partition::storage::dense::{
  DenseEdgeBackward, DenseEdgeEstimate, DenseEdgeForward, DenseNodeState, DenseSeqDistribution,
};
use crate::partition::storage::discrete::DiscreteStates;
use eyre::Report;
use ndarray::Array1;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::LogLh;
use treetime_utils::array::ndarray::argmax_first;

#[derive(Clone, Debug, Serialize)]
pub struct PartitionMarginalDiscrete {
  pub inputs: DenseInputs,
  pub states: DiscreteStates,
}

impl PartitionMarginalDiscrete {
  pub fn new(states: DiscreteStates, min_branch_length: f64, filter_uninformative_root: bool) -> Self {
    Self {
      inputs: DenseInputs {
        min_branch_length,
        filter_uninformative_root,
      },
      states,
    }
  }

  pub fn n_states(&self) -> usize {
    self.states.len()
  }

  pub fn get_sequence_length(&self) -> usize {
    1
  }

  pub(crate) fn attach_traits(
    &self,
    graph: &Graph,
    traits: &BTreeMap<String, String>,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
  ) -> Result<BTreeMap<GraphNodeKey, DenseNodeState>, Report> {
    let n_states = self.n_states();
    validate_trait_names(graph, traits, names)?;

    let mut node_states = BTreeMap::new();
    for leaf in graph.get_leaves() {
      let leaf_key = leaf.key();
      let leaf_name = names[&leaf_key].clone().unwrap_or_default();

      let profile = if let Some(trait_value) = traits.get(&leaf_name) {
        if let Some(index) = self.states.get_index(trait_value) {
          one_hot_profile(index, n_states)
        } else {
          uniform_profile(n_states)
        }
      } else {
        uniform_profile(n_states)
      };

      node_states.insert(
        leaf_key,
        DenseNodeState {
          seq: crate::partition::storage::dense::DenseSeqInfo::default(),
          profile: DenseSeqDistribution::new(profile, LogLh::ZERO),
        },
      );
    }

    Ok(node_states)
  }

  pub fn get_reconstructed_trait(
    &self,
    node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
    node_key: GraphNodeKey,
  ) -> Option<String> {
    let node = node_states.get(&node_key)?;
    let row = node.profile.dis.row(0);
    let argmax = argmax_first(&row)?;
    Some(self.states.get_name(argmax).to_owned())
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

impl MarginalPasses for PartitionMarginalDiscrete {
  type Node = DenseNodeState;
  type Backward = DenseEdgeBackward;
  type Forward = DenseEdgeForward;
  type Estimate = DenseEdgeEstimate;

  fn marginal_backward(
    &self,
    gtr: &GTR,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
  ) -> Result<MarginalBackward<DenseNodeState, DenseEdgeBackward>, Report> {
    indexed_backward(
      &self.inputs,
      gtr,
      None,
      1,
      IndexedKind::Discrete,
      graph,
      branch_lengths,
      node_states,
    )
  }

  fn marginal_forward(
    &self,
    gtr: &GTR,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
    backward: &BTreeMap<GraphEdgeKey, DenseEdgeBackward>,
  ) -> Result<MarginalForward<DenseNodeState, DenseEdgeForward, DenseEdgeEstimate>, Report> {
    indexed_forward(
      &self.inputs,
      gtr,
      None,
      IndexedKind::Discrete,
      graph,
      branch_lengths,
      node_states,
      backward,
    )
  }

  fn count_transitions(
    &self,
    gtr: &GTR,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
    backward: &BTreeMap<GraphEdgeKey, DenseEdgeBackward>,
    forward: &BTreeMap<GraphEdgeKey, DenseEdgeForward>,
  ) -> Result<MutationCounts, Report> {
    count_transitions_dense(&self.inputs, gtr, graph, branch_lengths, node_states, backward, forward)
  }
}
