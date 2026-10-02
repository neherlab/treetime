use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::MutationCounts;
use crate::partition::marginal::discrete::input::{missing_trait_profile, one_hot_profile, validate_trait_names};
use crate::partition::marginal::shared::data::{DenseInputs, count_transitions_dense};
use crate::partition::marginal::shared::pass::{IndexedKind, indexed_backward, indexed_forward};
use crate::partition::marginal::shared::update::{MarginalBackward, MarginalForward, MarginalPasses};
use crate::partition::storage::dense::{DenseEdgeBackward, DenseEdgeEstimate, DenseEdgeForward, DenseNodeState};
use crate::partition::storage::discrete::DiscreteStates;
use crate::progress::LogSink;
use eyre::Report;
use ndarray::{Array1, Array2};
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::array::ndarray::argmax_first;

#[derive(Clone, Debug, Serialize)]
pub struct PartitionMarginalDiscrete {
  pub(crate) inputs: DenseInputs,
  pub(crate) states: DiscreteStates,
  pub(crate) obs_leaves: BTreeMap<GraphNodeKey, Array2<f64>>,
}

impl PartitionMarginalDiscrete {
  pub(crate) fn new(
    states: DiscreteStates,
    graph: &Graph,
    traits: &BTreeMap<String, String>,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    min_branch_length: f64,
    filter_uninformative_root: bool,
    log: &dyn LogSink,
  ) -> Result<Self, Report> {
    validate_trait_names(graph, traits, names, log)?;
    let n_states = states.len();
    let obs_leaves = graph
      .get_leaves()
      .map(|leaf| {
        let leaf_key = leaf.key();
        let leaf_name = names[&leaf_key].clone().unwrap_or_default();
        let profile = traits
          .get(&leaf_name)
          .and_then(|trait_value| states.get_index(trait_value))
          .map_or_else(
            || missing_trait_profile(n_states),
            |index| one_hot_profile(index, n_states),
          );
        (leaf_key, profile)
      })
      .collect();
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

  fn indexed_kind(&self) -> IndexedKind<'_> {
    IndexedKind::Discrete {
      leaves: &self.obs_leaves,
    }
  }
}

impl MarginalPasses for PartitionMarginalDiscrete {
  type Node = DenseNodeState;
  type BackwardInput = ();
  type Backward = DenseEdgeBackward;
  type Forward = DenseEdgeForward;
  type Estimate = DenseEdgeEstimate;

  fn backward_input(_node_states: &BTreeMap<GraphNodeKey, DenseNodeState>) -> &Self::BackwardInput {
    &()
  }

  fn fresh_backward_input(_node_states: &BTreeMap<GraphNodeKey, DenseNodeState>) -> Self::BackwardInput {}

  fn marginal_backward(
    &self,
    gtr: &GTR,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    (): &(),
  ) -> Result<MarginalBackward<DenseNodeState, DenseEdgeBackward>, Report> {
    indexed_backward(&self.inputs, gtr, self.indexed_kind(), graph, branch_lengths)
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
      self.indexed_kind(),
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
