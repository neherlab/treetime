use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::MutationCounts;
use crate::partition::marginal::discrete::input::{one_hot_profile, uniform_profile, validate_trait_names};
use crate::partition::marginal::shared::data::{DenseInputs, count_transitions_dense};
use crate::partition::marginal::shared::pass::{IndexedKind, indexed_backward, indexed_forward};
use crate::partition::storage::dense::{
  DenseEdgeBackward, DenseEdgeEstimate, DenseEdgeForward, DenseNodeState, DenseSeqDistribution,
};
use crate::partition::storage::discrete::DiscreteStates;
use crate::partition::traits::HasGtr;
use eyre::Report;
use ndarray::Array1;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::LogLh;
use treetime_utils::array::ndarray::argmax_first;

/// The discrete marginal representation as durable, borrowed inputs: the model and the state alphabet.
/// The stage-filled node states, backward/forward messages, and edge estimates are owned separately by
/// the values the passes return.
#[derive(Clone, Debug, Serialize)]
pub struct PartitionMarginalDiscrete {
  pub inputs: DenseInputs,
  pub states: DiscreteStates,
}

impl PartitionMarginalDiscrete {
  pub fn new(gtr: GTR, states: DiscreteStates, min_branch_length: f64, filter_uninformative_root: bool) -> Self {
    Self {
      inputs: DenseInputs {
        gtr,
        min_branch_length,
        filter_uninformative_root,
      },
      states,
    }
  }

  pub fn n_states(&self) -> usize {
    self.states.len()
  }

  pub fn gtr(&self) -> &GTR {
    &self.inputs.gtr
  }

  pub fn gtr_mut(&mut self) -> &mut GTR {
    &mut self.inputs.gtr
  }

  pub fn get_sequence_length(&self) -> usize {
    1
  }

  pub fn weighted_rate(&self) -> f64 {
    self.inputs.gtr.mu
  }

  pub fn normalize_rate(&mut self, scale: f64) {
    self.inputs.gtr.mu /= scale;
  }

  /// Build the initial discrete node states by attaching each leaf's trait as a one-hot (or uniform)
  /// profile. Returns the leaf-seeded node-state map.
  pub fn attach_traits(
    &self,
    graph: &Graph,
    traits: &BTreeMap<String, String>,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
  ) -> Result<BTreeMap<GraphNodeKey, DenseNodeState>, Report> {
    let n_states = self.n_states();
    validate_trait_names(graph, traits, names)?;

    let mut node_states = BTreeMap::new();
    for leaf in graph.get_leaves() {
      let leaf_key = leaf.read_arc().key();
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

  pub fn marginal_backward(
    &self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
  ) -> Result<
    (
      BTreeMap<GraphNodeKey, DenseNodeState>,
      BTreeMap<GraphEdgeKey, DenseEdgeBackward>,
    ),
    Report,
  > {
    indexed_backward(
      &self.inputs,
      // discrete carries no residue alphabet; the indexed driver only uses the alphabet on the dense
      // leaf-profile branch, which discrete never takes.
      None,
      1,
      IndexedKind::Discrete,
      graph,
      branch_lengths,
      node_states,
    )
  }

  pub fn marginal_forward(
    &self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
    backward: &BTreeMap<GraphEdgeKey, DenseEdgeBackward>,
  ) -> Result<
    (
      BTreeMap<GraphNodeKey, DenseNodeState>,
      BTreeMap<GraphEdgeKey, DenseEdgeForward>,
      BTreeMap<GraphEdgeKey, DenseEdgeEstimate>,
    ),
    Report,
  > {
    indexed_forward(
      &self.inputs,
      None,
      IndexedKind::Discrete,
      graph,
      branch_lengths,
      node_states,
      backward,
    )
  }

  /// Run a full marginal update (backward, then forward) over the discrete representation, returning the
  /// updated node states, the per-edge backward messages, forward messages, and estimates as distinct
  /// owned values, plus the substitution log likelihood (root node likelihood after the backward pass).
  #[allow(clippy::type_complexity)]
  pub fn marginal_update(
    &self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: BTreeMap<GraphNodeKey, DenseNodeState>,
  ) -> Result<
    (
      BTreeMap<GraphNodeKey, DenseNodeState>,
      BTreeMap<GraphEdgeKey, DenseEdgeBackward>,
      BTreeMap<GraphEdgeKey, DenseEdgeForward>,
      BTreeMap<GraphEdgeKey, DenseEdgeEstimate>,
      LogLh,
    ),
    Report,
  > {
    let (node_states, backward) = self.marginal_backward(graph, branch_lengths, &node_states)?;
    let root_key = graph.get_exactly_one_root()?.read_arc().key();
    let log_lh = self.get_log_lh(&node_states, root_key);
    let (node_states, forward, estimates) = self.marginal_forward(graph, branch_lengths, &node_states, &backward)?;
    Ok((node_states, backward, forward, estimates, log_lh))
  }

  pub fn count_transitions(
    &self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
    backward: &BTreeMap<GraphEdgeKey, DenseEdgeBackward>,
    forward: &BTreeMap<GraphEdgeKey, DenseEdgeForward>,
  ) -> Result<MutationCounts, Report> {
    count_transitions_dense(&self.inputs, graph, branch_lengths, node_states, backward, forward)
  }

  pub fn get_log_lh(&self, node_states: &BTreeMap<GraphNodeKey, DenseNodeState>, node_key: GraphNodeKey) -> LogLh {
    node_states
      .get(&node_key)
      .map_or(LogLh::ZERO, |node| node.profile.log_lh)
  }
}

impl HasGtr for PartitionMarginalDiscrete {
  fn gtr(&self) -> &GTR {
    &self.inputs.gtr
  }
  fn gtr_mut(&mut self) -> &mut GTR {
    &mut self.inputs.gtr
  }
  fn sequence_length(&self) -> usize {
    1
  }
}
