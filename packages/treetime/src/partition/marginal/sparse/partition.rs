use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::MutationCounts;
use crate::partition::marginal::sample::{Resolve, SampleMode};
use crate::partition::marginal::sequences::sample_internal_sequences;
use crate::partition::marginal::shared::update::{MarginalBackward, MarginalEdges, MarginalForward, MarginalPasses};
use crate::partition::marginal::sparse::count::count_transitions_sparse;
use crate::partition::marginal::sparse::reconstruct::{
  map_seq, map_seq_overlay, map_seq_sampled, reconstruct_leaf_sequence,
};
use crate::partition::marginal::sparse::{backward, forward};
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::partition::storage::sparse::{
  SparseEdgeBackward, SparseEdgeForward, SparseEdgeObs, SparseNodeObs, SparseNodeState,
};
use crate::seq::mutation::{Mutation, MutationTrack, Sub, combine_edge_mutations};
use crate::seq::overlay::SeqOverlay;
use crate::{make_error, make_internal_report};
use eyre::{Report, WrapErr};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{LogLh, Seq, seq};
use treetime_utils::interval::range_union::range_union;

#[derive(Clone, Debug, Serialize, deser::Serialize)]
pub struct PartitionMarginalSparse {
  pub index: usize,
  pub alphabet: Alphabet,
  pub length: usize,
  pub root_sequence: Seq,
  pub obs_nodes: BTreeMap<GraphNodeKey, SparseNodeObs>,
  pub obs_edges: BTreeMap<GraphEdgeKey, SparseEdgeObs>,
}

impl PartitionMarginalSparse {
  pub(crate) fn edge_subs(
    &self,
    estimates: &BTreeMap<GraphEdgeKey, Vec<Sub>>,
    edge_key: GraphEdgeKey,
  ) -> Result<Vec<Sub>, Report> {
    match estimates.get(&edge_key) {
      Some(subs) => Ok(subs.clone()),
      None => make_error!("edge_subs() called before marginal inference populated subs_ml for edge {edge_key:?}"),
    }
  }

  pub fn edge_fitch_mutations(&self, edge_key: GraphEdgeKey, track: &MutationTrack) -> Result<Vec<Mutation>, Report> {
    let edge = &self.obs_edges[&edge_key];
    combine_edge_mutations(edge.fitch_subs().to_vec(), &edge.indels, track)
  }

  pub(crate) fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<crate::seq::indel::InDel> {
    self.obs_edges[&edge_key].indels.clone()
  }

  pub fn fitch_root_sequence(&self) -> Seq {
    self.root_sequence.clone()
  }

  pub(crate) fn root_sequence(
    &self,
    node_states: &BTreeMap<GraphNodeKey, SparseNodeState>,
    graph: &Graph,
  ) -> Result<Seq, Report> {
    Ok(map_seq(&node_states[&graph.root_key()?], &self.alphabet))
  }

  pub(crate) fn edge_effective_length(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<usize, Report> {
    let (parent_key, child_key) = graph.edge_endpoints(edge_key)?;
    let parent_non_char = &self.obs_nodes[&parent_key].non_char;
    let child_non_char = &self.obs_nodes[&child_key].non_char;

    let non_char_positions: usize = range_union(&[parent_non_char.clone(), child_non_char.clone()])
      .iter()
      .map(|(start, end)| end - start)
      .sum();

    Ok(self.length.saturating_sub(non_char_positions))
  }

  pub(crate) fn create_edge_contribution(
    &self,
    gtr: &GTR,
    backward: &BTreeMap<GraphEdgeKey, SparseEdgeBackward>,
    forward: &BTreeMap<GraphEdgeKey, SparseEdgeForward>,
    edge_key: GraphEdgeKey,
  ) -> Result<OptimizationContribution, Report> {
    OptimizationContribution::from_sparse(
      gtr,
      &backward[&edge_key],
      &forward[&edge_key],
      &self.obs_edges[&edge_key],
    )
  }

  pub(crate) fn edge_indel_count(&self, edge_key: GraphEdgeKey) -> usize {
    self.obs_edges[&edge_key].indels.len()
  }

  pub(crate) fn extract_ancestral_sequence(
    &self,
    node_states: &BTreeMap<GraphNodeKey, SparseNodeState>,
    node_key: GraphNodeKey,
  ) -> SeqOverlay {
    node_states.get(&node_key).map_or_else(
      || SeqOverlay::from(seq![]),
      |node| map_seq_overlay(node, &self.alphabet),
    )
  }

  pub(crate) fn reconcile_topology(&mut self, graph: &Graph) {
    let graph_node_keys: BTreeSet<GraphNodeKey> = graph.get_nodes().map(|n| n.key()).collect();
    let graph_edge_keys: BTreeSet<GraphEdgeKey> = graph.get_edges().map(|e| e.key()).collect();

    for &key in &graph_node_keys {
      self
        .obs_nodes
        .entry(key)
        .or_insert_with(|| SparseNodeObs::empty(&self.alphabet));
    }
    for &key in &graph_edge_keys {
      self.obs_edges.entry(key).or_default();
    }
    self.obs_nodes.retain(|k, _| graph_node_keys.contains(k));
    self.obs_edges.retain(|k, _| graph_edge_keys.contains(k));
  }

  pub(crate) fn node_sequence(
    &self,
    graph: &Graph,
    node_states: &BTreeMap<GraphNodeKey, SparseNodeState>,
    forward: &BTreeMap<GraphEdgeKey, SparseEdgeForward>,
    impute: bool,
    node_key: GraphNodeKey,
  ) -> Result<Seq, Report> {
    let node_data = &node_states[&node_key];
    let is_leaf = graph
      .get_node(node_key)
      .ok_or_else(|| make_internal_report!("Node {node_key} not found while reconstructing its sequence"))?
      .is_leaf();
    if !is_leaf {
      return Ok(map_seq(node_data, &self.alphabet));
    }
    let parent = graph
      .node_parent(node_key)
      .wrap_err_with(|| format!("When reconstructing the sequence of node {node_key}"))?;
    let (parent_state, msg_from_parent) = match parent {
      None => (None, None),
      Some((parent_key, edge_key)) => (
        Some(&node_states[&parent_key]),
        Some(&forward[&edge_key].msg_from_parent),
      ),
    };
    Ok(reconstruct_leaf_sequence(
      node_data,
      &self.obs_nodes[&node_key],
      msg_from_parent,
      parent_state,
      impute,
      &self.alphabet,
    ))
  }

  pub(crate) fn sample_sequences(
    &self,
    graph: &Graph,
    node_states: &BTreeMap<GraphNodeKey, SparseNodeState>,
    sample_mode: SampleMode,
    rng: &mut dyn rand::RngCore,
  ) -> Result<BTreeMap<GraphNodeKey, Seq>, Report> {
    sample_internal_sequences(graph, sample_mode, |key| {
      map_seq_sampled(&node_states[&key], &self.alphabet, &mut Resolve::Sample(&mut *rng))
    })
  }
}

impl MarginalPasses for PartitionMarginalSparse {
  type Node = SparseNodeState;
  type BackwardInput = BTreeMap<GraphNodeKey, SparseNodeState>;
  type Backward = SparseEdgeBackward;
  type Forward = SparseEdgeForward;
  type Estimate = Vec<Sub>;

  fn backward_input(node_states: &BTreeMap<GraphNodeKey, SparseNodeState>) -> &Self::BackwardInput {
    node_states
  }

  fn fresh_backward_input(node_states: &BTreeMap<GraphNodeKey, SparseNodeState>) -> Self::BackwardInput {
    let mut input = node_states.clone();
    for node in input.values_mut() {
      node.profile.log_lh = LogLh::ZERO;
    }
    input
  }

  fn marginal_backward(
    &self,
    gtr: &GTR,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    input: &Self::BackwardInput,
  ) -> Result<MarginalBackward<SparseNodeState, SparseEdgeBackward>, Report> {
    backward::process_backward_indexed(self, gtr, graph, branch_lengths, input)
  }

  fn marginal_forward(
    &self,
    gtr: &GTR,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: BTreeMap<GraphNodeKey, SparseNodeState>,
    backward: &BTreeMap<GraphEdgeKey, SparseEdgeBackward>,
  ) -> Result<MarginalForward<SparseNodeState, SparseEdgeForward, Vec<Sub>>, Report> {
    forward::process_forward_indexed(self, gtr, graph, branch_lengths, node_states, backward)
  }

  fn count_transitions(
    &self,
    gtr: &GTR,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: &BTreeMap<GraphNodeKey, SparseNodeState>,
    backward: &BTreeMap<GraphEdgeKey, SparseEdgeBackward>,
    forward: &BTreeMap<GraphEdgeKey, SparseEdgeForward>,
  ) -> Result<MutationCounts, Report> {
    count_transitions_sparse(gtr, self.length, graph, branch_lengths, node_states, backward, forward)
  }
}

pub type SparseMarginalEdges = MarginalEdges<SparseEdgeBackward, SparseEdgeForward, Vec<Sub>>;
