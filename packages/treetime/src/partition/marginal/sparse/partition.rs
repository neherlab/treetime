use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::sample::{Resolve, SampleMode};
use crate::ancestral::tip_states::TipStates;
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::MutationCounts;
use crate::make_error;
use crate::partition::marginal::shared::update::{MarginalBackward, MarginalEdges, MarginalForward, MarginalPasses};
use crate::partition::marginal::sparse::count::count_transitions_sparse;
use crate::partition::marginal::sparse::reconstruct::{map_seq, map_seq_sampled, reconstruct_leaf_sequence};
use crate::partition::marginal::sparse::{backward, forward};
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::partition::storage::sparse::{
  SparseEdgeBackward, SparseEdgeForward, SparseEdgeObs, SparseNodeObs, SparseNodeState,
};
use crate::seq::mutation::Sub;
use eyre::Report;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::GraphNodeForward;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{Seq, seq};
use treetime_utils::collections::container::get_exactly_one;
use treetime_utils::interval::range_union::range_union;

#[derive(Clone, Debug, Serialize)]
pub struct PartitionMarginalSparse {
  pub index: usize,
  pub alphabet: Alphabet,
  pub length: usize,
  pub root_sequence: Seq,
  pub obs_nodes: BTreeMap<GraphNodeKey, SparseNodeObs>,
  pub obs_edges: BTreeMap<GraphEdgeKey, SparseEdgeObs>,
}

impl PartitionMarginalSparse {
  pub fn get_sequence_length(&self) -> usize {
    self.length
  }

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

  pub(crate) fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<crate::seq::indel::InDel> {
    self.obs_edges[&edge_key].indels.clone()
  }

  pub(crate) fn root_sequence(&self) -> Seq {
    self.root_sequence.clone()
  }

  pub(crate) fn node_sequence(
    &self,
    node_states: &BTreeMap<GraphNodeKey, SparseNodeState>,
    node_key: GraphNodeKey,
  ) -> Seq {
    match node_states.get(&node_key) {
      Some(node) => node.emitted.clone().unwrap_or_else(|| map_seq(node, &self.alphabet)),
      None => seq![],
    }
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
  ) -> Seq {
    self.node_sequence(node_states, node_key)
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

  pub(crate) fn advance_node_state(
    &self,
    node_states: &mut BTreeMap<GraphNodeKey, SparseNodeState>,
    forward: &BTreeMap<GraphEdgeKey, SparseEdgeForward>,
    node: &GraphNodeForward,
    tips: TipStates,
    sample_mode: SampleMode,
    rng: &mut dyn rand::RngCore,
  ) -> Option<()> {
    let (parent_state, msg_from_parent) = if node.is_root {
      (None, None)
    } else {
      let (parent_key, edge_key) = get_exactly_one(&node.parent_keys).ok()?;
      (
        node_states.get(parent_key).cloned(),
        forward.get(edge_key).map(|f| f.msg_from_parent.clone()),
      )
    };

    let node_obs = self.obs_nodes.get(&node.key)?;
    let node_data = node_states.get(&node.key)?;
    let sample = sample_mode.samples_node(node.is_root);
    if node.is_leaf {
      let seq = reconstruct_leaf_sequence(
        node_data,
        node_obs,
        msg_from_parent.as_ref(),
        parent_state.as_ref(),
        tips.impute,
        &self.alphabet,
      );
      node_states.get_mut(&node.key)?.emitted = Some(seq);
    } else if sample {
      let seq = map_seq_sampled(node_data, &self.alphabet, &mut Resolve::Sample(rng));
      node_states.get_mut(&node.key)?.emitted = Some(seq);
    }

    if !tips.include_leaves && node.is_leaf {
      return None;
    }

    Some(())
  }

  pub(crate) fn reconstruct_node_sequence(
    &self,
    node_states: &mut BTreeMap<GraphNodeKey, SparseNodeState>,
    forward: &BTreeMap<GraphEdgeKey, SparseEdgeForward>,
    node: &GraphNodeForward,
    tips: TipStates,
    sample_mode: SampleMode,
    rng: &mut dyn rand::RngCore,
  ) -> Option<Seq> {
    self.advance_node_state(node_states, forward, node, tips, sample_mode, rng)?;
    Some(self.node_sequence(node_states, node.key))
  }
}

impl MarginalPasses for PartitionMarginalSparse {
  type Node = SparseNodeState;
  type Backward = SparseEdgeBackward;
  type Forward = SparseEdgeForward;
  type Estimate = Vec<Sub>;

  fn marginal_backward(
    &self,
    gtr: &GTR,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: &BTreeMap<GraphNodeKey, SparseNodeState>,
  ) -> Result<MarginalBackward<SparseNodeState, SparseEdgeBackward>, Report> {
    backward::process_backward_indexed(self, gtr, graph, branch_lengths, node_states)
  }

  fn marginal_forward(
    &self,
    gtr: &GTR,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: &BTreeMap<GraphNodeKey, SparseNodeState>,
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
