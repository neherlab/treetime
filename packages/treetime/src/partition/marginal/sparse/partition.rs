use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::sample::SampleMode;
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::MutationCounts;
use crate::make_error;
use crate::partition::marginal::shared::update::{MarginalBackward, MarginalEdges, MarginalForward, MarginalUpdate};
use crate::partition::marginal::sparse::count::count_transitions_sparse;
use crate::partition::marginal::sparse::reconstruct::{map_seq, map_seq_sampled, reconstruct_leaf_sequence};
use crate::partition::marginal::sparse::{backward, forward};
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::partition::storage::sparse::{
  SparseEdgeBackward, SparseEdgeForward, SparseEdgeObs, SparseNodeObs, SparseNodeState,
};
use crate::partition::traits::{HasGtr, PartitionBranchOps, PartitionOptimizeOps};
use crate::seq::mutation::Sub;
use eyre::Report;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::GraphNodeForward;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{LogLh, Seq, seq};
use treetime_utils::collections::container::get_exactly_one;
use treetime_utils::interval::range_union::range_union;

/// The sparse marginal representation as durable, borrowed inputs: the substitution model, the
/// alphabet/length/root metadata, and the Fitch-compressed per-node and per-edge observations. The
/// stage-filled node states, backward/forward messages, and edge estimates are owned separately by the
/// values the passes return.
#[derive(Clone, Debug, Serialize)]
pub struct PartitionMarginalSparse {
  pub index: usize,
  pub gtr: GTR,
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

  pub fn gtr(&self) -> &GTR {
    &self.gtr
  }

  pub fn gtr_mut(&mut self) -> &mut GTR {
    &mut self.gtr
  }

  pub fn weighted_rate(&self) -> f64 {
    self.length as f64 * self.gtr.mu
  }

  pub fn normalize_rate(&mut self, scale: f64) {
    self.gtr.mu /= scale;
  }

  pub fn edge_subs(
    &self,
    estimates: &BTreeMap<GraphEdgeKey, Vec<Sub>>,
    edge_key: GraphEdgeKey,
  ) -> Result<Vec<Sub>, Report> {
    match estimates.get(&edge_key) {
      Some(subs) => Ok(subs.clone()),
      None => make_error!("edge_subs() called before marginal inference populated subs_ml for edge {edge_key:?}"),
    }
  }

  pub fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<crate::seq::indel::InDel> {
    self.obs_edges[&edge_key].indels.clone()
  }

  pub fn root_sequence(&self) -> Seq {
    self.root_sequence.clone()
  }

  pub fn node_sequence(&self, node_states: &BTreeMap<GraphNodeKey, SparseNodeState>, node_key: GraphNodeKey) -> Seq {
    match node_states.get(&node_key) {
      Some(node) => node.emitted.clone().unwrap_or_else(|| map_seq(node, &self.alphabet)),
      None => seq![],
    }
  }

  pub fn edge_effective_length(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<usize, Report> {
    let (parent_key, child_key) = graph.edge_endpoints(edge_key)?;
    let parent_non_char = &self.obs_nodes[&parent_key].non_char;
    let child_non_char = &self.obs_nodes[&child_key].non_char;

    // non_char covers both gaps and unknowns (positions that do not evolve
    // under the substitution model). For internal nodes, non_char is the
    // intersection of children's non_char (Fitch backward pass), so a
    // position is excluded only when all descendants lack data there.
    let non_char_positions: usize = range_union(&[parent_non_char.clone(), child_non_char.clone()])
      .iter()
      .map(|(start, end)| end - start)
      .sum();

    Ok(self.length.saturating_sub(non_char_positions))
  }

  pub fn create_edge_contribution(
    &self,
    backward: &BTreeMap<GraphEdgeKey, SparseEdgeBackward>,
    forward: &BTreeMap<GraphEdgeKey, SparseEdgeForward>,
    edge_key: GraphEdgeKey,
  ) -> Result<OptimizationContribution, Report> {
    OptimizationContribution::from_sparse(
      &self.gtr,
      &backward[&edge_key],
      &forward[&edge_key],
      &self.obs_edges[&edge_key],
    )
  }

  pub fn edge_indel_count(&self, edge_key: GraphEdgeKey) -> usize {
    self.obs_edges[&edge_key].indels.len()
  }

  pub fn extract_ancestral_sequence(
    &self,
    node_states: &BTreeMap<GraphNodeKey, SparseNodeState>,
    node_key: GraphNodeKey,
  ) -> Seq {
    self.node_sequence(node_states, node_key)
  }

  /// Ensure the observation maps have entries for all nodes and edges in the graph, dropping stale
  /// entries. Placeholder observations are created for nodes and edges introduced by topology edits.
  pub fn reconcile_topology(&mut self, graph: &Graph) {
    let graph_node_keys: BTreeSet<GraphNodeKey> = graph.get_nodes().into_iter().map(|n| n.read_arc().key()).collect();
    let graph_edge_keys: BTreeSet<GraphEdgeKey> = graph.get_edges().into_iter().map(|e| e.read_arc().key()).collect();

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

  pub fn reconstruct_node_sequence(
    &self,
    node_states: &mut BTreeMap<GraphNodeKey, SparseNodeState>,
    forward: &BTreeMap<GraphEdgeKey, SparseEdgeForward>,
    node: &GraphNodeForward,
    include_leaves: bool,
    impute: bool,
    sample_mode: SampleMode,
    rng: &mut dyn rand::RngCore,
  ) -> Option<Seq> {
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
    let seq = if node.is_leaf {
      reconstruct_leaf_sequence(
        node_data,
        node_obs,
        msg_from_parent.as_ref(),
        parent_state.as_ref(),
        impute,
        &self.alphabet,
      )
    } else {
      map_seq_sampled(node_data, &self.alphabet, sample, rng)
    };

    // Record the result only where the accessors cannot derive it again: a posterior draw, or a tip,
    // whose observed ambiguity and optional imputation are not a function of the parsimony chain and
    // the posterior. Everything else stays derivable, keeping one source of truth.
    if sample || node.is_leaf {
      node_states.get_mut(&node.key)?.emitted = Some(seq.clone());
    }

    // A suppressed tip is still reconstructed above (so the node-data serializer reads the corrected
    // sequence), but is not emitted to the reconstructed-FASTA visitor.
    if !include_leaves && node.is_leaf {
      return None;
    }

    Some(seq)
  }
}

/// The marginal passes and the update built from them, over the sparse role-typed per-node and
/// per-edge element types. The two passes and the transition counts are the representation-specific
/// operations; the full update and the log-likelihood reads and reset are the same four-step skeleton
/// every representation runs. Stable inputs (`graph`, `branch_lengths`) are borrowed and every pass
/// returns new owned maps, so a failed pass leaves its inputs intact.
impl PartitionMarginalSparse {
  /// Run the marginal backward pass (children before parent).
  pub fn marginal_backward(
    &self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: &BTreeMap<GraphNodeKey, SparseNodeState>,
  ) -> Result<MarginalBackward<SparseNodeState, SparseEdgeBackward>, Report> {
    backward::process_backward_indexed(self, graph, branch_lengths, node_states)
  }

  /// Run the marginal forward pass (parent before children) over the backward messages.
  pub fn marginal_forward(
    &self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: &BTreeMap<GraphNodeKey, SparseNodeState>,
    backward: &BTreeMap<GraphEdgeKey, SparseEdgeBackward>,
  ) -> Result<MarginalForward<SparseNodeState, SparseEdgeForward, Vec<Sub>>, Report> {
    forward::process_forward_indexed(self, graph, branch_lengths, node_states, backward)
  }

  /// Count posterior-weighted state transitions over the tree, the input GTR inference reads.
  pub fn count_transitions(
    &self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    node_states: &BTreeMap<GraphNodeKey, SparseNodeState>,
    backward: &BTreeMap<GraphEdgeKey, SparseEdgeBackward>,
    forward: &BTreeMap<GraphEdgeKey, SparseEdgeForward>,
  ) -> Result<MutationCounts, Report> {
    count_transitions_sparse(
      &self.gtr,
      self.length,
      graph,
      branch_lengths,
      node_states,
      backward,
      forward,
    )
  }

  /// Run a full marginal update (backward, then forward) over the given node states.
  ///
  /// The substitution log likelihood is read at the root between the two passes: after the backward
  /// pass the root profile holds the likelihood of the observed data under the model, while the forward
  /// pass overwrites every node profile with its posterior.
  pub fn marginal_update(
    &self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: BTreeMap<GraphNodeKey, SparseNodeState>,
  ) -> Result<MarginalUpdate<SparseNodeState, SparseEdgeBackward, SparseEdgeForward, Vec<Sub>>, Report> {
    let MarginalBackward { node_states, backward } = self.marginal_backward(graph, branch_lengths, &node_states)?;
    let log_lh = self.root_log_lh(graph, &node_states)?;
    let MarginalForward {
      node_states,
      forward,
      estimates,
    } = self.marginal_forward(graph, branch_lengths, &node_states, &backward)?;
    Ok(MarginalUpdate {
      node_states,
      edges: MarginalEdges {
        backward,
        forward,
        estimates,
      },
      log_lh,
    })
  }

  /// The profile log likelihood recorded for one node, or zero when the node has no state.
  pub fn get_log_lh(&self, node_states: &BTreeMap<GraphNodeKey, SparseNodeState>, node_key: GraphNodeKey) -> LogLh {
    node_states
      .get(&node_key)
      .map_or(LogLh::ZERO, |node| node.profile.log_lh)
  }

  /// The profile log likelihood recorded at the tree root, or zero when the root has no state.
  pub fn root_log_lh(
    &self,
    graph: &Graph,
    node_states: &BTreeMap<GraphNodeKey, SparseNodeState>,
  ) -> Result<LogLh, Report> {
    let root_key = graph.get_exactly_one_root()?.read_arc().key();
    Ok(self.get_log_lh(node_states, root_key))
  }

  /// Zero every node's profile log likelihood before a backward pass whose result is read as a
  /// likelihood.
  ///
  /// The backward pass reads a leaf's profile log likelihood as its message to the parent and folds it
  /// up the tree into the root log likelihood. A leaf profile carried over from an earlier forward pass
  /// holds that pass's posterior log likelihood, which adds a large model-independent constant to the
  /// result. That constant does not move the optimum of a rate search in exact arithmetic, but its
  /// magnitude erodes precision in Brent's parabolic interpolation and shifts the selected rate.
  pub fn reset_node_log_lh(&self, node_states: &mut BTreeMap<GraphNodeKey, SparseNodeState>) {
    for node in node_states.values_mut() {
      node.profile.log_lh = LogLh::ZERO;
    }
  }
}

/// The per-edge results one sparse marginal update returns.
pub type SparseMarginalEdges = MarginalEdges<SparseEdgeBackward, SparseEdgeForward, Vec<Sub>>;

/// Short-lived read view over a completed sparse reconstruction: borrows the durable partition inputs
/// (including the Fitch observations) together with the node states and the per-edge results the
/// passes returned. Assembled at a consumer boundary purely to read; never stored.
pub struct SparseReadout<'a> {
  pub partition: &'a PartitionMarginalSparse,
  pub node_states: &'a BTreeMap<GraphNodeKey, SparseNodeState>,
  pub edges: &'a SparseMarginalEdges,
}

impl PartitionBranchOps for SparseReadout<'_> {
  fn sequence_length(&self) -> usize {
    self.partition.length
  }

  fn edge_subs(&self, _graph: &Graph, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    self.partition.edge_subs(&self.edges.estimates, edge_key)
  }

  fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<crate::seq::indel::InDel> {
    self.partition.edge_indels(edge_key)
  }

  fn root_sequence(&self, _graph: &Graph) -> Result<Seq, Report> {
    Ok(self.partition.root_sequence())
  }

  fn node_sequence(&self, node_key: GraphNodeKey) -> Seq {
    self.partition.node_sequence(self.node_states, node_key)
  }

  fn edge_effective_length(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<usize, Report> {
    self.partition.edge_effective_length(graph, edge_key)
  }
}

impl PartitionOptimizeOps for SparseReadout<'_> {
  fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report> {
    self
      .partition
      .create_edge_contribution(&self.edges.backward, &self.edges.forward, edge_key)
  }

  fn edge_indel_count(&self, edge_key: GraphEdgeKey) -> usize {
    self.partition.edge_indel_count(edge_key)
  }
}

impl HasGtr for PartitionMarginalSparse {
  fn gtr(&self) -> &GTR {
    &self.gtr
  }
  fn gtr_mut(&mut self) -> &mut GTR {
    &mut self.gtr
  }
  fn sequence_length(&self) -> usize {
    self.length
  }
}
