use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::sample::{SampleMode, resolve_profile};
use crate::constants::MIN_BRANCH_LENGTH_FRACTION;
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::MutationCounts;
use crate::make_report;
use crate::partition::marginal::shared::data::{DenseInputs, count_transitions_dense};
use crate::partition::marginal::shared::pass::{IndexedKind, indexed_backward, indexed_forward};
use crate::partition::marginal::shared::update::{
  MarginalBackward, MarginalEdges, MarginalForward, MarginalStates, MarginalUpdate,
};
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::partition::storage::dense::{
  DenseEdgeBackward, DenseEdgeEstimate, DenseEdgeForward, DenseNodeState, DenseSeqDistribution,
};
use crate::partition::traits::{HasGtr, PartitionBranchOps, PartitionOptimizeOps};
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::izip;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::GraphNodeForward;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::FastaRecord;
use treetime_primitives::{LogLh, Seq, seq};
use treetime_utils::array::ndarray::argmax_first;
use treetime_utils::interval::range::range_contains;
use treetime_utils::interval::range_union::range_union;

/// The dense marginal representation as durable, borrowed inputs: the substitution model and the
/// alphabet/length metadata. The stage-filled node states, backward/forward messages, and edge
/// estimates are owned separately by the values the passes return.
#[derive(Clone, Debug, Serialize)]
pub struct PartitionMarginalDense {
  pub inputs: DenseInputs,
  pub index: usize,
  pub alphabet: Alphabet,
  pub length: usize,
}

impl PartitionMarginalDense {
  pub fn new(index: usize, gtr: GTR, alphabet: Alphabet, length: usize) -> Self {
    let min_branch_length = MIN_BRANCH_LENGTH_FRACTION / length as f64;
    Self {
      inputs: DenseInputs {
        gtr,
        min_branch_length,
        // Nucleotide ancestral inference filters signal-free (gap-only) root
        // columns out of the equilibrium-frequency prior.
        filter_uninformative_root: true,
      },
      index,
      alphabet,
      length,
    }
  }

  pub fn gtr(&self) -> &GTR {
    &self.inputs.gtr
  }

  pub fn gtr_mut(&mut self) -> &mut GTR {
    &mut self.inputs.gtr
  }

  pub fn get_sequence_length(&self) -> usize {
    self.length
  }

  pub fn weighted_rate(&self) -> f64 {
    self.length as f64 * self.inputs.gtr.mu
  }

  pub fn normalize_rate(&mut self, scale: f64) {
    self.inputs.gtr.mu /= scale;
  }

  /// Build the initial dense node states by attaching each leaf's observed sequence. Internal-node
  /// entries are created lazily by the backward pass; this returns the leaf-seeded node-state map.
  pub fn attach_sequences(
    &self,
    graph: &Graph,
    aln: &[FastaRecord],
    names: &BTreeMap<GraphNodeKey, Option<String>>,
  ) -> Result<BTreeMap<GraphNodeKey, DenseNodeState>, Report> {
    let aln_by_name = aln.iter().fold(BTreeMap::new(), |mut records, record| {
      records.entry(record.seq_name.as_str()).or_insert(record);
      records
    });
    let mut node_states = BTreeMap::new();
    for leaf in graph.get_leaves() {
      let leaf_key = leaf.read_arc().key();

      let leaf_name = names[&leaf_key].clone().ok_or_else(|| {
        make_report!("Expected all leaf nodes to have names, such that they can be matched to their corresponding sequences. But found a leaf node that has no name.")
      })?;

      let leaf_fasta = aln_by_name
        .get(leaf_name.as_str())
        .copied()
        .ok_or_else(|| make_report!("Leaf sequence not found: '{leaf_name}'"))?;

      node_states.insert(leaf_key, DenseNodeState::new(&leaf_fasta.seq, &self.alphabet)?);
    }
    Ok(node_states)
  }

  pub fn edge_subs(
    &self,
    node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
    graph: &Graph,
    edge_key: GraphEdgeKey,
  ) -> Result<Vec<Sub>, Report> {
    let (parent_key, child_key) = graph.edge_endpoints(edge_key)?;
    let parent_non_char = &node_states[&parent_key].seq.non_char;
    let child_non_char = &node_states[&child_key].seq.non_char;

    let parent_profile = &node_states[&parent_key].profile.dis;
    let child_profile = &node_states[&child_key].profile.dis;
    let mut subs = Vec::new();

    for (pos, parent, child) in izip!(0..parent_profile.nrows(), parent_profile.rows(), child_profile.rows()) {
      let parent_state = self.alphabet.char(argmax_first(&parent).unwrap_or(0));
      let child_state = self.alphabet.char(argmax_first(&child).unwrap_or(0));
      if parent_state == child_state {
        continue;
      }
      if !self.alphabet.is_canonical(parent_state) || !self.alphabet.is_canonical(child_state) {
        continue;
      }
      if range_contains(parent_non_char, pos) || range_contains(child_non_char, pos) {
        continue;
      }

      subs.push(Sub::new(parent_state, pos, child_state)?);
    }

    Ok(subs)
  }

  pub fn edge_indels(
    &self,
    estimates: &BTreeMap<GraphEdgeKey, DenseEdgeEstimate>,
    edge_key: GraphEdgeKey,
  ) -> Vec<crate::seq::indel::InDel> {
    estimates[&edge_key].indels.clone()
  }

  pub fn root_sequence(
    &self,
    node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
    graph: &Graph,
  ) -> Result<Seq, Report> {
    Ok(assign_sequence(&node_states[&graph.root_key()?], &self.alphabet))
  }

  pub fn node_sequence(&self, node_states: &BTreeMap<GraphNodeKey, DenseNodeState>, node_key: GraphNodeKey) -> Seq {
    assign_sequence(&node_states[&node_key], &self.alphabet)
  }

  pub fn edge_effective_length(
    &self,
    node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
    graph: &Graph,
    edge_key: GraphEdgeKey,
  ) -> Result<usize, Report> {
    let (parent_key, child_key) = graph.edge_endpoints(edge_key)?;
    let parent_non_char = &node_states[&parent_key].seq.non_char;
    let child_non_char = &node_states[&child_key].seq.non_char;

    let non_char_positions: usize = range_union(&[parent_non_char.clone(), child_non_char.clone()])
      .iter()
      .map(|(start, end)| end - start)
      .sum();

    Ok(self.length.saturating_sub(non_char_positions))
  }

  pub fn create_edge_contribution(
    &self,
    backward: &BTreeMap<GraphEdgeKey, DenseEdgeBackward>,
    forward: &BTreeMap<GraphEdgeKey, DenseEdgeForward>,
    edge_key: GraphEdgeKey,
  ) -> OptimizationContribution {
    OptimizationContribution::from_dense(&self.inputs.gtr, &backward[&edge_key], &forward[&edge_key])
  }

  pub fn edge_indel_count(
    &self,
    estimates: &BTreeMap<GraphEdgeKey, DenseEdgeEstimate>,
    edge_key: GraphEdgeKey,
  ) -> usize {
    estimates[&edge_key].indels.len()
  }

  pub fn extract_ancestral_sequence(
    &self,
    node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
    node_key: GraphNodeKey,
  ) -> Seq {
    if let Some(seq_info) = node_states.get(&node_key) {
      // Convergence checks and other intermediate reads always use the deterministic most-likely
      // state, independent of the user's output sampling mode.
      assign_sequence(seq_info, &self.alphabet)
    } else {
      seq! {}
    }
  }

  /// Reconstruct the sequence for one node, recording it into the node state so the node-data serializer
  /// reads back the flag-aware sequence. Returns `None` for a suppressed tip.
  pub fn reconstruct_node_sequence(
    &self,
    node_states: &mut BTreeMap<GraphNodeKey, DenseNodeState>,
    node: &GraphNodeForward,
    include_leaves: bool,
    impute: bool,
    sample_mode: SampleMode,
    rng: &mut dyn rand::RngCore,
  ) -> Option<Seq> {
    let seq = {
      let seq_info = node_states.get(&node.key)?;
      if node.is_leaf {
        // The dense tip keeps its observed input (gaps and unknowns already stamped), so it never
        // suffers the sparse tip corruption. Imputation resolves ambiguous/unknown positions (N and
        // IUPAC codes, not gaps) to the argmax of the leaf marginal posterior. The forward pass folds
        // the observed ambiguity mask into the leaf profile, so its argmax is the parent-informed
        // most likely state, matching v0.
        let mut seq = seq_info.seq.sequence.clone();
        if impute && seq_info.profile.dis.nrows() == seq.len() {
          for pos in 0..seq.len() {
            let ch = seq[pos];
            if !self.alphabet.is_canonical(ch) && !self.alphabet.is_gap(ch) {
              if let Some(idx) = argmax_first(&seq_info.profile.dis.row(pos)) {
                seq[pos] = self.alphabet.char(idx);
              }
            }
          }
        }
        seq
      } else {
        let sample = sample_mode.samples_node(node.is_root);
        assign_sequence_sampled(seq_info, &self.alphabet, sample, rng)
      }
    };

    // Persist the reconstruction so the node-data serializer reads back this flag-aware sequence
    // (the dense augur path returns `seq.sequence`), keeping the JSON and the reconstructed FASTA
    // consistent and matching the sparse backend.
    if let Some(node_data) = node_states.get_mut(&node.key) {
      node_data.seq.sequence = seq.clone();
    }

    // A suppressed tip is still reconstructed above (so the node-data serializer reads the corrected
    // sequence), but is not emitted to the reconstructed-FASTA visitor.
    if !include_leaves && node.is_leaf {
      return None;
    }

    Some(seq)
  }
}

/// The marginal passes and the update built from them, over the dense role-typed per-node and per-edge
/// element types. The two passes and the transition counts are the representation-specific operations;
/// the full update and the log-likelihood reads and reset are the same four-step skeleton every
/// representation runs. Stable inputs (`graph`, `branch_lengths`) are borrowed and every pass returns
/// new owned maps, so a failed pass leaves its inputs intact.
impl PartitionMarginalDense {
  /// Run the marginal backward pass (children before parent).
  pub fn marginal_backward(
    &self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
  ) -> Result<MarginalBackward<DenseNodeState, DenseEdgeBackward>, Report> {
    indexed_backward(
      &self.inputs,
      Some(&self.alphabet),
      self.length,
      IndexedKind::Dense,
      graph,
      branch_lengths,
      node_states,
    )
  }

  /// Run the marginal forward pass (parent before children) over the backward messages.
  pub fn marginal_forward(
    &self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
    backward: &BTreeMap<GraphEdgeKey, DenseEdgeBackward>,
  ) -> Result<MarginalForward<DenseNodeState, DenseEdgeForward, DenseEdgeEstimate>, Report> {
    indexed_forward(
      &self.inputs,
      Some(&self.alphabet),
      IndexedKind::Dense,
      graph,
      branch_lengths,
      node_states,
      backward,
    )
  }

  /// Count posterior-weighted state transitions over the tree, the input GTR inference reads.
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

  /// Run a full marginal update (backward, then forward) over the given node states.
  ///
  /// The substitution log likelihood is read at the root between the two passes: after the backward
  /// pass the root profile holds the likelihood of the observed data under the model, while the forward
  /// pass overwrites every node profile with its posterior.
  pub fn marginal_update(
    &self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: BTreeMap<GraphNodeKey, DenseNodeState>,
  ) -> Result<MarginalUpdate<DenseNodeState, DenseEdgeBackward, DenseEdgeForward, DenseEdgeEstimate>, Report> {
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

  /// Run a full marginal update and return only the refreshed node states and the substitution log
  /// likelihood.
  ///
  /// For callers that read no per-edge value: the passes are identical to [`Self::marginal_update`], and
  /// the per-edge messages and estimates are dropped as the update returns instead of travelling through
  /// the caller.
  pub fn marginal_states(
    &self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: BTreeMap<GraphNodeKey, DenseNodeState>,
  ) -> Result<MarginalStates<DenseNodeState>, Report> {
    let MarginalUpdate {
      node_states, log_lh, ..
    } = self.marginal_update(graph, branch_lengths, node_states)?;
    Ok(MarginalStates { node_states, log_lh })
  }

  /// The profile log likelihood recorded for one node, or zero when the node has no state.
  pub fn get_log_lh(&self, node_states: &BTreeMap<GraphNodeKey, DenseNodeState>, node_key: GraphNodeKey) -> LogLh {
    node_states
      .get(&node_key)
      .map_or(LogLh::ZERO, |node| node.profile.log_lh)
  }

  /// The profile log likelihood recorded at the tree root, or zero when the root has no state.
  pub fn root_log_lh(
    &self,
    graph: &Graph,
    node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
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
  pub fn reset_node_log_lh(&self, node_states: &mut BTreeMap<GraphNodeKey, DenseNodeState>) {
    for node in node_states.values_mut() {
      node.profile.log_lh = LogLh::ZERO;
    }
  }
}

/// Deterministic most-likely-state sequence assignment. Used by the forward pass and convergence
/// reads, which must stay reproducible regardless of the user's output sampling mode.
pub(crate) fn assign_sequence(seq_info: &DenseNodeState, alphabet: &Alphabet) -> Seq {
  assign_sequence_sampled(seq_info, alphabet, false, &mut rand::thread_rng())
}

fn assign_sequence_sampled(
  seq_info: &DenseNodeState,
  alphabet: &Alphabet,
  sample: bool,
  rng: &mut dyn rand::RngCore,
) -> Seq {
  let mut seq = prof2seq_sampled(&seq_info.profile, alphabet, sample, rng);
  for gap in &seq_info.seq.gaps {
    seq[gap.0..gap.1].fill(alphabet.gap());
  }
  for unk in &seq_info.seq.unknown {
    seq[unk.0..unk.1].fill(alphabet.unknown());
  }
  seq
}

fn prof2seq_sampled(
  profile: &DenseSeqDistribution,
  alphabet: &Alphabet,
  sample: bool,
  rng: &mut dyn rand::RngCore,
) -> Seq {
  let mut seq = seq! {};
  for row in profile.dis.rows() {
    seq.push(alphabet.char(resolve_profile(row, sample, rng)));
  }
  seq
}

/// The per-edge results one dense marginal update returns.
pub type DenseMarginalEdges = MarginalEdges<DenseEdgeBackward, DenseEdgeForward, DenseEdgeEstimate>;

/// Short-lived read view over a completed dense reconstruction: borrows the durable partition inputs
/// together with the node states and the per-edge results the passes returned. Assembled at a
/// consumer boundary purely to read; never stored.
pub struct DenseReadout<'a> {
  pub partition: &'a PartitionMarginalDense,
  pub node_states: &'a BTreeMap<GraphNodeKey, DenseNodeState>,
  pub edges: &'a DenseMarginalEdges,
}

impl PartitionBranchOps for DenseReadout<'_> {
  fn sequence_length(&self) -> usize {
    self.partition.length
  }

  fn edge_subs(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    self.partition.edge_subs(self.node_states, graph, edge_key)
  }

  fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<crate::seq::indel::InDel> {
    self.partition.edge_indels(&self.edges.estimates, edge_key)
  }

  fn root_sequence(&self, graph: &Graph) -> Result<Seq, Report> {
    self.partition.root_sequence(self.node_states, graph)
  }

  fn node_sequence(&self, node_key: GraphNodeKey) -> Seq {
    self.partition.node_sequence(self.node_states, node_key)
  }

  fn edge_effective_length(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<usize, Report> {
    self.partition.edge_effective_length(self.node_states, graph, edge_key)
  }
}

impl PartitionOptimizeOps for DenseReadout<'_> {
  fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report> {
    Ok(
      self
        .partition
        .create_edge_contribution(&self.edges.backward, &self.edges.forward, edge_key),
    )
  }

  fn edge_indel_count(&self, edge_key: GraphEdgeKey) -> usize {
    self.partition.edge_indel_count(&self.edges.estimates, edge_key)
  }
}

impl HasGtr for PartitionMarginalDense {
  fn gtr(&self) -> &GTR {
    &self.inputs.gtr
  }
  fn gtr_mut(&mut self) -> &mut GTR {
    &mut self.inputs.gtr
  }
  fn sequence_length(&self) -> usize {
    self.length
  }
}
