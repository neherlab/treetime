use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::sample::{SampleMode, resolve_profile};
use crate::constants::MIN_BRANCH_LENGTH_FRACTION;
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::MutationCounts;
use crate::make_report;
use crate::partition::marginal::shared::data::{DenseInputs, count_transitions_dense};
use crate::partition::marginal::shared::pass::{IndexedKind, indexed_backward, indexed_forward};
use crate::partition::marginal::shared::update::{MarginalBackward, MarginalEdges, MarginalForward, MarginalPasses};
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::partition::storage::dense::{
  DenseEdgeBackward, DenseEdgeEstimate, DenseEdgeForward, DenseNodeState, DenseSeqDistribution,
};
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::izip;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::GraphNodeForward;
use treetime_graph::node::GraphNodeKey;
use treetime_io::nwk::NwkFastaNodeInput;
use treetime_primitives::{Seq, seq};
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
  pub fn new(index: usize, alphabet: Alphabet, length: usize) -> Self {
    let min_branch_length = MIN_BRANCH_LENGTH_FRACTION / length as f64;
    Self {
      inputs: DenseInputs {
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

  pub fn get_sequence_length(&self) -> usize {
    self.length
  }

  /// Build the initial dense node states by attaching each leaf's observed sequence. Internal-node
  /// entries are created lazily by the backward pass; this returns the leaf-seeded node-state map.
  pub fn attach_sequences(
    &self,
    graph: &Graph,
    node_inputs: &BTreeMap<GraphNodeKey, NwkFastaNodeInput>,
  ) -> Result<BTreeMap<GraphNodeKey, DenseNodeState>, Report> {
    let mut node_states = BTreeMap::new();
    for leaf in graph.get_leaves() {
      let leaf_key = leaf.read_arc().key();
      let node = &node_inputs[&leaf_key];

      let seq = node
        .aln
        .as_ref()
        .ok_or_else(|| make_report!("Leaf sequence not found: '{}'", node.name.as_deref().unwrap_or("")))?;

      node_states.insert(leaf_key, DenseNodeState::new(seq, &self.alphabet)?);
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
    gtr: &GTR,
    backward: &BTreeMap<GraphEdgeKey, DenseEdgeBackward>,
    forward: &BTreeMap<GraphEdgeKey, DenseEdgeForward>,
    edge_key: GraphEdgeKey,
  ) -> OptimizationContribution {
    OptimizationContribution::from_dense(gtr, &backward[&edge_key], &forward[&edge_key])
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

impl MarginalPasses for PartitionMarginalDense {
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
      Some(&self.alphabet),
      self.length,
      IndexedKind::Dense,
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
      Some(&self.alphabet),
      IndexedKind::Dense,
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
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
    backward: &BTreeMap<GraphEdgeKey, DenseEdgeBackward>,
    forward: &BTreeMap<GraphEdgeKey, DenseEdgeForward>,
  ) -> Result<MutationCounts, Report> {
    count_transitions_dense(&self.inputs, gtr, graph, branch_lengths, node_states, backward, forward)
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
