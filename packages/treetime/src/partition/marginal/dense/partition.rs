use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::reconstruction::{ReconstructedSequences, reconstruct_preorder};
use crate::ancestral::sample::{Resolve, SampleMode, resolve_profile};
use crate::ancestral::tip_states::TipStates;
use crate::constants::MIN_BRANCH_LENGTH_FRACTION;
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::MutationCounts;
use crate::make_report;
use crate::partition::marginal::shared::data::{DenseInputs, count_transitions_dense};
use crate::partition::marginal::shared::pass::{IndexedKind, indexed_backward, indexed_forward};
use crate::partition::marginal::shared::update::{MarginalBackward, MarginalEdges, MarginalForward, MarginalPasses};
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::partition::storage::dense::{
  DenseEdgeBackward, DenseEdgeEstimate, DenseEdgeForward, DenseLeafObs, DenseNodeState, DenseSeqDistribution,
};
use crate::seq::alignment::{NodeSeqInput, get_common_length_of_node_inputs};
use crate::seq::indel::InDel;
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::izip;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{Seq, seq};
use treetime_utils::array::ndarray::argmax_first;
use treetime_utils::interval::range::range_contains;
use treetime_utils::interval::range_union::range_union;

#[derive(Clone, Debug, Serialize)]
pub struct PartitionMarginalDense {
  pub(crate) inputs: DenseInputs,
  pub(crate) index: usize,
  pub(crate) alphabet: Alphabet,
  pub(crate) length: usize,
  pub(crate) obs_leaves: BTreeMap<GraphNodeKey, Option<DenseLeafObs>>,
}

impl PartitionMarginalDense {
  #[allow(
    clippy::as_conversions,
    reason = "count/index numeric cast is exact for the domain range"
  )]
  pub(crate) fn new(
    index: usize,
    alphabet: Alphabet,
    graph: &Graph,
    node_inputs: &BTreeMap<GraphNodeKey, NodeSeqInput>,
  ) -> Result<Self, Report> {
    let length = get_common_length_of_node_inputs(node_inputs)?;
    let obs_leaves = graph
      .get_leaves()
      .map(|leaf| {
        let leaf_key = leaf.key();
        let node = &node_inputs[&leaf_key];
        let seq = node
          .seq
          .as_ref()
          .ok_or_else(|| make_report!("Leaf sequence not found: '{}'", node.name.as_deref().unwrap_or("")))?;
        Ok((leaf_key, Some(DenseLeafObs::new(seq, &alphabet))))
      })
      .collect::<Result<BTreeMap<_, _>, Report>>()?;
    let min_branch_length = MIN_BRANCH_LENGTH_FRACTION / length as f64;
    Ok(Self {
      inputs: DenseInputs {
        min_branch_length,
        filter_uninformative_root: true,
      },
      index,
      alphabet,
      length,
      obs_leaves,
    })
  }

  pub(crate) fn edge_subs(
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

  pub(crate) fn edge_indels(
    &self,
    estimates: &BTreeMap<GraphEdgeKey, DenseEdgeEstimate>,
    edge_key: GraphEdgeKey,
  ) -> Vec<InDel> {
    estimates[&edge_key].indels.clone()
  }

  pub(crate) fn root_sequence(
    &self,
    node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
    graph: &Graph,
  ) -> Result<Seq, Report> {
    Ok(assign_sequence(&node_states[&graph.root_key()?], &self.alphabet))
  }

  pub(crate) fn edge_effective_length(
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

  pub(crate) fn create_edge_contribution(
    &self,
    gtr: &GTR,
    backward: &BTreeMap<GraphEdgeKey, DenseEdgeBackward>,
    forward: &BTreeMap<GraphEdgeKey, DenseEdgeForward>,
    edge_key: GraphEdgeKey,
  ) -> OptimizationContribution {
    OptimizationContribution::from_dense(gtr, &backward[&edge_key], &forward[&edge_key])
  }

  pub(crate) fn edge_indel_count(
    &self,
    estimates: &BTreeMap<GraphEdgeKey, DenseEdgeEstimate>,
    edge_key: GraphEdgeKey,
  ) -> usize {
    estimates[&edge_key].indels.len()
  }

  pub(crate) fn extract_ancestral_sequence(
    &self,
    node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
    node_key: GraphNodeKey,
  ) -> Seq {
    if let Some(seq_info) = node_states.get(&node_key) {
      assign_sequence(seq_info, &self.alphabet)
    } else {
      seq! {}
    }
  }

  pub(crate) fn reconstruct_sequences(
    &self,
    graph: &Graph,
    node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
    tips: TipStates,
    sample_mode: SampleMode,
    rng: &mut dyn rand::RngCore,
  ) -> Result<ReconstructedSequences, Report> {
    reconstruct_preorder(graph, tips.include_leaves, |node| {
      let seq_info = &node_states[&node.key];
      if node.is_leaf {
        return Ok(self.reconstruct_leaf_sequence(seq_info, tips.impute));
      }
      let mut resolve = if sample_mode.samples_node(node.is_root) {
        Resolve::Sample(&mut *rng)
      } else {
        Resolve::Argmax
      };
      Ok(assign_sequence_sampled(seq_info, &self.alphabet, &mut resolve))
    })
  }

  fn reconstruct_leaf_sequence(&self, seq_info: &DenseNodeState, impute: bool) -> Seq {
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
  }

  fn indexed_kind(&self) -> IndexedKind<'_> {
    IndexedKind::Dense {
      alphabet: &self.alphabet,
      length: self.length,
      leaves: &self.obs_leaves,
    }
  }
}

impl MarginalPasses for PartitionMarginalDense {
  type Node = DenseNodeState;
  type BackwardInput = ();
  type Backward = DenseEdgeBackward;
  type Forward = DenseEdgeForward;
  type Estimate = DenseEdgeEstimate;

  fn backward_input(_node_states: &BTreeMap<GraphNodeKey, DenseNodeState>) -> &Self::BackwardInput {
    &()
  }

  fn backward_input_with_reset_log_lh(_node_states: &BTreeMap<GraphNodeKey, DenseNodeState>) -> Self::BackwardInput {}

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

pub(crate) fn assign_sequence(seq_info: &DenseNodeState, alphabet: &Alphabet) -> Seq {
  assign_sequence_sampled(seq_info, alphabet, &mut Resolve::Argmax)
}

fn assign_sequence_sampled(seq_info: &DenseNodeState, alphabet: &Alphabet, resolve: &mut Resolve) -> Seq {
  let mut seq = prof2seq_sampled(&seq_info.profile, alphabet, resolve);
  for gap in &seq_info.seq.gaps {
    seq[gap.0..gap.1].fill(alphabet.gap());
  }
  for unk in &seq_info.seq.unknown {
    seq[unk.0..unk.1].fill(alphabet.unknown());
  }
  seq
}

fn prof2seq_sampled(profile: &DenseSeqDistribution, alphabet: &Alphabet, resolve: &mut Resolve) -> Seq {
  let mut seq = seq! {};
  for row in profile.dis.rows() {
    seq.push(alphabet.char(resolve_profile(row, resolve)));
  }
  seq
}

pub type DenseMarginalEdges = MarginalEdges<DenseEdgeBackward, DenseEdgeForward, DenseEdgeEstimate>;
