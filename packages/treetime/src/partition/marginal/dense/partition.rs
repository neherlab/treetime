use crate::alphabet::alphabet::Alphabet;
use crate::constants::MIN_BRANCH_LENGTH_FRACTION;
use crate::gtr::gtr::GTR;
use crate::partition::marginal::posterior_states::{AmbiguityThreshold, most_likely_state};
use crate::partition::marginal::sample::{Resolve, SampleMode, resolve_profile};
use crate::partition::marginal::sequences::sample_internal_sequences;
use crate::partition::marginal::shared::data::DenseInputs;
use crate::partition::marginal::shared::pass::{IndexedKind, IndexedPartition};
use crate::partition::marginal::shared::update::MarginalEdges;
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::partition::storage::dense::{
  DenseEdgeBackward, DenseEdgeEstimate, DenseEdgeForward, DenseLeafObs, DenseNodeState, DenseSeqDistribution,
};
use crate::seq::alignment::{NodeSeqInput, get_common_length_of_node_inputs};
use crate::seq::indel::InDel;
use crate::seq::mutation::Sub;
use crate::seq::overlay::SeqOverlay;
use crate::seq::sink::AmbiguousSites;
use crate::{make_internal_report, make_report};
use eyre::Report;
use itertools::{Itertools, izip};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{Seq, seq};
use treetime_utils::cast::ToF64Exact;
use treetime_utils::interval::range::range_contains;
use treetime_utils::interval::range_union::range_union;

#[derive(Clone, Debug)]
pub struct PartitionMarginalDense {
  pub(crate) inputs: DenseInputs,
  pub(crate) alphabet: Alphabet,
  pub(crate) length: usize,
  pub(crate) obs_leaves: NodeMap<DenseLeafObs>,
}

impl PartitionMarginalDense {
  pub(crate) fn new(alphabet: Alphabet, graph: &Graph, node_inputs: &NodeMap<NodeSeqInput>) -> Result<Self, Report> {
    let length = get_common_length_of_node_inputs(node_inputs)?;
    let obs_leaves = graph
      .get_leaves()
      .map(|leaf| {
        let leaf_key = leaf.key();
        let node = &node_inputs[&leaf_key];
        let seq = node
          .seq
          .as_ref()
          .ok_or_else(|| make_report!("Leaf sequence not found: '{}'", node.name))?;
        Ok((leaf_key, DenseLeafObs::new(seq, &alphabet)))
      })
      .collect::<Result<NodeMap<_>, Report>>()?;
    let min_branch_length = MIN_BRANCH_LENGTH_FRACTION / length.to_f64_exact();
    Ok(Self {
      inputs: DenseInputs {
        min_branch_length,
        filter_uninformative_root: true,
      },
      alphabet,
      length,
      obs_leaves,
    })
  }

  #[must_use]
  pub(crate) fn with_total_length(self, total_length: usize) -> Self {
    Self {
      inputs: DenseInputs {
        min_branch_length: MIN_BRANCH_LENGTH_FRACTION / total_length.to_f64_exact(),
        ..self.inputs
      },
      ..self
    }
  }

  pub(crate) fn edge_subs(
    &self,
    node_states: &NodeMap<DenseNodeState>,
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
      let parent_state = self.alphabet.char(most_likely_state(&parent).unwrap_or(0));
      let child_state = self.alphabet.char(most_likely_state(&child).unwrap_or(0));
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

  pub(crate) fn edge_indels(&self, estimates: &EdgeMap<DenseEdgeEstimate>, edge_key: GraphEdgeKey) -> Vec<InDel> {
    estimates[&edge_key].indels.clone()
  }

  pub(crate) fn root_sequence(&self, node_states: &NodeMap<DenseNodeState>, graph: &Graph) -> Result<Seq, Report> {
    Ok(assign_sequence(&node_states[&graph.root_key()?], &self.alphabet))
  }

  pub(crate) fn edge_effective_length(
    &self,
    node_states: &NodeMap<DenseNodeState>,
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
    backward: &EdgeMap<DenseEdgeBackward>,
    forward: &EdgeMap<DenseEdgeForward>,
    edge_key: GraphEdgeKey,
  ) -> Result<OptimizationContribution, Report> {
    OptimizationContribution::from_dense(gtr, &backward[&edge_key], &forward[&edge_key])
  }

  pub(crate) fn edge_indel_count(&self, estimates: &EdgeMap<DenseEdgeEstimate>, edge_key: GraphEdgeKey) -> usize {
    estimates[&edge_key].indels.len()
  }

  pub(crate) fn extract_ancestral_sequence(
    &self,
    node_states: &NodeMap<DenseNodeState>,
    node_key: GraphNodeKey,
  ) -> SeqOverlay {
    SeqOverlay::from(
      node_states
        .get(node_key)
        .map_or_else(|| seq! {}, |seq_info| assign_sequence(seq_info, &self.alphabet)),
    )
  }

  pub(crate) fn node_sequence(
    &self,
    graph: &Graph,
    node_states: &NodeMap<DenseNodeState>,
    impute: bool,
    node_key: GraphNodeKey,
  ) -> Result<Seq, Report> {
    let seq_info = &node_states[&node_key];
    let is_leaf = graph
      .get_node(node_key)
      .ok_or_else(|| make_internal_report!("Node {node_key} not found while reconstructing its sequence"))?
      .is_leaf();
    Ok(if is_leaf {
      self.reconstruct_leaf_sequence(seq_info, impute)
    } else {
      assign_sequence(seq_info, &self.alphabet)
    })
  }

  pub(crate) fn sample_sequences(
    &self,
    graph: &Graph,
    node_states: &NodeMap<DenseNodeState>,
    sample_mode: SampleMode,
    rng: &mut dyn rand::Rng,
  ) -> Result<NodeMap<Seq>, Report> {
    sample_internal_sequences(graph, sample_mode, |key| {
      assign_sequence_sampled(&node_states[&key], &self.alphabet, &mut Resolve::Sample(&mut *rng))
    })
  }

  pub(crate) fn node_ambiguity(
    &self,
    graph: &Graph,
    node_states: &NodeMap<DenseNodeState>,
    impute: bool,
    node_key: GraphNodeKey,
    threshold: AmbiguityThreshold,
  ) -> Result<AmbiguousSites, Report> {
    let seq_info = &node_states[&node_key];
    let is_leaf = graph
      .get_node(node_key)
      .ok_or_else(|| make_internal_report!("Node {node_key} not found while reporting its ambiguous sites"))?
      .is_leaf();
    let profile = &seq_info.profile.dis;
    let sites = if is_leaf {
      self
        .imputed_positions(seq_info, impute)
        .filter_map(|pos| {
          threshold
            .ambiguity_char(&profile.row(pos), &self.alphabet)
            .map(|code| (pos, code))
        })
        .collect()
    } else {
      (0..profile.nrows())
        .filter(|&pos| !range_contains(&seq_info.seq.gaps, pos) && !range_contains(&seq_info.seq.unknown, pos))
        .filter_map(|pos| {
          threshold
            .ambiguity_char(&profile.row(pos), &self.alphabet)
            .map(|code| (pos, code))
        })
        .collect()
    };
    Ok(sites)
  }

  fn imputed_positions<'a>(&'a self, seq_info: &'a DenseNodeState, impute: bool) -> impl Iterator<Item = usize> + 'a {
    let seq = &seq_info.seq.sequence;
    let imputes = impute && seq_info.profile.dis.nrows() == seq.len();
    seq
      .iter()
      .positions(move |&ch| imputes && !self.alphabet.is_canonical(ch) && !self.alphabet.is_gap(ch))
  }

  fn reconstruct_leaf_sequence(&self, seq_info: &DenseNodeState, impute: bool) -> Seq {
    let mut seq = seq_info.seq.sequence.clone();
    for pos in self.imputed_positions(seq_info, impute) {
      if let Some(idx) = most_likely_state(&seq_info.profile.dis.row(pos)) {
        seq[pos] = self.alphabet.char(idx);
      }
    }
    seq
  }
}

impl IndexedPartition for PartitionMarginalDense {
  fn inputs(&self) -> &DenseInputs {
    &self.inputs
  }

  fn indexed_kind(&self) -> IndexedKind<'_> {
    IndexedKind::Dense {
      alphabet: &self.alphabet,
      length: self.length,
      leaves: &self.obs_leaves,
    }
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
