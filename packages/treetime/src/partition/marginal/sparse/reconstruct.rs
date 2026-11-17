use crate::alphabet::alphabet::Alphabet;
use crate::partition::marginal::mixture::site_log_lh::{INLINE_CATEGORIES, SparseMixture};
use crate::partition::marginal::mixture::sparse::SparseCategory;
use crate::partition::marginal::posterior_states::{AmbiguityThreshold, most_likely_state};
use crate::partition::marginal::sample::{Resolve, resolve_profile};
use crate::partition::marginal::sparse::globally_missing::SitePosterior;
use crate::partition::storage::sparse::{
  SparseEdgeForward, SparseEdgeObs, SparseNodeObs, SparseNodeState, SparseSeqDistribution,
};
use crate::seq::overlay::SeqOverlay;
use crate::seq::sink::AmbiguousSites;
use itertools::izip;
use ndarray::{Array1, aview1};
use smallvec::{SmallVec, smallvec};
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{AsciiChar, Seq};
use treetime_utils::interval::range_difference::range_difference;

pub(crate) fn parsimony_seq(
  parent_seq: &Seq,
  edge_obs: &SparseEdgeObs,
  node_obs: &SparseNodeObs,
  alphabet: &Alphabet,
) -> Seq {
  let mut seq = parent_seq.clone();

  for sub in edge_obs.fitch_subs() {
    seq[sub.pos()] = sub.qry();
  }

  for indel in &edge_obs.indels {
    if indel.is_deletion() {
      seq[indel.range.0..indel.range.1].fill(alphabet.gap());
    } else {
      seq[indel.range.0..indel.range.1].copy_from_slice(&indel.seq);
    }
  }

  for r in range_difference(&node_obs.unknown, &node_obs.gaps) {
    seq[r.0..r.1].fill(alphabet.unknown());
  }

  seq
}

pub(crate) fn map_seq(node: &SparseNodeState, alphabet: &Alphabet) -> Seq {
  map_seq_sampled(node, alphabet, &mut Resolve::Argmax)
}

pub(crate) fn map_seq_overlay(node: &SparseNodeState, alphabet: &Alphabet) -> SeqOverlay {
  let overrides = variable_states(node, alphabet, &mut Resolve::Argmax).collect();
  SeqOverlay::new(Arc::clone(&node.sequence), overrides)
}

pub(crate) fn map_seq_sampled(node: &SparseNodeState, alphabet: &Alphabet, resolve: &mut Resolve) -> Seq {
  let mut seq = (*node.sequence).clone();

  if matches!(resolve, Resolve::Argmax) {
    node.globally_missing.write_most_likely(&mut seq, alphabet.gap());
    for (pos, state) in explicit_states(node, alphabet, resolve) {
      seq[pos] = state;
    }
  } else {
    for (pos, state) in variable_states(node, alphabet, resolve) {
      seq[pos] = state;
    }
  }

  if matches!(resolve, Resolve::Sample(_)) {
    for pos in 0..seq.len() {
      if node.has_posterior_at(pos) {
        continue;
      }
      if let Some(fixed) = node.profile.fixed.get(seq[pos]) {
        seq[pos] = alphabet.char(resolve_profile(aview1(&fixed.dis), resolve));
      }
    }
  }

  seq
}

fn variable_states<'a>(
  node: &'a SparseNodeState,
  alphabet: &'a Alphabet,
  resolve: &'a mut Resolve<'_>,
) -> impl Iterator<Item = (usize, AsciiChar)> + 'a {
  node
    .posteriors()
    .filter(|&(pos, _)| !is_gap(node, pos, alphabet))
    .map(move |(pos, posterior)| {
      let state = match posterior {
        SitePosterior::Explicit(var) => alphabet.char(resolve_profile(aview1(var.dis), resolve)),
        SitePosterior::Equilibrium => node.globally_missing.resolve(alphabet, resolve),
      };
      (pos, state)
    })
}

fn explicit_states<'a>(
  node: &'a SparseNodeState,
  alphabet: &'a Alphabet,
  resolve: &'a mut Resolve<'_>,
) -> impl Iterator<Item = (usize, AsciiChar)> + 'a {
  node
    .profile
    .variable
    .iter()
    .filter(|&(&pos, _)| !is_gap(node, pos, alphabet))
    .map(move |(&pos, var)| (pos, alphabet.char(resolve_profile(aview1(var.dis), resolve))))
}

pub(crate) fn reconstruct_leaf_sequence(
  node: &SparseNodeState,
  node_obs: &SparseNodeObs,
  down: Option<&LeafDown<'_>>,
  impute: bool,
  alphabet: &Alphabet,
) -> Seq {
  let mut seq = (*node.sequence).clone();

  for (&pos, states) in &node_obs.fitch.variable {
    seq[pos] = alphabet.set_to_char(*states);
  }

  let (true, Some(down)) = (impute, down) else {
    return seq;
  };

  let unknown_positions = node_obs.unknown.iter().flat_map(|&(start, end)| start..end);
  let ambiguous_positions = node_obs.fitch.variable.keys().copied();
  for pos in unknown_positions.chain(ambiguous_positions) {
    seq[pos] = impute_state(seq[pos], pos, down, alphabet);
  }

  seq
}

pub(crate) fn internal_ambiguity(
  node: &SparseNodeState,
  alphabet: &Alphabet,
  threshold: AmbiguityThreshold,
) -> AmbiguousSites {
  let equilibrium_code = threshold.ambiguity_char(&aview1(node.globally_missing.profile()), alphabet);
  let mut sites: AmbiguousSites = node
    .posteriors()
    .filter(|&(pos, _)| !is_gap(node, pos, alphabet))
    .filter_map(|(pos, posterior)| {
      match posterior {
        SitePosterior::Explicit(var) => threshold.ambiguity_char(&aview1(var.dis), alphabet),
        SitePosterior::Equilibrium => equilibrium_code,
      }
      .map(|code| (pos, code))
    })
    .collect();
  let fixed_codes: BTreeMap<AsciiChar, AsciiChar> = node
    .profile
    .fixed
    .iter()
    .filter_map(|fixed| {
      threshold
        .ambiguity_char(&aview1(&fixed.dis), alphabet)
        .map(|code| (fixed.state, code))
    })
    .collect();
  if !fixed_codes.is_empty() {
    let fixed_sites = node
      .sequence
      .iter()
      .enumerate()
      .filter(|&(pos, _)| !node.has_posterior_at(pos))
      .filter_map(|(pos, state)| fixed_codes.get(state).map(|&code| (pos, code)));
    sites.extend(fixed_sites);
  }
  sites
}

pub(crate) fn leaf_ambiguity(
  node: &SparseNodeState,
  node_obs: &SparseNodeObs,
  down: &LeafDown<'_>,
  alphabet: &Alphabet,
  threshold: AmbiguityThreshold,
) -> AmbiguousSites {
  let unknown_positions = node_obs.unknown.iter().flat_map(|&(start, end)| start..end);
  let ambiguous_positions = node_obs.fitch.variable.keys().copied();
  unknown_positions
    .chain(ambiguous_positions)
    .filter_map(|pos| {
      let observed = node_obs
        .fitch
        .variable
        .get(&pos)
        .map_or(node.sequence[pos], |&states| alphabet.set_to_char(states));
      let posterior = down.posterior(observed, pos, alphabet)?;
      threshold.ambiguity_char(&posterior, alphabet).map(|code| (pos, code))
    })
    .collect()
}

pub(crate) fn impute_state(observed: AsciiChar, pos: usize, down: &LeafDown<'_>, alphabet: &Alphabet) -> AsciiChar {
  down
    .posterior(observed, pos, alphabet)
    .and_then(|posterior| most_likely_state(&posterior))
    .map_or(observed, |idx| alphabet.char(idx))
}

pub(crate) struct LeafDown<'a> {
  categories: SmallVec<[(&'a SparseSeqDistribution, &'a SparseNodeState); INLINE_CATEGORIES]>,
  mixture: Option<&'a SparseMixture>,
}

impl<'a> LeafDown<'a> {
  pub(crate) fn single(down: &'a SparseSeqDistribution, parent: &'a SparseNodeState) -> Self {
    Self {
      categories: smallvec![(down, parent)],
      mixture: None,
    }
  }

  pub(crate) fn mixture(
    categories: SmallVec<[(&'a SparseSeqDistribution, &'a SparseNodeState); INLINE_CATEGORIES]>,
    mixture: &'a SparseMixture,
  ) -> Self {
    Self {
      categories,
      mixture: Some(mixture),
    }
  }

  pub(crate) fn posterior(&self, observed: AsciiChar, pos: usize, alphabet: &Alphabet) -> Option<Array1<f64>> {
    let Some(mixture) = self.mixture else {
      let &(down, parent) = self.categories.first()?;
      return leaf_posterior(observed, pos, down, parent, alphabet);
    };
    let (_, parent) = self.categories.first()?;
    let weights = mixture.posterior(pos, parent.sequence[pos]);
    izip!(&self.categories, weights)
      .map(|(&(down, parent), &weight)| leaf_posterior(observed, pos, down, parent, alphabet).map(|post| post * weight))
      .reduce(|total, term| match (total, term) {
        (Some(total), Some(term)) => Some(total + term),
        (total, term) => total.or(term),
      })
      .flatten()
  }
}

pub(crate) fn leaf_posterior(
  observed: AsciiChar,
  pos: usize,
  down: &SparseSeqDistribution,
  parent: &SparseNodeState,
  alphabet: &Alphabet,
) -> Option<Array1<f64>> {
  let from_parent = down
    .variable
    .get(pos)
    .or_else(|| parent.globally_missing.get(pos))
    .map(|var| var.dis)
    .or_else(|| {
      down
        .fixed
        .get(map_state(parent, pos, alphabet))
        .map(|fixed| fixed.dis.as_slice())
    })?;
  let Ok(mask) = alphabet.get_profile(observed) else {
    return None;
  };
  let mut posterior = &aview1(from_parent) * mask;
  let total = posterior.sum();
  if total > 0.0 && total.is_finite() {
    posterior /= total;
  }
  Some(posterior)
}

fn is_gap(node: &SparseNodeState, pos: usize, alphabet: &Alphabet) -> bool {
  node.sequence[pos] == alphabet.gap()
}

pub(crate) fn map_state(node: &SparseNodeState, pos: usize, alphabet: &Alphabet) -> AsciiChar {
  match node.profile.variable.get(pos) {
    Some(var) => alphabet.char(most_likely_state(&aview1(var.dis)).unwrap_or(0)),
    None => node
      .globally_missing
      .most_likely_state_at(pos)
      .or_else(|| node.sequence.get(pos).copied())
      .unwrap_or_else(|| alphabet.char(0)),
  }
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum LeafDowns<'a> {
  Single {
    forward: &'a EdgeMap<SparseEdgeForward>,
    node_states: &'a NodeMap<SparseNodeState>,
  },
  Mixture {
    categories: &'a [SparseCategory],
    mixture: &'a SparseMixture,
  },
}

impl<'a> LeafDowns<'a> {
  pub(crate) fn down(self, edge_key: GraphEdgeKey, parent_key: GraphNodeKey) -> LeafDown<'a> {
    match self {
      Self::Single { forward, node_states } => {
        LeafDown::single(&forward[&edge_key].msg_from_parent, &node_states[&parent_key])
      },
      Self::Mixture { categories, mixture } => LeafDown::mixture(
        categories
          .iter()
          .map(|category| {
            (
              &category.edges.forward[&edge_key].msg_from_parent,
              &category.node_states[&parent_key],
            )
          })
          .collect(),
        mixture,
      ),
    }
  }
}
