use crate::alphabet::alphabet::Alphabet;
use crate::partition::marginal::posterior_states::{AmbiguityThreshold, most_likely_state};
use crate::partition::marginal::sample::{Resolve, resolve_profile};
use crate::partition::storage::sparse::{SparseEdgeObs, SparseNodeObs, SparseNodeState, SparseSeqDistribution};
use crate::seq::overlay::SeqOverlay;
use crate::seq::sink::AmbiguousSites;
use ndarray::{Array1, aview1};
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_primitives::{AsciiChar, Seq};

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

  for r in &node_obs.unknown {
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

  for (pos, state) in variable_states(node, alphabet, resolve) {
    seq[pos] = state;
  }

  if matches!(resolve, Resolve::Sample(_)) {
    for pos in 0..seq.len() {
      if node.profile.variable.contains_key(pos) {
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
    .profile
    .variable
    .iter()
    .filter(|&(&pos, _)| node.sequence[pos] != alphabet.gap())
    .map(move |(&pos, var)| (pos, alphabet.char(resolve_profile(aview1(var.dis), resolve))))
}

pub(crate) fn reconstruct_leaf_sequence(
  node: &SparseNodeState,
  node_obs: &SparseNodeObs,
  msg_from_parent: Option<&SparseSeqDistribution>,
  parent: Option<&SparseNodeState>,
  impute: bool,
  alphabet: &Alphabet,
) -> Seq {
  let mut seq = (*node.sequence).clone();

  for (&pos, states) in &node_obs.fitch.variable {
    seq[pos] = alphabet.set_to_char(*states);
  }

  let (true, Some(down), Some(parent)) = (impute, msg_from_parent, parent) else {
    return seq;
  };

  let unknown_positions = node_obs.unknown.iter().flat_map(|&(start, end)| start..end);
  let ambiguous_positions = node_obs.fitch.variable.keys().copied();
  for pos in unknown_positions.chain(ambiguous_positions) {
    seq[pos] = impute_state(seq[pos], pos, down, parent, alphabet);
  }

  seq
}

pub(crate) fn internal_ambiguity(
  node: &SparseNodeState,
  alphabet: &Alphabet,
  threshold: AmbiguityThreshold,
) -> AmbiguousSites {
  let mut sites: AmbiguousSites = node
    .profile
    .variable
    .iter()
    .filter(|&(&pos, _)| node.sequence[pos] != alphabet.gap())
    .filter_map(|(&pos, var)| {
      threshold
        .ambiguity_char(&aview1(var.dis), alphabet)
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
      .filter(|&(pos, _)| !node.profile.variable.contains_key(pos))
      .filter_map(|(pos, state)| fixed_codes.get(state).map(|&code| (pos, code)));
    sites.extend(fixed_sites);
  }
  sites
}

pub(crate) fn leaf_ambiguity(
  node: &SparseNodeState,
  node_obs: &SparseNodeObs,
  down: &SparseSeqDistribution,
  parent: &SparseNodeState,
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
      let posterior = leaf_posterior(observed, pos, down, parent, alphabet)?;
      threshold.ambiguity_char(&posterior, alphabet).map(|code| (pos, code))
    })
    .collect()
}

pub(crate) fn impute_state(
  observed: AsciiChar,
  pos: usize,
  down: &SparseSeqDistribution,
  parent: &SparseNodeState,
  alphabet: &Alphabet,
) -> AsciiChar {
  leaf_posterior(observed, pos, down, parent, alphabet)
    .and_then(|posterior| most_likely_state(&posterior))
    .map_or(observed, |idx| alphabet.char(idx))
}

fn leaf_posterior(
  observed: AsciiChar,
  pos: usize,
  down: &SparseSeqDistribution,
  parent: &SparseNodeState,
  alphabet: &Alphabet,
) -> Option<Array1<f64>> {
  let from_parent = down.variable.get(pos).map(|var| var.dis).or_else(|| {
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

pub(crate) fn map_state(node: &SparseNodeState, pos: usize, alphabet: &Alphabet) -> AsciiChar {
  match node.profile.variable.get(pos) {
    Some(var) => alphabet.char(most_likely_state(&aview1(var.dis)).unwrap_or(0)),
    None => node.sequence.get(pos).copied().unwrap_or_else(|| alphabet.char(0)),
  }
}
