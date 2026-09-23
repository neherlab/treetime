use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::sample::{Resolve, resolve_profile};
use crate::partition::storage::sparse::{SparseEdgeObs, SparseNodeObs, SparseNodeState, SparseSeqDistribution};
use treetime_primitives::{AsciiChar, Seq};
use treetime_utils::array::ndarray::argmax_first;

pub fn parsimony_seq(parent_seq: &Seq, edge_obs: &SparseEdgeObs, node_obs: &SparseNodeObs, alphabet: &Alphabet) -> Seq {
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

pub fn map_seq(node: &SparseNodeState, alphabet: &Alphabet) -> Seq {
  map_seq_sampled(node, alphabet, &mut Resolve::Argmax)
}

pub fn map_seq_sampled(node: &SparseNodeState, alphabet: &Alphabet, resolve: &mut Resolve) -> Seq {
  let mut seq = node.sequence.clone();

  for (&pos, var) in &node.profile.variable {
    if seq[pos] != alphabet.gap() {
      seq[pos] = alphabet.char(resolve_profile(var.dis.view(), resolve));
    }
  }

  if matches!(resolve, Resolve::Sample(_)) {
    for pos in 0..seq.len() {
      if node.profile.variable.contains_key(&pos) {
        continue;
      }
      if let Some(fixed) = node.profile.fixed.get(&seq[pos]) {
        seq[pos] = alphabet.char(resolve_profile(fixed.view(), resolve));
      }
    }
  }

  seq
}

pub fn reconstruct_leaf_sequence(
  node: &SparseNodeState,
  node_obs: &SparseNodeObs,
  msg_from_parent: Option<&SparseSeqDistribution>,
  parent: Option<&SparseNodeState>,
  impute: bool,
  alphabet: &Alphabet,
) -> Seq {
  let mut seq = node.sequence.clone();

  for (&pos, states) in &node_obs.fitch.variable {
    seq[pos] = alphabet.set_to_char(*states);
  }

  let (true, Some(down), Some(parent)) = (impute, msg_from_parent, parent) else {
    return seq;
  };

  let unknown_positions = node_obs.unknown.iter().flat_map(|&(start, end)| start..end);
  let ambiguous_positions = node_obs.fitch.variable.keys().copied();
  for pos in unknown_positions.chain(ambiguous_positions) {
    let observed = seq[pos];
    let posterior = down
      .variable
      .get(&pos)
      .map(|var| &var.dis)
      .or_else(|| down.fixed.get(&map_state(parent, pos, alphabet)));
    let (Some(posterior), Ok(mask)) = (posterior, alphabet.get_profile(observed)) else {
      continue;
    };
    let combined = posterior * mask;
    if let Some(idx) = argmax_first(&combined.view()) {
      seq[pos] = alphabet.char(idx);
    }
  }

  seq
}

pub fn map_state(node: &SparseNodeState, pos: usize, alphabet: &Alphabet) -> AsciiChar {
  match node.profile.variable.get(&pos) {
    Some(var) => alphabet.char(argmax_first(&var.dis.view()).unwrap_or(0)),
    None => node.sequence.get(pos).copied().unwrap_or_else(|| alphabet.char(0)),
  }
}
