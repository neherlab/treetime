//! Turning a node's posterior into a sequence.
//!
//! `SparseSeqInfo::sequence` always holds the *parsimony* sequence: the root sequence with each
//! edge's Fitch substitutions and indels applied, masked by the node's own missing data. It is built
//! once by the marginal forward pass and never rewritten afterwards.
//!
//! The MAP sequence is a view over that: the parsimony sequence with every position in
//! `profile.variable` resolved from its posterior. Positions absent from `profile.variable` resolved
//! to a single state during message combination, and that state is by construction the parsimony
//! state, so taking them from the parsimony chain is exact.
//!
//! Keeping the stored chain free of MAP states is what makes the two agree. Chaining a node's MAP
//! sequence off its parent's would silently propagate any node whose argmax disagrees with Fitch
//! into every descendant that resolved the position.

use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::sample::resolve_profile;
use crate::partition::storage::sparse::{SparseEdgePartition, SparseNodePartition};
use treetime_primitives::{AsciiChar, Seq};
use treetime_utils::array::ndarray::argmax_first;

/// Extend the parsimony chain from the parent across one edge.
pub(crate) fn parsimony_seq(
  parent_seq: &Seq,
  edge: &SparseEdgePartition,
  node: &SparseNodePartition,
  alphabet: &Alphabet,
) -> Seq {
  let mut seq = parent_seq.clone();

  for sub in edge.fitch_subs() {
    seq[sub.pos()] = sub.qry();
  }

  for indel in &edge.indels {
    if indel.is_deletion() {
      seq[indel.range.0..indel.range.1].fill(alphabet.gap());
    } else {
      seq[indel.range.0..indel.range.1].copy_from_slice(&indel.seq);
    }
  }

  for r in &node.seq.unknown {
    seq[r.0..r.1].fill(alphabet.unknown());
  }

  seq
}

/// The node's most likely sequence.
pub(crate) fn map_seq(node: &SparseNodePartition, alphabet: &Alphabet) -> Seq {
  map_seq_sampled(node, alphabet, false, &mut rand::thread_rng())
}

/// The node's most likely sequence, or one draw from its posterior when `sample` is set.
///
/// Gaps stay gaps: a gap marks missing data, not an uncertain base, so a variable position on a
/// deletion is left alone. Unknown (`N`) positions are resolved from the posterior, which is the
/// inference that fills them in.
pub(crate) fn map_seq_sampled(
  node: &SparseNodePartition,
  alphabet: &Alphabet,
  sample: bool,
  rng: &mut dyn rand::RngCore,
) -> Seq {
  let mut seq = node.seq.sequence.clone();

  for (&pos, var) in &node.profile.variable {
    if seq[pos] != alphabet.gap() {
      seq[pos] = alphabet.char(resolve_profile(var.dis.view(), sample, rng));
    }
  }

  // Sampling draws every position, not only the variable ones: an invariant position is still a
  // distribution, just one shared across all positions holding that character. `fixed` is keyed by
  // canonical state, so gaps and unknowns find no entry and stay as they are.
  if sample {
    for pos in 0..seq.len() {
      if node.profile.variable.contains_key(&pos) {
        continue;
      }
      if let Some(fixed) = node.profile.fixed.get(&seq[pos]) {
        seq[pos] = alphabet.char(resolve_profile(fixed.view(), true, rng));
      }
    }
  }

  seq
}

/// The node's most likely state at one position.
pub(crate) fn map_state(node: &SparseNodePartition, pos: usize, alphabet: &Alphabet) -> AsciiChar {
  match node.profile.variable.get(&pos) {
    Some(var) => alphabet.char(argmax_first(&var.dis.view()).unwrap_or(0)),
    None => node.seq.sequence.get(pos).copied().unwrap_or_else(|| alphabet.char(0)),
  }
}

/// Reconstruct a leaf output sequence from its own observed data.
///
/// Mirrors the dense backend, which keeps the observed tip sequence and never chains a leaf through
/// its parent. Without imputation the tip emits exactly its observed input. With imputation every
/// ambiguous or unknown position (`N` and IUPAC codes, but not gaps) is resolved to the argmax of the
/// leaf marginal posterior. That posterior is the parent marginal evolved across the branch
/// (`msg_from_parent`) restricted by the observed ambiguity mask, matching v0's per-leaf marginal
/// profile; on long tip branches this differs from simply copying the parent MAP state.
pub(crate) fn reconstruct_leaf_sequence(
  node: &SparseNodePartition,
  edge: Option<&SparseEdgePartition>,
  parent: Option<&SparseNodePartition>,
  impute: bool,
  alphabet: &Alphabet,
) -> Seq {
  let mut seq = node.seq.sequence.clone();

  // Fitch compression stores each observed IUPAC ambiguity as a single resolved canonical state in
  // `seq.sequence`; restore the observed ambiguity code so a non-imputing tip echoes its true input.
  // Unknown (`N`) positions already hold the unknown character in `seq.sequence`.
  for (&pos, states) in &node.seq.fitch.variable {
    seq[pos] = alphabet.set_to_char(*states);
  }

  // A tip with no parent edge (a single-node tree) has no down-message to impute from.
  let (true, Some(edge), Some(parent)) = (impute, edge, parent) else {
    return seq;
  };
  let down = &edge.msg_from_parent;

  // Imputable positions are the unknown (`N`) ranges and the IUPAC ambiguity positions detected from
  // the partition structure, not from the character being non-canonical (Fitch may already have
  // resolved an IUPAC code to a canonical state). Gaps are inferred deletions and are left untouched.
  let unknown_positions = node.seq.unknown.iter().flat_map(|&(start, end)| start..end);
  let ambiguous_positions = node.seq.fitch.variable.keys().copied();
  for pos in unknown_positions.chain(ambiguous_positions) {
    let observed = seq[pos];
    // Parent posterior at this site: an explicit variable distribution when the parent varies here,
    // otherwise the fixed-column distribution keyed by the parent's MAP state.
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
