use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::sample::resolve_profile;
use crate::partition::storage::sparse::{SparseEdgePartition, SparseNodePartition};
use treetime_primitives::Seq;
use treetime_utils::array::ndarray::argmax_first;

pub(crate) fn reconstruct_map_seq(
  base_seq: &Seq,
  edge: Option<&SparseEdgePartition>,
  node: &SparseNodePartition,
  alphabet: &Alphabet,
) -> Seq {
  let mut rng = rand::thread_rng();
  reconstruct_map_seq_sampled(base_seq, edge, node, alphabet, false, &mut rng)
}

pub(crate) fn reconstruct_map_seq_sampled(
  base_seq: &Seq,
  edge: Option<&SparseEdgePartition>,
  node: &SparseNodePartition,
  alphabet: &Alphabet,
  sample: bool,
  rng: &mut dyn rand::RngCore,
) -> Seq {
  let mut seq = if let Some(edge) = edge {
    let mut seq = base_seq.clone();
    for m in edge.fitch_subs() {
      seq[m.pos()] = m.qry();
    }
    for indel in &edge.indels {
      if indel.is_deletion() {
        seq[indel.range.0..indel.range.1].fill(alphabet.gap());
      } else {
        seq[indel.range.0..indel.range.1].copy_from_slice(&indel.seq);
      }
    }
    seq
  } else {
    base_seq.clone()
  };

  for r in &node.seq.unknown {
    seq[r.0..r.1].fill(alphabet.unknown());
  }

  for (pos, states) in &node.profile.variable {
    seq[*pos] = alphabet.char(resolve_profile(states.dis.view(), sample, rng));
  }

  if sample {
    let variable_positions = &node.profile.variable;
    let fixes: Vec<(usize, usize)> = seq
      .iter()
      .enumerate()
      .filter(|(pos, _)| !variable_positions.contains_key(pos))
      .filter_map(|(pos, &ch)| {
        node
          .profile
          .fixed
          .get(&ch)
          .map(|fp| (pos, resolve_profile(fp.view(), true, rng)))
      })
      .collect();
    for (pos, idx) in fixes {
      seq[pos] = alphabet.char(idx);
    }
  }

  // Re-apply deletions last. The profile writes above are keyed by position and would otherwise
  // resurrect a residue at a site the edge reports as deleted, leaving the reconstructed sequence
  // disagreeing with its own indel track.
  if let Some(edge) = edge {
    for indel in &edge.indels {
      if indel.is_deletion() {
        seq[indel.range.0..indel.range.1].fill(alphabet.gap());
      }
    }
  }

  seq
}

/// Reconstruct a leaf output sequence from its own observed data.
///
/// Mirrors the dense backend, which keeps the observed tip sequence and never chains a leaf through
/// its parent. The forward pass leaves a leaf's `seq.sequence` equal to its observed input, so
/// without imputation the tip emits exactly that. With imputation every ambiguous or unknown
/// position (`N` and IUPAC codes, but not gaps) is resolved to the argmax of the leaf marginal
/// posterior. That posterior is the parent marginal evolved across the branch (`msg_from_parent`)
/// restricted by the observed ambiguity mask, matching v0's per-leaf marginal profile; on long tip
/// branches this differs from simply copying the parent MAP state.
pub(crate) fn reconstruct_leaf_sequence(
  node: &SparseNodePartition,
  edge: Option<&SparseEdgePartition>,
  base_seq: &Seq,
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
  let (true, Some(edge)) = (impute, edge) else {
    return seq;
  };
  let down = &edge.msg_from_parent;

  // Impute every ambiguous or unknown position to the leaf marginal argmax. Imputable positions are
  // the unknown (`N`) ranges and the IUPAC ambiguity positions detected from the partition structure,
  // not from the character being non-canonical (Fitch may already have resolved an IUPAC code to a
  // canonical state). Gaps are inferred deletions and are left untouched.
  let unknown_positions = node.seq.unknown.iter().flat_map(|&(start, end)| start..end);
  let ambiguous_positions = node.seq.fitch.variable.keys().copied();
  for pos in unknown_positions.chain(ambiguous_positions) {
    let observed = seq[pos];
    // Parent posterior at this site: an explicit variable distribution when the parent varies here,
    // otherwise the fixed-column distribution keyed by the parent MAP state.
    let posterior = down
      .variable
      .get(&pos)
      .map(|var| &var.dis)
      .or_else(|| base_seq.get(pos).and_then(|parent_char| down.fixed.get(parent_char)));
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
