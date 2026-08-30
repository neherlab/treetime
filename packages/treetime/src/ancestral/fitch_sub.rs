use crate::alphabet::alphabet::{Alphabet, FILL_CHAR, NON_CHAR, VARIABLE_CHAR};
use crate::partition::storage::sparse::{SparseEdgePartition, SparseSeqInfo};
use crate::seq::composition::Composition;
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use treetime_primitives::{AlphabetLike, AsciiChar, Seq, StateSet, StateSetStatus};
use treetime_utils::interval::range::range_contains;

/// Backward pass: resolve candidate positions using minimum-change parsimony.
///
/// This is the authoritative resolution pass. For every candidate position it collects each
/// informative child's state set exactly once and combines them, writing the result into
/// `sequence` and the returned variable map.
///
/// Candidate positions come from two sources, which must both be supplied:
///
/// - positions that are variable in at least one child, i.e. where a child carries an IUPAC
///   ambiguity code or was itself left ambiguous by its own backward pass;
/// - `discovered`, the positions where children hold differing canonical states, found by
///   [`discover_fixed_disagreements_backward`]. These appear in no child's `fitch.variable`.
///
/// The combination rule is the intersection of the child state sets when that is non-empty, and
/// otherwise the plurality: the states held by the greatest number of children. Intersection is
/// already minimal when non-empty, since its members are contained in every child and so attain
/// the largest possible count. See [`StateSet::from_plurality`] for why the plurality, rather
/// than Fitch's union fallback, is required once a node has three or more children.
///
/// Positions transmitted along an edge or in `non_char` ranges are skipped.
pub fn resolve_variable_positions_backward(
  children: &[(&SparseSeqInfo, &SparseEdgePartition)],
  discovered: &[usize],
  non_char: &[(usize, usize)],
  sequence: &mut Seq,
) -> BTreeMap<usize, StateSet> {
  let variable_positions = children
    .iter()
    .flat_map(|(c, _)| c.fitch.variable.keys().copied())
    .chain(discovered.iter().copied())
    .unique()
    .collect_vec();

  let mut variable = BTreeMap::new();

  for pos in variable_positions {
    if range_contains(non_char, pos) {
      continue; // this node has no character state here, so nothing is variable
    }

    // Collect child profiles (1D vectors)
    let child_profiles = children
      .iter()
      .filter_map(|(child, edge)| {
        if let Some(transmission) = &edge.transmission {
          if range_contains(transmission, pos) {
            return None; // transmission field is not currently used
          }
        }
        if range_contains(&child.non_char, pos) {
          return None; // this position does not have character state information
        }
        let state = match child.fitch.variable.get(&pos) {
          Some(var_pos) => *var_pos,
          None => StateSet::from_char(child.sequence[pos]),
        };
        Some(state)
      })
      .collect_vec();

    if child_profiles.is_empty() {
      // No child carries character state here. The node's `non_char` is derived from the
      // children's, so this is unreachable via that route; guard rather than fall through and
      // store an empty state set, which no consumer can resolve.
      continue;
    }

    // Calculate minimum-change parsimony.
    // If we save the states of the children for each position that is variable in the node,
    // then we would not need the full sequences in the forward pass.
    let intersection = StateSet::from_intersection(&child_profiles);

    match intersection.get() {
      StateSetStatus::Unambiguous(state) => {
        // intersection has a single state, write it
        sequence[pos] = state;
      },
      StateSetStatus::Ambiguous(_) => {
        // more than one possible states
        variable.insert(pos, intersection);
        sequence[pos] = VARIABLE_CHAR;
      },
      StateSetStatus::Empty => {
        // No state is shared by every child, so retain those shared by the most. With one or
        // two children the plurality is exactly the union, which is Fitch's fallback; with three
        // or more the union can retain non-minimal states.
        let resolved = if child_profiles.len() <= 2 {
          StateSet::from_union(&child_profiles)
        } else {
          StateSet::from_plurality(&child_profiles)
        };
        variable.insert(pos, resolved);
        sequence[pos] = VARIABLE_CHAR;
      },
    }
  }

  variable
}

/// Backward pass: find positions where children hold differing canonical states.
///
/// This is a discovery pass only. It writes the agreed state into `sequence` where all children
/// agree, and marks a position with `VARIABLE_CHAR` as soon as any child disagrees. It does not
/// combine state sets: resolution belongs to [`resolve_variable_positions_backward`], which runs
/// afterwards and recomputes each candidate position from every child.
///
/// The returned positions are those marked `VARIABLE_CHAR`, in ascending order.
///
/// Detecting a disagreement only requires noticing that *some* child differs, never all children
/// at once, so this keeps the child-major loop over whole sequences. Comparing every child against
/// every position is `O(children · length)` and irreducible while nodes store dense sequences, so
/// the loop deliberately performs no map lookups or state-set construction. Positions are recorded
/// as they are flagged rather than by rescanning `sequence` afterwards, which would add another
/// full-length pass per node; the `VARIABLE_CHAR` skip guarantees each is recorded at most once.
pub fn discover_fixed_disagreements_backward(
  children: &[(&SparseSeqInfo, &SparseEdgePartition)],
  alphabet: &Alphabet,
  sequence: &mut Seq,
) -> Vec<usize> {
  let mut discovered = vec![];
  for &(child, _) in children {
    for (pos, parent_state) in sequence.iter_mut().enumerate() {
      let child_state = child.sequence[pos];
      if *parent_state == child_state || *parent_state == NON_CHAR || *parent_state == VARIABLE_CHAR {
        continue; // agrees with the parent, known non-char, or already flagged as disagreeing
      }
      if alphabet.is_canonical(child_state) {
        if *parent_state == FILL_CHAR {
          // first child to claim this position: adopt its state as the provisional agreement
          *parent_state = child_state;
        } else {
          // a canonical state differing from the provisional agreement: flag for resolution
          *parent_state = VARIABLE_CHAR;
          discovered.push(pos);
        }
      }
    }
  }

  // Flags are emitted in ascending order within each child but interleaved across children.
  // Sorting keeps the returned contract stable for callers and tests; it is over the handful of
  // disagreeing positions, not the sequence length.
  discovered.sort_unstable();
  discovered
}

/// Commit one state from a set of equally parsimonious alternatives.
///
/// Every member of a resolved state set yields the same subtree cost, so this choice is free with
/// respect to parsimony. It is kept deterministic so that output is reproducible between runs, and
/// ordered by the alphabet rather than by byte value; see [`Alphabet::first_canonical`].
///
/// Falls back to the lowest byte if the set holds no canonical state, which preserves the previous
/// behaviour for sets that consumers cannot interpret anyway.
fn choose_state(states: StateSet, alphabet: &Alphabet) -> AsciiChar {
  alphabet.first_canonical(states).unwrap_or_else(|| states.get_one())
}

/// Forward pass: resolve variable substitution states at the root.
///
/// The root has no parent to inherit from, so each variable position is committed by the
/// tie-break policy in [`choose_state`].
///
/// Variable indels at the root default to present (no gap). Direction is
/// resolved in the forward pass on children via parent state.
pub fn resolve_root_forward(
  sequence: &mut Seq,
  variable: &BTreeMap<usize, StateSet>,
  chosen_state: &mut BTreeMap<usize, AsciiChar>,
  alphabet: &Alphabet,
) {
  for (pos, states) in variable {
    let chosen = choose_state(*states, alphabet);
    sequence[*pos] = chosen;
    chosen_state.insert(*pos, chosen);
  }
}

/// Forward pass: resolve non-root variable positions, detecting substitutions.
///
/// For each variable position, picks the parent state if present in the child's
/// state set, otherwise picks one state and records a substitution. Also detects
/// parent-only variable positions that introduce mutations.
///
/// Returns the list of substitutions for this edge.
pub fn resolve_nonroot_substitutions_forward(
  sequence: &mut Seq,
  gaps: &[(usize, usize)],
  variable: &mut BTreeMap<usize, StateSet>,
  chosen_state: &mut BTreeMap<usize, AsciiChar>,
  composition: &mut Composition,
  parent_seq: &SparseSeqInfo,
  alphabet: &Alphabet,
) -> Result<Vec<Sub>, Report> {
  let mut subs = vec![];

  // for each variable position, pick a state or a mutation
  for (pos, states) in variable.iter_mut() {
    if range_contains(gaps, *pos) {
      continue; // deleted on this edge: there is no character state to resolve or mutate
    }
    let pnuc = parent_seq.sequence[*pos];
    if alphabet.is_canonical(pnuc) {
      // check whether parent is in child profile (sum>0 --> parent state is in profile)
      if states.contains(pnuc) {
        // Preferring the parent's state is not a tie-break heuristic: given a minimal state set it
        // is provably optimal, since matching the parent avoids a change on this edge while any
        // other member of the set costs exactly one more. This is what pushes mutations toward the
        // tips.
        sequence[*pos] = pnuc;
      } else {
        let cnuc = choose_state(*states, alphabet);
        sequence[*pos] = cnuc;
        let m = Sub::new(pnuc, *pos, cnuc)?;
        m.check_determined(alphabet)?;
        composition.add_sub(&m);
        subs.push(m);
      }
    } else if alphabet.is_gap(pnuc) && !range_contains(gaps, *pos) {
      // if parent is gap, but child isn't, we need to resolve variable states
      sequence[*pos] = choose_state(*states, alphabet);
    }
    chosen_state.insert(*pos, sequence[*pos]);
  }

  for &pos in parent_seq.fitch.variable.keys() {
    if variable.contains_key(&pos) || range_contains(&parent_seq.gaps, pos) || range_contains(gaps, pos) {
      continue;
    }

    // NOTE: access to full_seq would not be necessary if we had saved the
    // child state of variable positions in the backward pass
    let node_nuc = sequence[pos];
    if alphabet.is_canonical(node_nuc) && parent_seq.sequence[pos] != node_nuc {
      let m = Sub::new(parent_seq.sequence[pos], pos, node_nuc)?;
      m.check_determined(alphabet)?;
      composition.add_sub(&m);
      subs.push(m);
    }
  }

  // Sort by position: the two loops above (child variable, then parent-only variable)
  // can emit positions out of order when parent-only positions precede child positions.
  subs.sort();
  Ok(subs)
}

/// Forward pass: fill gaps and unknown positions in the sequence, finalize root composition.
pub fn finalize_sequence_forward(
  sequence: &mut Seq,
  gaps: &[(usize, usize)],
  unknown: &[(usize, usize)],
  composition: &mut Composition,
  alphabet: &Alphabet,
  is_root: bool,
) {
  // fill in the gapped positions. this is done for all nodes, including the root, the composition of non-root nodes is already correct
  for r in gaps {
    sequence[r.0..r.1].fill(alphabet.gap());
  }
  for r in unknown {
    // composition is already adjusted
    sequence[r.0..r.1].fill(alphabet.unknown());
  }
  if is_root {
    // if the node is the root, the composition is calculated from the full sequence
    *composition = Composition::with_seq(sequence.as_slice(), alphabet.chars(), alphabet.gap());
  }
}
