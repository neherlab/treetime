use crate::alphabet::alphabet::{Alphabet, NON_CHAR, VARIABLE_CHAR};
use crate::partition::sparse::{SparseEdgePartition, SparseSeqInfo};
use crate::seq::composition::Composition;
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use treetime_primitives::{AlphabetLike, AsciiChar, Seq, StateSet, StateSetStatus};
use treetime_utils::interval::range::range_contains;

/// Resolve substitution-informative positions during the Fitch backward pass.
///
/// The candidate list is complete for substitutions because invariant canonical
/// columns have one state across all leaves. At each candidate, intersect all
/// informative child state sets; an empty intersection retains their union,
/// preserving the existing unordered-character recurrence.
pub fn resolve_informative_positions_backward(
  children: &[(&SparseSeqInfo, &SparseEdgePartition)],
  informative_positions: &[usize],
  sequence: &mut Seq,
) -> BTreeMap<usize, StateSet> {
  let mut variable = BTreeMap::new();

  for &pos in informative_positions {
    if sequence[pos] == NON_CHAR {
      continue;
    }

    let child_profiles = children
      .iter()
      .filter_map(|(child, edge)| {
        if let Some(transmission) = &edge.transmission {
          if range_contains(transmission, pos) {
            return None;
          }
        }
        if range_contains(&child.non_char, pos) {
          return None;
        }
        let state = match child.fitch.variable.get(&pos) {
          Some(var_pos) => *var_pos,
          None => StateSet::from_char(child.sequence[pos]),
        };
        Some(state)
      })
      .collect_vec();

    let intersection = StateSet::from_intersection(&child_profiles);

    match intersection.get() {
      StateSetStatus::Unambiguous(state) => {
        sequence[pos] = state;
      },
      StateSetStatus::Ambiguous(_) => {
        variable.insert(pos, intersection);
        sequence[pos] = VARIABLE_CHAR;
      },
      StateSetStatus::Empty => {
        let union = StateSet::from_union(&child_profiles);
        variable.insert(pos, union);
        sequence[pos] = VARIABLE_CHAR;
      },
    }
  }

  variable
}

/// Forward pass: resolve variable substitution states at the root.
///
/// Variable indels at the root default to present (no gap). Direction is
/// resolved in the forward pass on children via parent state.
pub fn resolve_root_forward(
  sequence: &mut Seq,
  variable: &BTreeMap<usize, StateSet>,
  chosen_state: &mut BTreeMap<usize, AsciiChar>,
) {
  for (pos, states) in variable {
    let chosen = states.get_one();
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
    let pnuc = parent_seq.sequence[*pos];
    if alphabet.is_canonical(pnuc) {
      // check whether parent is in child profile (sum>0 --> parent state is in profile)
      if states.contains(pnuc) {
        sequence[*pos] = pnuc;
      } else {
        let cnuc = states.get_one();
        sequence[*pos] = cnuc;
        let m = Sub::new(pnuc, *pos, cnuc)?;
        m.check_determined(alphabet)?;
        composition.add_sub(&m);
        subs.push(m);
      }
    } else if alphabet.is_gap(pnuc) && !range_contains(gaps, *pos) {
      // if parent is gap, but child isn't, we need to resolve variable states
      sequence[*pos] = states.get_one();
    }
    chosen_state.insert(*pos, sequence[*pos]);
  }

  for &pos in parent_seq.fitch.variable.keys() {
    if variable.contains_key(&pos) || range_contains(&parent_seq.gaps, pos) {
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
