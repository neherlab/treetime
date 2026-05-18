use crate::alphabet::alphabet::{Alphabet, FILL_CHAR, NON_CHAR, VARIABLE_CHAR};
use crate::representation::payload::sparse::{Deletion, SparseEdgePartition, SparseSeqInfo};
use crate::seq::composition::Composition;
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use treetime_primitives::{AlphabetLike, AsciiChar, BitSet128, Seq, StateSet, StateSetStatus, stateset};
use treetime_utils::interval::range::range_contains;

/// Backward pass: resolve variable positions using Fitch parsimony.
///
/// For each position that is variable in at least one child, collects the child
/// state sets, computes the intersection (or union when the intersection is
/// empty), and writes the result into `sequence` and the returned variable map.
///
/// Positions transmitted along an edge or in `non_char` ranges are skipped.
pub fn resolve_variable_positions_backward(
  children: &[(&SparseSeqInfo, &SparseEdgePartition)],
  non_char: &[(usize, usize)],
  sequence: &mut Seq,
) -> BTreeMap<usize, StateSet> {
  let variable_positions = children
    .iter()
    .flat_map(|(c, _)| c.fitch.variable.keys().copied())
    .unique()
    .collect_vec();

  let mut variable = BTreeMap::new();

  for pos in variable_positions {
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

    // Calculate Fitch parsimony.
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
        let union = StateSet::from_union(&child_profiles);
        variable.insert(pos, union);
        sequence[pos] = VARIABLE_CHAR;
      },
    }
  }

  variable
}

/// Backward pass: resolve fixed positions where children disagree with current parent state.
///
/// Scans each child's fixed sequence positions. When a child has a canonical state
/// that differs from the current parent state, either sets the parent (if still
/// FILL_CHAR) or promotes the position to variable.
pub fn resolve_fixed_positions_backward(
  children: &[(&SparseSeqInfo, &SparseEdgePartition)],
  alphabet: &Alphabet,
  sequence: &mut Seq,
  variable: &mut BTreeMap<usize, StateSet>,
) {
  for &(child, _) in children {
    for (pos, parent_state) in sequence.iter_mut().enumerate() {
      let child_state = child.sequence[pos];
      if *parent_state == child_state || *parent_state == NON_CHAR {
        continue; // if parent is equal to child state or we know it's a non-char, skip
      }
      if alphabet.is_canonical(child_state) {
        if *parent_state == FILL_CHAR {
          // if child state is canonical and parent is still FILL_CHAR, set parent_state
          *parent_state = child_state;
        } else {
          // otherwise set or update the variable state
          *variable.entry(pos).or_insert_with(|| stateset! {*parent_state}) += child_state;
          *parent_state = VARIABLE_CHAR;
        }
      }
    }
  }
}

/// Forward pass: resolve variable states and indels at the root by majority rule.
pub fn resolve_root_forward(
  sequence: &mut Seq,
  gaps: &mut Vec<(usize, usize)>,
  variable: &BTreeMap<usize, StateSet>,
  variable_indel: &BTreeMap<(usize, usize), Deletion>,
  chosen_state: &mut BTreeMap<usize, AsciiChar>,
) {
  for (pos, states) in variable {
    let chosen = states.get_one();
    sequence[*pos] = chosen;
    chosen_state.insert(*pos, chosen);
  }
  // process indels as majority rule at the root
  for (r, indel) in variable_indel {
    if indel.deleted > indel.present {
      gaps.push(*r);
    }
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
    *composition = Composition::with_sequence(sequence.iter().copied(), alphabet.chars(), alphabet.gap());
  }
}
