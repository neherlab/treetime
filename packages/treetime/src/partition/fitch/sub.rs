use crate::alphabet::alphabet::{Alphabet, FILL_CHAR, NON_CHAR, VARIABLE_CHAR};
use crate::partition::storage::sparse::{FitchSeqInfo, SparseEdgeObs};
use crate::seq::composition::Composition;
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use treetime_primitives::{AlphabetLike, AsciiChar, Seq, StateSet, StateSetStatus};
use treetime_utils::interval::range::range_contains;

pub(crate) fn resolve_variable_positions_backward(
  children: &[(&FitchSeqInfo, &SparseEdgeObs)],
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

    if child_profiles.is_empty() {
      continue;
    }

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

pub(crate) fn discover_fixed_disagreements_backward(
  children: &[(&FitchSeqInfo, &SparseEdgeObs)],
  alphabet: &Alphabet,
  sequence: &mut Seq,
) -> Vec<usize> {
  let mut discovered = vec![];
  for &(child, _) in children {
    for (pos, parent_state) in sequence.iter_mut().enumerate() {
      let child_state = child.sequence[pos];
      if *parent_state == child_state || *parent_state == NON_CHAR || *parent_state == VARIABLE_CHAR {
        continue;
      }
      if alphabet.is_canonical(child_state) {
        if *parent_state == FILL_CHAR {
          *parent_state = child_state;
        } else {
          *parent_state = VARIABLE_CHAR;
          discovered.push(pos);
        }
      }
    }
  }

  discovered.sort_unstable();
  discovered
}

pub(crate) fn resolve_root_forward(
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

pub(crate) fn resolve_nonroot_substitutions_forward(
  sequence: &mut Seq,
  gaps: &[(usize, usize)],
  variable: &mut BTreeMap<usize, StateSet>,
  chosen_state: &mut BTreeMap<usize, AsciiChar>,
  composition: &mut Composition,
  parent_seq: &FitchSeqInfo,
  alphabet: &Alphabet,
) -> Result<Vec<Sub>, Report> {
  let mut subs = vec![];

  for (pos, states) in variable.iter_mut() {
    if range_contains(gaps, *pos) {
      continue;
    }
    let pnuc = parent_seq.sequence[*pos];
    if alphabet.is_canonical(pnuc) {
      if states.contains(pnuc) {
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
      sequence[*pos] = choose_state(*states, alphabet);
    }
    chosen_state.insert(*pos, sequence[*pos]);
  }

  for &pos in parent_seq.fitch.variable.keys() {
    if variable.contains_key(&pos) || range_contains(&parent_seq.gaps, pos) || range_contains(gaps, pos) {
      continue;
    }

    let node_nuc = sequence[pos];
    if alphabet.is_canonical(node_nuc) && parent_seq.sequence[pos] != node_nuc {
      let m = Sub::new(parent_seq.sequence[pos], pos, node_nuc)?;
      m.check_determined(alphabet)?;
      composition.add_sub(&m);
      subs.push(m);
    }
  }

  subs.sort();
  Ok(subs)
}

fn choose_state(states: StateSet, alphabet: &Alphabet) -> AsciiChar {
  alphabet.first_canonical(states).unwrap_or_else(|| states.get_one())
}

pub(crate) fn finalize_sequence_forward(
  sequence: &mut Seq,
  gaps: &[(usize, usize)],
  unknown: &[(usize, usize)],
  composition: &mut Composition,
  alphabet: &Alphabet,
  is_root: bool,
) {
  for r in gaps {
    sequence[r.0..r.1].fill(alphabet.gap());
  }
  for r in unknown {
    sequence[r.0..r.1].fill(alphabet.unknown());
  }
  if is_root {
    *composition = Composition::with_seq(sequence.as_slice(), alphabet.chars(), alphabet.gap());
  }
}
