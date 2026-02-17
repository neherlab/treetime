use crate::alphabet::alphabet::Alphabet;
use crate::representation::payload::sparse::{MarginalSparseSeqDistribution, VarPos};
use crate::seq::composition::Composition;
use eyre::Report;
use maplit::btreemap;
use ndarray::{Array1, Array2};
use ndarray_stats::QuantileExt;
use std::collections::BTreeMap;
use std::iter::zip;
use treetime_primitives::AsciiChar;
use treetime_utils::interval::range::range_contains;

pub const EPS: f64 = 1e-4;

pub fn combine_messages(
  composition: &Composition,
  messages: &[MarginalSparseSeqDistribution],
  variable_pos: &BTreeMap<usize, AsciiChar>,
  reference_states: &[BTreeMap<usize, AsciiChar>],
  alphabet: &Alphabet,
  gtr_weight: Option<&Array1<f64>>,
) -> Result<MarginalSparseSeqDistribution, Report> {
  let mut seq_dis = MarginalSparseSeqDistribution {
    variable: btreemap! {},
    variable_indel: btreemap! {},
    fixed: btreemap! {},
    fixed_counts: composition.clone(),
    log_lh: messages.iter().map(|m| m.log_lh).sum(),
  };

  let mut fixed_counts = composition
    .counts()
    .iter()
    .map(|(k, v)| (*k, *v as f64))
    .collect::<BTreeMap<_, _>>();

  for (&pos, &state) in variable_pos {
    let mut all_states_equal = true;
    let mut vec = if let Some(gtr_weight) = gtr_weight {
      gtr_weight.clone()
    } else {
      Array1::from_elem(alphabet.n_canonical(), 1.0)
    };

    for (msg, states) in zip(messages, reference_states) {
      if let Some(var) = msg.variable.get(&pos) {
        vec *= &var.dis;
        if var.state != state {
          all_states_equal = false;
        }
      } else if let Some(ref_state) = states.get(&pos) {
        if alphabet.is_canonical(*ref_state) {
          vec *= &msg.fixed[ref_state];
        }
        if ref_state != &state {
          all_states_equal = false;
        }
      } else {
        vec *= &msg.fixed[&state];
      }
    }

    let vec_norm = vec.sum();
    seq_dis.log_lh += vec_norm.ln();
    if let Some(count) = fixed_counts.get_mut(&state) {
      *count -= 1.0;
    }

    if (*vec.max()? < (1.0 - EPS) * vec_norm) || !all_states_equal {
      seq_dis.fixed_counts.adjust_count(state, -1);
      let dis = vec / vec_norm;
      seq_dis.variable.insert(pos, VarPos { dis, state });
    }
  }

  for state in alphabet.canonical() {
    let mut vec = if let Some(gtr_weight) = gtr_weight {
      gtr_weight.clone()
    } else {
      Array1::from_elem(alphabet.n_canonical(), 1.0)
    };

    for msg in messages {
      vec *= &msg.fixed[&state];
    }
    let vec_norm = vec.sum();

    seq_dis.log_lh += fixed_counts[&state] * vec_norm.ln();
    seq_dis.fixed.insert(state, vec / vec_norm);
  }
  Ok(seq_dis)
}

pub fn propagate_raw(
  exp_qt: &Array2<f64>,
  seq_dis: &MarginalSparseSeqDistribution,
  transmission: Option<&[(usize, usize)]>,
) -> MarginalSparseSeqDistribution {
  let mut message = MarginalSparseSeqDistribution {
    variable: btreemap! {},
    variable_indel: btreemap! {},
    fixed: btreemap! {},
    fixed_counts: seq_dis.fixed_counts.clone(),
    log_lh: seq_dis.log_lh,
  };
  for (pos, state) in &seq_dis.variable {
    if let Some(transmission) = &transmission {
      if !range_contains(transmission, *pos) {
        continue;
      }
    }

    let dis = exp_qt.dot(&state.dis);
    let child_state = state.state;
    message.variable.insert(
      *pos,
      VarPos {
        dis,
        state: child_state,
      },
    );
  }

  for (&s, p) in &seq_dis.fixed {
    message.fixed.insert(s, exp_qt.dot(p));
  }

  message
}
