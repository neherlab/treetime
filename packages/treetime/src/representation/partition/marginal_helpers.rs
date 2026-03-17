use crate::alphabet::alphabet::Alphabet;
use crate::representation::payload::sparse::{MarginalSparseSeqDistribution, VarPos};
use crate::seq::composition::Composition;
use eyre::Report;
use maplit::btreemap;
use ndarray::{Array1, Array2};
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

  let n_states = alphabet.n_canonical();
  let initial_log = gtr_weight.map_or_else(|| Array1::zeros(n_states), |w| w.mapv(f64::ln));

  for (&pos, &state) in variable_pos {
    let mut all_states_equal = true;
    let mut log_vec = initial_log.clone();

    for (msg, states) in zip(messages, reference_states) {
      if let Some(var) = msg.variable.get(&pos) {
        log_vec.zip_mut_with(&var.dis, |lv, &p| *lv += p.ln());
        if var.state != state {
          all_states_equal = false;
        }
      } else if let Some(ref_state) = states.get(&pos) {
        if alphabet.is_canonical(*ref_state) {
          log_vec.zip_mut_with(&msg.fixed[ref_state], |lv, &p| *lv += p.ln());
        }
        if ref_state != &state {
          all_states_equal = false;
        }
      } else {
        log_vec.zip_mut_with(&msg.fixed[&state], |lv, &p| *lv += p.ln());
      }
    }

    let (dis, log_norm) = logsumexp_normalize(&log_vec);
    seq_dis.log_lh += log_norm;
    if let Some(count) = fixed_counts.get_mut(&state) {
      *count -= 1.0;
    }

    let max_prob = dis.iter().copied().fold(0.0_f64, f64::max);
    if max_prob < (1.0 - EPS) || !all_states_equal {
      seq_dis.fixed_counts.adjust_count(state, -1);
      seq_dis.variable.insert(pos, VarPos { dis, state });
    }
  }

  for state in alphabet.canonical() {
    let mut log_vec = initial_log.clone();

    for msg in messages {
      log_vec.zip_mut_with(&msg.fixed[&state], |lv, &p| *lv += p.ln());
    }

    let (dis, log_norm) = logsumexp_normalize(&log_vec);
    seq_dis.log_lh += fixed_counts[&state] * log_norm;
    seq_dis.fixed.insert(state, dis);
  }
  Ok(seq_dis)
}

/// Normalize a log-probability vector using the logsumexp trick.
///
/// Subtracts the maximum log-value before exponentiating to prevent underflow
/// when combining many small probabilities. Returns the normalized probability
/// vector and log(sum(exp(log_vec))).
///
/// When all entries are -inf (all states have zero probability), returns a
/// uniform distribution with log_norm = NEG_INFINITY.
fn logsumexp_normalize(log_vec: &Array1<f64>) -> (Array1<f64>, f64) {
  let max_val = log_vec.iter().copied().fold(f64::NEG_INFINITY, f64::max);

  if !max_val.is_finite() {
    let n = log_vec.len();
    return (Array1::from_elem(n, 1.0 / n as f64), f64::NEG_INFINITY);
  }

  let shifted = log_vec.mapv(|v| (v - max_val).exp());
  let shifted_sum = shifted.sum();
  let log_norm = max_val + shifted_sum.ln();
  let normalized = shifted / shifted_sum;

  (normalized, log_norm)
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
