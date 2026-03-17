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

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_abs_diff_eq;
  use ndarray::array;

  #[test]
  fn test_logsumexp_normalize_equal_log_probs() {
    // ln(0.25) for each of 4 states: already-normalized probabilities in log space.
    // sum(exp(log_vec)) = 4 * 0.25 = 1.0, so log_norm = ln(1.0) = 0.0.
    let log_vec = Array1::from_elem(4, 0.25_f64.ln());
    let (normalized, log_norm) = logsumexp_normalize(&log_vec);

    helpers::assert_valid_distribution(&normalized);
    assert_abs_diff_eq!(normalized, array![0.25, 0.25, 0.25, 0.25], epsilon = 1e-10);
    assert_abs_diff_eq!(log_norm, 0.0, epsilon = 1e-10);
  }

  #[test]
  fn test_logsumexp_normalize_descending() {
    // Unnormalized log-values [0, -1, -2, -3]. The logsumexp trick subtracts max=0,
    // exponentiates to [1, e^-1, e^-2, e^-3], and divides by the sum S.
    let log_vec = array![0.0, -1.0, -2.0, -3.0];
    let (normalized, log_norm) = logsumexp_normalize(&log_vec);

    let sum = 1.0 + (-1.0_f64).exp() + (-2.0_f64).exp() + (-3.0_f64).exp();
    let expected = array![
      1.0 / sum,
      (-1.0_f64).exp() / sum,
      (-2.0_f64).exp() / sum,
      (-3.0_f64).exp() / sum
    ];

    helpers::assert_valid_distribution(&normalized);
    assert_abs_diff_eq!(normalized, expected, epsilon = 1e-10);
    assert_abs_diff_eq!(log_norm, sum.ln(), epsilon = 1e-10);
  }

  #[test]
  fn test_logsumexp_normalize_single_dominant_state() {
    // One state at 0.0, others at -100.0. exp(-100) ~ 3.7e-44, so the dominant
    // state concentrates nearly all probability mass.
    let log_vec = array![0.0, -100.0, -100.0, -100.0];
    let (normalized, log_norm) = logsumexp_normalize(&log_vec);

    helpers::assert_valid_distribution(&normalized);
    assert_abs_diff_eq!(normalized, array![1.0, 0.0, 0.0, 0.0], epsilon = 1e-10);
    assert_abs_diff_eq!(log_norm, 0.0, epsilon = 1e-10);
  }

  #[test]
  fn test_logsumexp_normalize_all_neg_inf() {
    // All states have zero probability (log = -inf). The function falls back
    // to a uniform distribution with log_norm = -inf.
    let log_vec = Array1::from_elem(4, f64::NEG_INFINITY);
    let (normalized, log_norm) = logsumexp_normalize(&log_vec);

    assert_abs_diff_eq!(normalized, array![0.25, 0.25, 0.25, 0.25], epsilon = 1e-10);
    assert!(
      log_norm == f64::NEG_INFINITY,
      "log_norm should be NEG_INFINITY, got {log_norm}"
    );
  }

  #[test]
  fn test_logsumexp_normalize_all_neg_inf_three_states() {
    // All-inf fallback with 3 states: uniform = 1/3 each.
    let log_vec = Array1::from_elem(3, f64::NEG_INFINITY);
    let (normalized, log_norm) = logsumexp_normalize(&log_vec);

    assert_abs_diff_eq!(normalized, Array1::from_elem(3, 1.0 / 3.0), epsilon = 1e-10);
    assert!(
      log_norm == f64::NEG_INFINITY,
      "log_norm should be NEG_INFINITY, got {log_norm}"
    );
  }

  #[test]
  fn test_logsumexp_normalize_mixed_finite_neg_inf() {
    // States 0 and 2 have finite log-probs, states 1 and 3 are -inf.
    // exp(-inf) = 0, so only states 0 and 2 contribute.
    // Equivalent to softmax([0, -1]) distributed across 4 positions.
    let log_vec = array![0.0, f64::NEG_INFINITY, -1.0, f64::NEG_INFINITY];
    let (normalized, log_norm) = logsumexp_normalize(&log_vec);

    let sum = 1.0 + (-1.0_f64).exp();
    let expected = array![1.0 / sum, 0.0, (-1.0_f64).exp() / sum, 0.0];

    helpers::assert_valid_distribution(&normalized);
    assert_abs_diff_eq!(normalized, expected, epsilon = 1e-10);
    assert_abs_diff_eq!(log_norm, sum.ln(), epsilon = 1e-10);
  }

  #[test]
  fn test_logsumexp_normalize_single_finite_rest_neg_inf() {
    // Only one finite state: gets all probability mass.
    let log_vec = array![f64::NEG_INFINITY, -3.0, f64::NEG_INFINITY, f64::NEG_INFINITY];
    let (normalized, log_norm) = logsumexp_normalize(&log_vec);

    helpers::assert_valid_distribution(&normalized);
    assert_abs_diff_eq!(normalized, array![0.0, 1.0, 0.0, 0.0], epsilon = 1e-10);
    assert_abs_diff_eq!(log_norm, -3.0, epsilon = 1e-10);
  }

  #[test]
  fn test_logsumexp_normalize_large_negative_values() {
    // All values around -1000. Without the logsumexp trick, exp(-1000) underflows
    // to 0. The trick subtracts max=-1000, leaving [0, -1, -2, -3], producing
    // the same relative probabilities as the descending case.
    let log_vec = array![-1000.0, -1001.0, -1002.0, -1003.0];
    let (normalized, log_norm) = logsumexp_normalize(&log_vec);

    let reference = array![0.0, -1.0, -2.0, -3.0];
    let (reference_normalized, reference_log_norm) = logsumexp_normalize(&reference);

    helpers::assert_valid_distribution(&normalized);
    assert_abs_diff_eq!(normalized, reference_normalized, epsilon = 1e-10);
    assert_abs_diff_eq!(log_norm, reference_log_norm - 1000.0, epsilon = 1e-10);
  }

  #[test]
  fn test_logsumexp_normalize_large_positive_values() {
    // Large positive values: exp(1000) overflows without the logsumexp trick.
    // Same relative ratios as [0, -1, -2, -3] shifted by +1003.
    let log_vec = array![1003.0, 1002.0, 1001.0, 1000.0];
    let (normalized, log_norm) = logsumexp_normalize(&log_vec);

    let reference = array![0.0, -1.0, -2.0, -3.0];
    let (reference_normalized, reference_log_norm) = logsumexp_normalize(&reference);

    helpers::assert_valid_distribution(&normalized);
    assert_abs_diff_eq!(normalized, reference_normalized, epsilon = 1e-10);
    assert_abs_diff_eq!(log_norm, reference_log_norm + 1003.0, epsilon = 1e-10);
  }

  #[test]
  fn test_logsumexp_normalize_single_state() {
    // Single-element vector: the only state gets probability 1.0,
    // log_norm equals the input value.
    let log_vec = array![5.0];
    let (normalized, log_norm) = logsumexp_normalize(&log_vec);

    assert_abs_diff_eq!(normalized, array![1.0], epsilon = 1e-10);
    assert_abs_diff_eq!(log_norm, 5.0, epsilon = 1e-10);
  }

  #[test]
  fn test_logsumexp_normalize_two_equal_states() {
    // Two equal log-values: each state gets 0.5, log_norm = value + ln(2).
    let log_vec = array![0.0, 0.0];
    let (normalized, log_norm) = logsumexp_normalize(&log_vec);

    assert_abs_diff_eq!(normalized, array![0.5, 0.5], epsilon = 1e-10);
    assert_abs_diff_eq!(log_norm, 2.0_f64.ln(), epsilon = 1e-10);
  }

  mod helpers {
    use approx::assert_abs_diff_eq;
    use ndarray::Array1;

    pub fn assert_valid_distribution(normalized: &Array1<f64>) {
      assert!(
        normalized.iter().all(|&v| v >= 0.0),
        "all probabilities must be non-negative: {normalized}"
      );
      assert!(
        normalized.iter().all(|&v| v.is_finite()),
        "all probabilities must be finite: {normalized}"
      );
      assert_abs_diff_eq!(normalized.sum(), 1.0, epsilon = 1e-10);
    }
  }
}
