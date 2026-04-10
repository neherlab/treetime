use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::representation::payload::sparse::{MarginalSparseSeqDistribution, VarPos};
use crate::seq::composition::Composition;
use eyre::Report;
use maplit::btreemap;
use ndarray::{Array1, Array2, ArrayView1};
use std::collections::BTreeMap;
use std::iter::zip;
#[cfg(test)]
use treetime_primitives::AlphabetLike;
use treetime_primitives::AsciiChar;
use treetime_utils::array::ndarray::argmax_first;
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

    let (dis, log_norm) = logsumexp_normalize(log_vec.view());
    seq_dis.log_lh += log_norm;
    if let Some(count) = fixed_counts.get_mut(&state) {
      *count -= 1.0;
    }

    let map_state = alphabet.char(argmax_first(&dis.view()).unwrap_or(0));
    let max_prob = dis.iter().copied().fold(0.0_f64, f64::max);
    if max_prob < (1.0 - EPS) || !all_states_equal || map_state != state {
      seq_dis.fixed_counts.adjust_count(state, -1);
      seq_dis.variable.insert(pos, VarPos { dis, state });
    }
  }

  for state in alphabet.canonical() {
    let mut log_vec = initial_log.clone();

    for msg in messages {
      log_vec.zip_mut_with(&msg.fixed[&state], |lv, &p| *lv += p.ln());
    }

    let (dis, log_norm) = logsumexp_normalize(log_vec.view());
    let multiplicity = fixed_counts[&state];
    // Zero fixed-site multiplicity means this shared row corresponds to no
    // actual sites. Its contribution is therefore exactly zero, even if the
    // row itself underflowed and `log_norm` is `NEG_INFINITY`.
    if multiplicity != 0.0 {
      seq_dis.log_lh += multiplicity * log_norm;
    }
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
pub fn logsumexp_normalize(log_vec: ArrayView1<'_, f64>) -> (Array1<f64>, f64) {
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

/// Propagate sparse message with per-site rate variation.
///
/// For variable positions, computes a position-specific P(t) using the site rate
/// at that position. For fixed positions, uses the default scalar mu rate. This is
/// an approximation: all fixed positions with the same state share one propagated
/// profile, so per-site rate variation cannot be represented for them. The
/// approximation is exact when rates are uniform. See
/// `docs/port-intentional-changes/sparse-fixed-position-scalar-rate-approximation.md`.
///
/// `transpose`: when true, uses P(t)^T (backward pass: child -> parent).
///              when false, uses P(t) directly (forward pass: parent -> child).
pub fn propagate_raw_per_site(
  gtr: &GTR,
  branch_length: f64,
  transpose: bool,
  seq_dis: &MarginalSparseSeqDistribution,
  transmission: Option<&[(usize, usize)]>,
) -> MarginalSparseSeqDistribution {
  let site_rates = gtr
    .site_rates
    .as_ref()
    .expect("propagate_raw_per_site requires site_rates");
  let default_exp_qt = if transpose {
    gtr.expQt(branch_length).t().to_owned()
  } else {
    gtr.expQt(branch_length)
  };

  let mut message = MarginalSparseSeqDistribution {
    variable: btreemap! {},
    variable_indel: btreemap! {},
    fixed: btreemap! {},
    fixed_counts: seq_dis.fixed_counts.clone(),
    log_lh: seq_dis.log_lh,
  };

  for (&pos, state) in &seq_dis.variable {
    if let Some(transmission) = &transmission {
      if !range_contains(transmission, pos) {
        continue;
      }
    }

    let rate = site_rates[pos];
    let exp_qt_pos = if transpose {
      gtr.expQt_with_rate(branch_length, rate).t().to_owned()
    } else {
      gtr.expQt_with_rate(branch_length, rate)
    };
    let dis = exp_qt_pos.dot(&state.dis);
    message.variable.insert(
      pos,
      VarPos {
        dis,
        state: state.state,
      },
    );
  }

  // Fixed positions: use default scalar mu rate
  for (&s, p) in &seq_dis.fixed {
    message.fixed.insert(s, default_exp_qt.dot(p));
  }

  message
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::seq::composition::Composition;
  use approx::assert_abs_diff_eq;
  use ndarray::array;

  #[test]
  fn test_logsumexp_normalize_equal_log_probs() {
    // ln(0.25) for each of 4 states: already-normalized probabilities in log space.
    // sum(exp(log_vec)) = 4 * 0.25 = 1.0, so log_norm = ln(1.0) = 0.0.
    let log_vec = Array1::from_elem(4, 0.25_f64.ln());
    let (normalized, log_norm) = logsumexp_normalize(log_vec.view());

    helpers::assert_valid_distribution(&normalized);
    assert_abs_diff_eq!(normalized, array![0.25, 0.25, 0.25, 0.25], epsilon = 1e-10);
    assert_abs_diff_eq!(log_norm, 0.0, epsilon = 1e-10);
  }

  #[test]
  fn test_logsumexp_normalize_descending() {
    // Unnormalized log-values [0, -1, -2, -3]. The logsumexp trick subtracts max=0,
    // exponentiates to [1, e^-1, e^-2, e^-3], and divides by the sum S.
    let log_vec = array![0.0, -1.0, -2.0, -3.0];
    let (normalized, log_norm) = logsumexp_normalize(log_vec.view());

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
    let (normalized, log_norm) = logsumexp_normalize(log_vec.view());

    helpers::assert_valid_distribution(&normalized);
    assert_abs_diff_eq!(normalized, array![1.0, 0.0, 0.0, 0.0], epsilon = 1e-10);
    assert_abs_diff_eq!(log_norm, 0.0, epsilon = 1e-10);
  }

  #[test]
  fn test_logsumexp_normalize_all_neg_inf() {
    // All states have zero probability (log = -inf). The function falls back
    // to a uniform distribution with log_norm = -inf.
    let log_vec = Array1::from_elem(4, f64::NEG_INFINITY);
    let (normalized, log_norm) = logsumexp_normalize(log_vec.view());

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
    let (normalized, log_norm) = logsumexp_normalize(log_vec.view());

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
    let (normalized, log_norm) = logsumexp_normalize(log_vec.view());

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
    let (normalized, log_norm) = logsumexp_normalize(log_vec.view());

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
    let (normalized, log_norm) = logsumexp_normalize(log_vec.view());

    let reference = array![0.0, -1.0, -2.0, -3.0];
    let (reference_normalized, reference_log_norm) = logsumexp_normalize(reference.view());

    helpers::assert_valid_distribution(&normalized);
    assert_abs_diff_eq!(normalized, reference_normalized, epsilon = 1e-10);
    assert_abs_diff_eq!(log_norm, reference_log_norm - 1000.0, epsilon = 1e-10);
  }

  #[test]
  fn test_logsumexp_normalize_large_positive_values() {
    // Large positive values: exp(1000) overflows without the logsumexp trick.
    // Same relative ratios as [0, -1, -2, -3] shifted by +1003.
    let log_vec = array![1003.0, 1002.0, 1001.0, 1000.0];
    let (normalized, log_norm) = logsumexp_normalize(log_vec.view());

    let reference = array![0.0, -1.0, -2.0, -3.0];
    let (reference_normalized, reference_log_norm) = logsumexp_normalize(reference.view());

    helpers::assert_valid_distribution(&normalized);
    assert_abs_diff_eq!(normalized, reference_normalized, epsilon = 1e-10);
    assert_abs_diff_eq!(log_norm, reference_log_norm + 1003.0, epsilon = 1e-10);
  }

  #[test]
  fn test_logsumexp_normalize_single_state() {
    // Single-element vector: the only state gets probability 1.0,
    // log_norm equals the input value.
    let log_vec = array![5.0];
    let (normalized, log_norm) = logsumexp_normalize(log_vec.view());

    assert_abs_diff_eq!(normalized, array![1.0], epsilon = 1e-10);
    assert_abs_diff_eq!(log_norm, 5.0, epsilon = 1e-10);
  }

  #[test]
  fn test_logsumexp_normalize_two_equal_states() {
    // Two equal log-values: each state gets 0.5, log_norm = value + ln(2).
    let log_vec = array![0.0, 0.0];
    let (normalized, log_norm) = logsumexp_normalize(log_vec.view());

    assert_abs_diff_eq!(normalized, array![0.5, 0.5], epsilon = 1e-10);
    assert_abs_diff_eq!(log_norm, 2.0_f64.ln(), epsilon = 1e-10);
  }

  #[test]
  fn test_propagate_raw_per_site_forward() {
    use crate::gtr::get_gtr::{JC69Params, jc69};

    let a = AsciiChar::from_byte_unchecked(b'A');
    let c = AsciiChar::from_byte_unchecked(b'C');
    let g = AsciiChar::from_byte_unchecked(b'G');

    let mut gtr = jc69(JC69Params::default()).unwrap();
    gtr.set_site_rates(Array1::from_vec(vec![
      0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 0.1,
    ]));

    let t = 0.3;
    let variable = btreemap! {
      0_usize  => VarPos { dis: array![0.9, 0.03, 0.04, 0.03], state: a },
      5_usize  => VarPos { dis: array![0.1, 0.6,  0.2,  0.1 ], state: c },
      10_usize => VarPos { dis: array![0.05, 0.05, 0.8, 0.1 ], state: g },
    };

    let seq_dis = MarginalSparseSeqDistribution {
      variable,
      variable_indel: btreemap! {},
      fixed: btreemap! {},
      fixed_counts: Composition::new(std::iter::empty::<AsciiChar>(), AsciiChar::from_byte_unchecked(b'-')),
      log_lh: 0.0,
    };

    let result = propagate_raw_per_site(&gtr, t, false, &seq_dis, None);

    for (&pos, var_pos) in &result.variable {
      let rate = gtr.site_rates.as_ref().unwrap()[pos];
      let exp_qt = gtr.expQt_with_rate(t, rate);
      let expected = exp_qt.dot(&seq_dis.variable[&pos].dis);
      assert_abs_diff_eq!(var_pos.dis, expected, epsilon = 1e-14);
    }
  }

  #[test]
  fn test_propagate_raw_per_site_backward() {
    use crate::gtr::get_gtr::{JC69Params, jc69};

    let a = AsciiChar::from_byte_unchecked(b'A');
    let c = AsciiChar::from_byte_unchecked(b'C');
    let g = AsciiChar::from_byte_unchecked(b'G');

    let mut gtr = jc69(JC69Params::default()).unwrap();
    gtr.set_site_rates(Array1::from_vec(vec![
      0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 0.1,
    ]));

    let t = 0.3;
    let variable = btreemap! {
      0_usize  => VarPos { dis: array![0.9, 0.03, 0.04, 0.03], state: a },
      5_usize  => VarPos { dis: array![0.1, 0.6,  0.2,  0.1 ], state: c },
      10_usize => VarPos { dis: array![0.05, 0.05, 0.8, 0.1 ], state: g },
    };

    let seq_dis = MarginalSparseSeqDistribution {
      variable,
      variable_indel: btreemap! {},
      fixed: btreemap! {},
      fixed_counts: Composition::new(std::iter::empty::<AsciiChar>(), AsciiChar::from_byte_unchecked(b'-')),
      log_lh: 0.0,
    };

    let result = propagate_raw_per_site(&gtr, t, true, &seq_dis, None);

    for (&pos, var_pos) in &result.variable {
      let rate = gtr.site_rates.as_ref().unwrap()[pos];
      let exp_qt = gtr.expQt_with_rate(t, rate);
      let expected = exp_qt.t().dot(&seq_dis.variable[&pos].dis);
      assert_abs_diff_eq!(var_pos.dis, expected, epsilon = 1e-14);
    }
  }

  #[test]
  fn test_combine_messages_keeps_deterministic_state_change_explicit() -> Result<(), Report> {
    use crate::alphabet::alphabet::Alphabet;
    use treetime_primitives::AsciiChar;

    let alphabet = Alphabet::default();
    let a = AsciiChar::from_byte_unchecked(b'A');
    let g = AsciiChar::from_byte_unchecked(b'G');
    let composition = Composition::with_sequence([g], alphabet.chars(), alphabet.gap());
    let fixed = alphabet
      .determined()
      .map(|state| Ok((state, alphabet.get_profile(state)?.clone())))
      .collect::<Result<BTreeMap<_, _>, Report>>()?;

    let message = MarginalSparseSeqDistribution {
      variable: btreemap! {},
      variable_indel: btreemap! {},
      fixed,
      fixed_counts: composition.clone(),
      log_lh: 0.0,
    };
    let variable_pos = btreemap! { 0_usize => g };
    let reference_states = vec![btreemap! { 0_usize => a }];

    let combined = combine_messages(
      &composition,
      &[message],
      &variable_pos,
      &reference_states,
      &alphabet,
      None,
    )?;

    let explicit = combined
      .variable
      .get(&0)
      .expect("state-changing deterministic sites must stay explicit");
    assert_eq!(explicit.state, g);
    assert_abs_diff_eq!(explicit.dis, alphabet.get_profile(a)?.clone(), epsilon = 1e-14);
    assert_eq!(combined.fixed_counts.get(g), Some(0));

    Ok(())
  }
  #[test]
  fn test_combine_messages_zero_multiplicity_does_not_create_nan_log_lh() -> Result<(), Report> {
    use crate::alphabet::alphabet::Alphabet;

    let alphabet = Alphabet::default();
    let a = AsciiChar::from_byte_unchecked(b'A');
    let c = AsciiChar::from_byte_unchecked(b'C');
    let g = AsciiChar::from_byte_unchecked(b'G');
    let t = AsciiChar::from_byte_unchecked(b'T');
    let composition = Composition::with_sequence([c], alphabet.chars(), alphabet.gap());
    let message = MarginalSparseSeqDistribution {
      variable: btreemap! {},
      variable_indel: btreemap! {},
      fixed: btreemap! {
        a => array![0.0, 0.0, 0.0, 0.0],
        c => array![0.0, 1.0, 0.0, 0.0],
        g => array![0.0, 0.0, 1.0, 0.0],
        t => array![0.0, 0.0, 0.0, 1.0],
      },
      fixed_counts: composition.clone(),
      log_lh: 0.0,
    };

    let combined = combine_messages(&composition, &[message], &btreemap! {}, &[], &alphabet, None)?;

    assert!(
      combined.log_lh.is_finite(),
      "zero-multiplicity fixed rows must not poison log-likelihood"
    );
    assert_eq!(combined.fixed_counts.get(a), Some(0));

    Ok(())
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
