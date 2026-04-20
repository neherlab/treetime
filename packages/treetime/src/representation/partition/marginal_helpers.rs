use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::representation::payload::sparse::{MarginalSparseSeqDistribution, VarPos};
use crate::seq::composition::Composition;
use eyre::Report;
use maplit::btreemap;
use ndarray::{Array1, Array2, ArrayView1};
use std::collections::BTreeMap;
use std::iter::zip;
use treetime_primitives::AsciiChar;
use treetime_utils::array::ndarray::{is_max_above, max_or};
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

    if !is_site_resolved(&dis, EPS) || !all_states_equal {
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
    seq_dis.log_lh += fixed_counts[&state] * log_norm;
    seq_dis.fixed.insert(state, dis);
  }
  Ok(seq_dis)
}

/// Fused LSE + softmax: returns `(softmax(x), LSE(x))` from unnormalized log-probabilities.
///
/// Naive `exp()` of log-probabilities overflows or underflows for large magnitudes.
/// The LSE shift (`exp(x - max)`) keeps all exponents in a safe range. LSE and
/// softmax share the intermediate `exp(x - max)` in a single fused pass. Direct
/// division `exp(x - max) / sum` preserves tighter sum-to-one than the decomposed
/// `exp(x - LSE(x))` form, which would introduce an extra `ln` to `exp` round-trip.
///
/// Degenerate input (all -inf): returns uniform distribution, LSE = -inf.
pub fn logsumexp_normalize(log_vec: ArrayView1<'_, f64>) -> (Array1<f64>, f64) {
  // LSE shift constant: subtract max before exp to prevent overflow
  let max_val = max_or(&log_vec, f64::NEG_INFINITY);

  // All states have zero probability: fall back to uniform
  if !max_val.is_finite() {
    let n = log_vec.len();
    return (Array1::from_elem(n, 1.0 / n as f64), f64::NEG_INFINITY);
  }

  // Shared intermediate: exp(x - max) used by both LSE and softmax
  let shifted = log_vec.mapv(|v| (v - max_val).exp());
  let shifted_sum = shifted.sum();

  // LSE(x) = max + ln(sum(exp(x - max)))
  let log_norm = max_val + shifted_sum.ln();

  // softmax(x) = exp(x - max) / sum(exp(x - max))
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

/// Whether a posterior is dominated by a single state (peak probability >= 1 - epsilon),
/// meaning the site can be demoted from variable to fixed in the sparse representation.
fn is_site_resolved(dis: &Array1<f64>, epsilon: f64) -> bool {
  is_max_above(dis, 1.0 - epsilon)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::seq::composition::Composition;
  use approx::{assert_abs_diff_eq, assert_ulps_eq};
  use ndarray::array;
  use rstest::rstest;

  const NEG_INF: f64 = f64::NEG_INFINITY;

  /// All-finite inputs: softmax output has no zeros.
  /// Expected values computed from the mathematical definition of LSE + softmax.
  #[rustfmt::skip]
  #[rstest]
  #[case::uniform_4(             Array1::from_elem(4, 0.25_f64.ln()))]
  #[case::descending(            array![0.0, -1.0, -2.0, -3.0])]
  #[case::ascending(             array![-3.0, -2.0, -1.0, 0.0])]
  #[case::all_same_nonzero(      array![3.0, 3.0, 3.0, 3.0])]
  #[case::close_values(          array![-0.001, -0.002, -0.003, -0.004])]
  #[case::wide_spread(           array![0.0, -10.0, -20.0, -30.0])]
  #[case::single_state(          array![5.0])]
  #[case::two_equal(             array![0.0, 0.0])]
  #[case::two_asymmetric(        array![0.0, -5.0])]
  #[case::three_ascending(       array![1.0, 2.0, 3.0])]
  #[case::twenty_uniform(        Array1::from_elem(20, 0.0))]
  #[case::nuc_realistic_profile( array![-0.1, -2.3, -4.5, -6.7])]
  #[trace]
  fn test_logsumexp_normalize_finite(#[case] input: Array1<f64>) {
    let (normalized, log_norm) = logsumexp_normalize(input.view());

    let (expected_normalized, expected_log_norm) = helpers::reference_lse_softmax(&input);

    helpers::assert_valid_distribution(&normalized);
    assert_ulps_eq!(normalized, expected_normalized, max_ulps = 2);
    assert_ulps_eq!(log_norm, expected_log_norm, max_ulps = 2);
  }

  /// Inputs containing -inf: some output probabilities are exactly 0.0.
  /// ULPs undefined at 0.0, so abs_diff for the probability vector, ULPs for log_norm.
  #[rustfmt::skip]
  #[rstest]
  #[case::dominant_rest_tiny(      array![0.0, -100.0, -100.0, -100.0])]
  #[case::mixed_finite_neg_inf(    array![0.0, NEG_INF, -1.0, NEG_INF])]
  #[case::single_finite_rest_inf(  array![NEG_INF, -3.0, NEG_INF, NEG_INF])]
  #[case::first_finite_rest_inf(   array![0.0, NEG_INF, NEG_INF, NEG_INF])]
  #[case::last_finite_rest_inf(    array![NEG_INF, NEG_INF, NEG_INF, 7.0])]
  #[case::two_finite_two_inf(      array![-1.0, NEG_INF, -2.0, NEG_INF])]
  #[trace]
  fn test_logsumexp_normalize_with_neg_inf(#[case] input: Array1<f64>) {
    let (normalized, log_norm) = logsumexp_normalize(input.view());

    let (expected_normalized, expected_log_norm) = helpers::reference_lse_softmax(&input);

    helpers::assert_valid_distribution(&normalized);
    assert_abs_diff_eq!(normalized, expected_normalized, epsilon = 1e-15);
    assert_ulps_eq!(log_norm, expected_log_norm, max_ulps = 2);
  }

  /// All-inf degenerate input: uniform fallback, log_norm = -inf.
  #[rustfmt::skip]
  #[rstest]
  #[case::n1(  1)]
  #[case::n2(  2)]
  #[case::n3(  3)]
  #[case::n4(  4)]
  #[case::n20( 20)]
  #[trace]
  fn test_logsumexp_normalize_degenerate(#[case] n_states: usize) {
    let log_vec = Array1::from_elem(n_states, NEG_INF);
    let (normalized, log_norm) = logsumexp_normalize(log_vec.view());

    let expected_uniform = Array1::from_elem(n_states, 1.0 / n_states as f64);
    assert_ulps_eq!(normalized, expected_uniform, max_ulps = 0);
    assert!(log_norm == NEG_INF, "expected NEG_INFINITY, got {log_norm}");
  }

  /// Shift invariance: adding a constant to all inputs preserves softmax,
  /// shifts log_norm by the same constant.
  #[rustfmt::skip]
  #[rstest]
  #[case::shift_neg_1000( -1000.0)]
  #[case::shift_neg_100(   -100.0)]
  #[case::shift_pos_100(    100.0)]
  #[case::shift_pos_1000(  1000.0)]
  #[trace]
  fn test_logsumexp_normalize_shift_invariant(#[case] offset: f64) {
    let base = array![0.0, -1.0, -2.0, -3.0];
    let shifted = base.mapv(|v| v + offset);

    let (base_norm, base_log_norm) = logsumexp_normalize(base.view());
    let (shifted_norm, shifted_log_norm) = logsumexp_normalize(shifted.view());

    helpers::assert_valid_distribution(&shifted_norm);
    assert_ulps_eq!(shifted_norm, base_norm, max_ulps = 2);
    assert_ulps_eq!(shifted_log_norm, base_log_norm + offset, max_ulps = 2);
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
      assert_abs_diff_eq!(normalized.sum(), 1.0, epsilon = 1e-15);
    }

    /// Reference LSE + softmax from mathematical definition (decomposed form).
    /// Uses the same max-shift for LSE stability but computes softmax via
    /// `exp(x - LSE(x))` rather than the fused `exp(x - max) / sum`.
    pub fn reference_lse_softmax(input: &Array1<f64>) -> (Array1<f64>, f64) {
      let max_val = input.iter().copied().fold(f64::NEG_INFINITY, f64::max);
      if !max_val.is_finite() {
        return (Array1::from_elem(input.len(), 1.0 / input.len() as f64), f64::NEG_INFINITY);
      }
      let lse = max_val + input.iter().map(|&v| (v - max_val).exp()).sum::<f64>().ln();
      let softmax = input.mapv(|v| (v - lse).exp());
      (softmax, lse)
    }
  }
}
