use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::representation::payload::sparse::{MarginalSparseSeqDistribution, VarPos};
use crate::seq::composition::Composition;
use eyre::Report;
use maplit::btreemap;
use ndarray::{Array1, Array2};
use std::collections::BTreeMap;
use std::iter::zip;
use treetime_primitives::AsciiChar;
use treetime_utils::array::ndarray::is_max_above;
use treetime_utils::array::softmax_with_log_norm::softmax_with_log_norm;
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

    let (dis, log_norm) = softmax_with_log_norm(log_vec.view());
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

    let (dis, log_norm) = softmax_with_log_norm(log_vec.view());
    seq_dis.log_lh += fixed_counts[&state] * log_norm;
    seq_dis.fixed.insert(state, dis);
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
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::seq::composition::Composition;
  use approx::assert_abs_diff_eq;
  use ndarray::array;

  #[test]
  fn test_propagate_raw_per_site_forward() {
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

    // JC69 closed-form: P(t)_ii = (1/4)(1 + 3*exp(-4/3 * mu * rate * t)),
    // P(t)_ij = (1/4)(1 - exp(-4/3 * mu * rate * t)). mu=0.75 after GTR normalization.
    let expected_0 = array![0.8094601846762877, 0.06064424518648728, 0.06925132495073787, 0.06064424518648728];
    let expected_5 = array![0.16767825458589602, 0.4420840726329092, 0.22255941819529867, 0.16767825458589602];
    let expected_10 = array![0.05591089329029837, 0.05591089329029837, 0.7837450434516795, 0.10443316996772378];

    assert_abs_diff_eq!(result.variable[&0].dis, expected_0, epsilon = 1e-14);
    assert_abs_diff_eq!(result.variable[&5].dis, expected_5, epsilon = 1e-14);
    assert_abs_diff_eq!(result.variable[&10].dis, expected_10, epsilon = 1e-14);
  }

  #[test]
  fn test_propagate_raw_per_site_backward() {
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

    // JC69 P(t) is symmetric, so P(t)^T = P(t). Same expected values as forward.
    let expected_0 = array![0.8094601846762877, 0.06064424518648728, 0.06925132495073787, 0.06064424518648728];
    let expected_5 = array![0.16767825458589602, 0.4420840726329092, 0.22255941819529867, 0.16767825458589602];
    let expected_10 = array![0.05591089329029837, 0.05591089329029837, 0.7837450434516795, 0.10443316996772378];

    assert_abs_diff_eq!(result.variable[&0].dis, expected_0, epsilon = 1e-14);
    assert_abs_diff_eq!(result.variable[&5].dis, expected_5, epsilon = 1e-14);
    assert_abs_diff_eq!(result.variable[&10].dis, expected_10, epsilon = 1e-14);
  }
}
