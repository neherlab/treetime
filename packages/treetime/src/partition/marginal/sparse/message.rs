#[cfg(test)]
mod __tests__;

use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::partition::storage::sparse::{SparseSeqDistribution, VarPos};
use crate::seq::composition::Composition;
use eyre::Report;
use maplit::btreemap;
use ndarray::{Array1, Array2};
use std::collections::{BTreeMap, BTreeSet};
use std::iter::zip;
use treetime_primitives::{AsciiChar, LogLh};
use treetime_utils::array::ndarray::is_max_above;
use treetime_utils::array::softmax_with_log_norm::softmax_with_log_norm;
use treetime_utils::interval::range::range_contains;

pub const EPS: f64 = 1e-4;

pub fn combine_messages(
  composition: &Composition,
  messages: &[SparseSeqDistribution],
  variable_pos: &BTreeMap<usize, AsciiChar>,
  reference_states: &[BTreeMap<usize, AsciiChar>],
  alphabet: &Alphabet,
  gtr_weight: Option<&Array1<f64>>,
) -> Result<SparseSeqDistribution, Report> {
  let mut seq_dis = SparseSeqDistribution {
    variable: btreemap! {},
    variable_indel: BTreeSet::new(),
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
    seq_dis.log_lh += LogLh::new(log_norm);
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
    seq_dis.log_lh += LogLh::new(fixed_counts[&state] * log_norm);
    seq_dis.fixed.insert(state, dis);
  }
  Ok(seq_dis)
}

pub fn propagate_raw(
  exp_qt: &Array2<f64>,
  seq_dis: &SparseSeqDistribution,
  transmission: Option<&[(usize, usize)]>,
) -> SparseSeqDistribution {
  let mut message = SparseSeqDistribution {
    variable: btreemap! {},
    variable_indel: BTreeSet::new(),
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
  seq_dis: &SparseSeqDistribution,
  transmission: Option<&[(usize, usize)]>,
) -> SparseSeqDistribution {
  let site_rates = gtr
    .site_rates
    .as_ref()
    .expect("propagate_raw_per_site requires site_rates");
  let default_exp_qt = if transpose {
    gtr.expQt(branch_length).t().to_owned()
  } else {
    gtr.expQt(branch_length)
  };

  let mut message = SparseSeqDistribution {
    variable: btreemap! {},
    variable_indel: BTreeSet::new(),
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

/// Normalize a 1D sparse-site distribution in place.
///
/// Normalizes `dis` to sum to 1 and returns the log-likelihood contribution
/// `weight * ln(norm)`, where `weight` is the number of sites sharing this
/// distribution (`1.0` for a single variable site, the fixed-state count for a
/// fixed block). When the norm is non-positive or non-finite, falls back to a
/// uniform distribution and contributes `f64::NEG_INFINITY` (unweighted,
/// matching the dense 2D normalization fallback).
pub fn normalize_1d_inplace(dis: &mut Array1<f64>, weight: f64) -> f64 {
  let norm = dis.sum();
  if norm > 0.0 && norm.is_finite() {
    *dis /= norm;
    weight * norm.ln()
  } else {
    dis.fill(1.0 / dis.len() as f64);
    f64::NEG_INFINITY
  }
}
