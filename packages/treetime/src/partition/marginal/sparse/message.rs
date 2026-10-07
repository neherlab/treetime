use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::partition::storage::sparse::{SparseSeqDistribution, VarPos};
use crate::partition::storage::var_pos_map::VarPosMap;
use crate::seq::composition::Composition;
use eyre::Report;
use itertools::izip;
use maplit::btreemap;
use ndarray::{Array1, Array2};
use std::collections::BTreeMap;
use treetime_primitives::{AsciiChar, LogLh};
use treetime_utils::array::ndarray::is_max_above;
use treetime_utils::array::softmax_with_log_norm::softmax_with_log_norm_owned;
use treetime_utils::interval::range::range_contains;

const EPS: f64 = 1e-4;

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub(crate) fn combine_messages(
  composition: &Composition,
  messages: &[&SparseSeqDistribution],
  variable_pos: &BTreeMap<usize, AsciiChar>,
  reference_states: &[BTreeMap<usize, AsciiChar>],
  alphabet: &Alphabet,
  gtr_weight: Option<&Array1<f64>>,
) -> Result<SparseSeqDistribution, Report> {
  let mut seq_dis = SparseSeqDistribution {
    variable: VarPosMap::default(),
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
  let log_fixed: Vec<BTreeMap<AsciiChar, Array1<f64>>> = messages
    .iter()
    .map(|msg| {
      msg
        .fixed
        .iter()
        .map(|(&state, profile)| (state, profile.mapv(f64::ln)))
        .collect()
    })
    .collect();

  for (&pos, &state) in variable_pos {
    let mut all_states_equal = true;
    let mut log_vec = initial_log.clone();

    for (msg, states, log_fixed) in izip!(messages, reference_states, &log_fixed) {
      if let Some(var) = msg.variable.get(pos) {
        log_vec.zip_mut_with(&var.dis, |lv, &p| *lv += p.ln());
        if var.state != state {
          all_states_equal = false;
        }
      } else if let Some(ref_state) = states.get(&pos) {
        if alphabet.is_canonical(*ref_state) {
          log_vec.zip_mut_with(&log_fixed[ref_state], |lv, &lp| *lv += lp);
        }
        if ref_state != &state {
          all_states_equal = false;
        }
      } else {
        log_vec.zip_mut_with(&log_fixed[&state], |lv, &lp| *lv += lp);
      }
    }

    let (dis, log_norm) = softmax_with_log_norm_owned(log_vec);
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

    for log_fixed in &log_fixed {
      log_vec.zip_mut_with(&log_fixed[&state], |lv, &lp| *lv += lp);
    }

    let (dis, log_norm) = softmax_with_log_norm_owned(log_vec);
    seq_dis.log_lh += LogLh::new(fixed_counts[&state] * log_norm);
    seq_dis.fixed.insert(state, dis);
  }
  Ok(seq_dis)
}

pub(crate) fn propagate_raw(
  exp_qt: &Array2<f64>,
  seq_dis: &SparseSeqDistribution,
  transmission: Option<&[(usize, usize)]>,
) -> SparseSeqDistribution {
  let mut message = SparseSeqDistribution {
    variable: VarPosMap::default(),
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

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
pub(crate) fn propagate_raw_per_site(
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
    variable: VarPosMap::default(),
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

  for (&s, p) in &seq_dis.fixed {
    message.fixed.insert(s, default_exp_qt.dot(p));
  }

  message
}

fn is_site_resolved(dis: &Array1<f64>, epsilon: f64) -> bool {
  is_max_above(dis, 1.0 - epsilon)
}

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub(crate) fn normalize_1d_inplace(dis: &mut Array1<f64>, weight: f64) -> f64 {
  let norm = dis.sum();
  if norm > 0.0 && norm.is_finite() {
    *dis /= norm;
    weight * norm.ln()
  } else {
    dis.fill(1.0 / dis.len() as f64);
    f64::NEG_INFINITY
  }
}
