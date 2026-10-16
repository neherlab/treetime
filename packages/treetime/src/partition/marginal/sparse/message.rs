use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::make_internal_report;
use crate::partition::marginal::sparse::candidate_sites::CandidateSites;
use crate::partition::storage::sparse::{FixedProfiles, INLINE_STATES, SiteProfile, SparseSeqDistribution};
use crate::partition::storage::var_pos_map::VarPosMap;
use crate::seq::composition::Composition;
use eyre::Report;
use itertools::izip;
use ndarray::linalg::general_mat_vec_mul;
use ndarray::{Array1, Array2, ArrayViewMut1, aview_mut1, aview1};
use smallvec::{SmallVec, smallvec};
use treetime_primitives::{AsciiChar, LogLh};
use treetime_utils::array::ndarray::is_max_above;
use treetime_utils::array::softmax_with_log_norm::softmax_with_log_norm_inplace;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::interval::range::range_contains;

const EPS: f64 = 1e-4;

pub(crate) fn combine_messages(
  composition: &Composition,
  messages: &[&SparseSeqDistribution],
  sites: &CandidateSites,
  alphabet: &Alphabet,
  gtr_weight: Option<&Array1<f64>>,
) -> Result<SparseSeqDistribution, Report> {
  let n_states = alphabet.n_canonical();
  let mut seq_dis = SparseSeqDistribution {
    variable: VarPosMap::with_capacity(sites.len()),
    fixed: FixedProfiles::default(),
    log_lh: messages.iter().map(|m| m.log_lh).sum(),
  };

  let canonical: SmallVec<[AsciiChar; INLINE_STATES]> = alphabet.canonical().collect();
  let mut kept_counts = canonical
    .iter()
    .map(|&state| {
      composition
        .get(state)
        .ok_or_else(|| make_internal_report!("Canonical state {state} has no count"))
    })
    .collect::<Result<SmallVec<[usize; INLINE_STATES]>, Report>>()?;
  let mut weights: SmallVec<[f64; INLINE_STATES]> = kept_counts.iter().map(|count| count.to_f64_exact()).collect();

  let initial_log: SiteProfile =
    gtr_weight.map_or_else(|| smallvec![0.0; n_states], |w| w.iter().map(|p| p.ln()).collect());
  let log_fixed: SmallVec<[FixedProfiles; 2]> = messages
    .iter()
    .map(|msg| msg.fixed.map_dis(|profile| profile.iter().map(|p| p.ln()).collect()))
    .collect();

  for (pos, state, references) in sites.iter() {
    let mut all_states_equal = true;
    let mut log_vec = initial_log.clone();

    for (msg, reference, log_fixed) in izip!(messages, references, &log_fixed) {
      if let Some(var) = msg.variable.get(pos) {
        add_log(&mut log_vec, var.dis.iter().map(|p| p.ln()));
        if var.state != state {
          all_states_equal = false;
        }
      } else if let Some(ref_state) = *reference {
        if alphabet.is_canonical(ref_state) {
          add_log(&mut log_vec, log_fixed[ref_state].dis.iter().copied());
        }
        if ref_state != state {
          all_states_equal = false;
        }
      } else {
        add_log(&mut log_vec, log_fixed[state].dis.iter().copied());
      }
    }

    let log_norm = softmax_with_log_norm_inplace(aview_mut1(&mut log_vec));
    let dis = log_vec;
    seq_dis.log_lh += LogLh::new(log_norm);
    let index = canonical.iter().position(|&c| c == state);
    if let Some(index) = index {
      weights[index] -= 1.0;
    }

    if !is_site_resolved(&dis, EPS) || !all_states_equal {
      if let Some(index) = index {
        kept_counts[index] = kept_counts[index].saturating_sub(1);
      }
      seq_dis.variable.insert(pos, state, &dis);
    }
  }

  for (index, &state) in canonical.iter().enumerate() {
    let mut log_vec = initial_log.clone();

    for log_fixed in &log_fixed {
      add_log(&mut log_vec, log_fixed[state].dis.iter().copied());
    }

    let log_norm = softmax_with_log_norm_inplace(aview_mut1(&mut log_vec));
    seq_dis.log_lh += LogLh::new(weights[index] * log_norm);
    seq_dis.fixed.insert(state, kept_counts[index], log_vec);
  }
  Ok(seq_dis)
}

pub(crate) fn propagate_raw(
  exp_qt: &Array2<f64>,
  seq_dis: &SparseSeqDistribution,
  transmission: Option<&[(usize, usize)]>,
) -> SparseSeqDistribution {
  let n_states = exp_qt.nrows();
  let mut message = SparseSeqDistribution {
    variable: VarPosMap::with_capacity(seq_dis.variable.len()),
    fixed: FixedProfiles::default(),
    log_lh: seq_dis.log_lh,
  };
  for (pos, var) in &seq_dis.variable {
    if let Some(transmission) = &transmission {
      if !range_contains(transmission, *pos) {
        continue;
      }
    }
    message
      .variable
      .insert_with(*pos, var.state, n_states, |row| propagate_into(exp_qt, var.dis, row));
  }

  for fixed in seq_dis.fixed.iter() {
    message
      .fixed
      .insert(fixed.state, fixed.count, propagate_site(exp_qt, &fixed.dis));
  }

  message
}

pub(crate) fn propagate_raw_per_site(
  gtr: &GTR,
  site_rates: &Array1<f64>,
  branch_length: f64,
  transpose: bool,
  seq_dis: &SparseSeqDistribution,
  transmission: Option<&[(usize, usize)]>,
) -> SparseSeqDistribution {
  let default_exp_qt = if transpose {
    gtr.expQt(branch_length).t().to_owned()
  } else {
    gtr.expQt(branch_length)
  };

  let mut message = SparseSeqDistribution {
    variable: VarPosMap::with_capacity(seq_dis.variable.len()),
    fixed: FixedProfiles::default(),
    log_lh: seq_dis.log_lh,
  };

  for (&pos, var) in &seq_dis.variable {
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
    message.variable.insert_with(pos, var.state, exp_qt_pos.nrows(), |row| {
      propagate_into(&exp_qt_pos, var.dis, row);
    });
  }

  for fixed in seq_dis.fixed.iter() {
    message
      .fixed
      .insert(fixed.state, fixed.count, propagate_site(&default_exp_qt, &fixed.dis));
  }

  message
}

pub(crate) fn divide_out_site(profile: &[f64], message: &[f64]) -> SiteProfile {
  izip!(profile, message)
    .map(|(&value, &divisor)| value / divisor.max(f64::MIN_POSITIVE))
    .collect()
}

fn propagate_site(exp_qt: &Array2<f64>, profile: &[f64]) -> SiteProfile {
  let mut result = smallvec![0.0; exp_qt.nrows()];
  propagate_into(exp_qt, profile, &mut result);
  result
}

fn propagate_into(exp_qt: &Array2<f64>, profile: &[f64], out: &mut [f64]) {
  general_mat_vec_mul(1.0, exp_qt, &aview1(profile), 0.0, &mut aview_mut1(out));
}

fn add_log(log_vec: &mut SiteProfile, log_values: impl Iterator<Item = f64>) {
  for (lv, lp) in izip!(log_vec.iter_mut(), log_values) {
    *lv += lp;
  }
}

fn is_site_resolved(dis: &[f64], epsilon: f64) -> bool {
  is_max_above(&aview1(dis), 1.0 - epsilon)
}

pub(crate) fn normalize_1d_inplace(mut dis: ArrayViewMut1<'_, f64>, weight: f64) -> f64 {
  let norm = dis.sum();
  if norm > 0.0 && norm.is_finite() {
    dis /= norm;
    weight * norm.ln()
  } else {
    dis.fill(1.0 / dis.len().to_f64_exact());
    f64::NEG_INFINITY
  }
}
