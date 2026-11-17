use crate::alphabet::alphabet::Alphabet;
use crate::make_internal_report;
use crate::partition::marginal::posterior_states::most_likely_state;
use crate::partition::marginal::shared::normalize::divide_out_value;
use crate::partition::marginal::sparse::candidate_sites::CandidateSites;
use crate::partition::storage::sparse::{FixedProfiles, INLINE_STATES, SiteProfile, SiteScales, SparseSeqDistribution};
use crate::partition::storage::var_pos_map::VarPosMap;
use crate::seq::composition::Composition;
use eyre::Report;
use itertools::{Itertools, izip};
use ndarray::linalg::general_mat_vec_mul;
use ndarray::{Array1, Array2, ArrayViewMut1, aview_mut1, aview1};
use smallvec::{SmallVec, smallvec};
use treetime_primitives::{AsciiChar, LogLh, Seq};
use treetime_utils::array::ndarray::is_max_above;
use treetime_utils::array::softmax_with_log_norm::softmax_with_log_norm_inplace;
use treetime_utils::cast::ToF64Exact;

pub(crate) const EPS: f64 = 1e-4;

pub(crate) fn combine_messages(
  composition: &Composition,
  messages: &[&SparseSeqDistribution],
  sites: &CandidateSites,
  alphabet: &Alphabet,
  gtr_weight: Option<&Array1<f64>>,
  scale_inputs: Option<ScaleInputs<'_>>,
) -> Result<SparseSeqDistribution, Report> {
  let n_states = alphabet.n_canonical();
  let mut seq_dis = SparseSeqDistribution {
    variable: VarPosMap::with_capacity(sites.len()),
    fixed: FixedProfiles::default(),
    log_lh: messages.iter().map(|m| m.log_lh).sum(),
    scales: SiteScales::default(),
  };
  let mut candidate_scales = Vec::with_capacity(if scale_inputs.is_some() { sites.len() } else { 0 });

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
  let shared: SmallVec<[SharedProfile; INLINE_STATES]> = canonical
    .iter()
    .map(|&state| {
      let mut dis = initial_log.clone();
      for log_fixed in &log_fixed {
        add_log(&mut dis, log_fixed[state].dis.iter().copied());
      }
      let log_norm = softmax_with_log_norm_inplace(aview_mut1(&mut dis));
      SharedProfile { dis, log_norm }
    })
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
    if let Some(inputs) = &scale_inputs {
      let child_scales: f64 = izip!(messages, inputs.child_sequences)
        .map(|(msg, child_sequence)| message_site_scale(msg, child_sequence, pos, alphabet))
        .sum();
      candidate_scales.push((pos, log_norm + child_scales));
    }
    let dis = log_vec;
    seq_dis.log_lh += LogLh::new(log_norm);
    let index = canonical.iter().position(|&c| c == state);
    if let Some(index) = index {
      weights[index] -= 1.0;
    }

    let is_shared = is_site_resolved(&dis, EPS)
      && (all_states_equal || index.is_some_and(|index| is_shared_profile(&dis, &shared[index].dis, index)));
    if !is_shared {
      let label = if let Some(index) = index {
        kept_counts[index] = kept_counts[index].saturating_sub(1);
        state
      } else {
        alphabet.char(most_likely_state(&aview1(&dis)).unwrap_or(0))
      };
      seq_dis.variable.insert(pos, label, &dis);
    }
  }

  let shared_log_norms: SmallVec<[(AsciiChar, f64); INLINE_STATES]> = izip!(&canonical, &shared)
    .map(|(&state, profile)| (state, profile.log_norm))
    .collect();
  for (index, (&state, SharedProfile { dis, log_norm })) in izip!(&canonical, shared).enumerate() {
    seq_dis.log_lh += LogLh::new(weights[index] * log_norm);
    seq_dis.fixed.insert(state, kept_counts[index], dis);
  }
  if let Some(inputs) = scale_inputs {
    seq_dis.scales = combine_scales(
      messages,
      &inputs,
      &seq_dis,
      sites,
      &shared_log_norms,
      candidate_scales,
      alphabet,
    );
  }
  Ok(seq_dis)
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct ScaleInputs<'a> {
  pub(crate) sequence: &'a Seq,
  pub(crate) child_sequences: &'a [&'a Seq],
}

fn message_site_scale(message: &SparseSeqDistribution, sequence: &Seq, pos: usize, alphabet: &Alphabet) -> f64 {
  let state = sequence[pos];
  if message.variable.contains_key(pos) || !alphabet.is_canonical(state) {
    message.scales.site(pos)
  } else {
    message.scales.group(state) + message.scales.site(pos)
  }
}

fn combine_scales(
  messages: &[&SparseSeqDistribution],
  inputs: &ScaleInputs<'_>,
  seq_dis: &SparseSeqDistribution,
  sites: &CandidateSites,
  shared_log_norms: &[(AsciiChar, f64)],
  candidate_scales: Vec<(usize, f64)>,
  alphabet: &Alphabet,
) -> SiteScales {
  let groups: SmallVec<[(AsciiChar, f64); INLINE_STATES]> = shared_log_norms
    .iter()
    .map(|&(state, log_norm)| {
      let children: f64 = messages.iter().map(|msg| msg.scales.group(state)).sum();
      (state, children + log_norm)
    })
    .collect();
  let group = |state: AsciiChar| {
    groups
      .iter()
      .find(|(group_state, _)| *group_state == state)
      .map_or(0.0, |(_, scale)| *scale)
  };
  let own = candidate_scales.into_iter().map(|(pos, scale)| {
    let state = inputs.sequence[pos];
    let base = if seq_dis.variable.contains_key(pos) || !alphabet.is_canonical(state) {
      0.0
    } else {
      group(state)
    };
    (pos, scale - base)
  });
  let inherited = messages
    .iter()
    .flat_map(|msg| msg.scales.sites().iter().copied())
    .filter(|&(pos, _)| !sites.contains(pos));
  let sites = own
    .chain(inherited)
    .sorted_by_key(|&(pos, _)| pos)
    .coalesce(|(pos_a, a), (pos_b, b)| {
      if pos_a == pos_b {
        Ok((pos_a, a + b))
      } else {
        Err(((pos_a, a), (pos_b, b)))
      }
    })
    .collect();
  SiteScales::new(groups, sites)
}

pub(crate) fn propagate_raw(exp_qt: &Array2<f64>, seq_dis: &SparseSeqDistribution) -> SparseSeqDistribution {
  let n_states = exp_qt.nrows();
  let mut message = SparseSeqDistribution {
    variable: VarPosMap::with_capacity(seq_dis.variable.len()),
    fixed: FixedProfiles::default(),
    log_lh: seq_dis.log_lh,
    scales: seq_dis.scales.clone(),
  };
  for (pos, var) in &seq_dis.variable {
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

pub(crate) fn divide_out_site(profile: &[f64], message: &[f64]) -> SiteProfile {
  izip!(profile, message)
    .map(|(&value, &divisor)| divide_out_value(value, divisor))
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

fn is_shared_profile(dis: &[f64], shared: &[f64], state_index: usize) -> bool {
  dis == shared && most_likely_state(&aview1(dis)) == Some(state_index)
}

pub(crate) fn normalize_1d_inplace(mut dis: ArrayViewMut1<'_, f64>, weight: f64) -> f64 {
  let norm = dis.sum();
  if norm > 0.0 && norm.is_finite() {
    dis /= norm;
    weight * norm.ln()
  } else {
    let uniform = 1.0 / dis.len().to_f64_exact();
    dis.fill(uniform);
    f64::NEG_INFINITY
  }
}

struct SharedProfile {
  dis: SiteProfile,
  log_norm: f64,
}
