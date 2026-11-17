use crate::alphabet::alphabet::Alphabet;
use crate::make_internal_report;
use crate::partition::marginal::posterior_states::most_likely_state;
use crate::partition::marginal::shared::normalize::divide_out_value;
use crate::partition::marginal::sparse::candidate_sites::CandidateSites;
use crate::partition::storage::sparse::{FixedProfiles, INLINE_STATES, SiteProfile, SiteScales, SparseSeqDistribution};
use crate::partition::storage::var_pos_map::VarPosMap;
use crate::seq::composition::Composition;
use eyre::Report;
use itertools::izip;
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
  let mut scales =
    scale_inputs.map(|inputs| SiteScaleBuilder::new(inputs, messages, &canonical, &shared, sites.len()));

  for (pos, state, references) in sites.iter() {
    let mut all_states_equal = true;
    let mut log_vec = initial_log.clone();
    let mut child_variable = SmallVec::<[bool; 4]>::with_capacity(messages.len());

    for (msg, reference, log_fixed) in izip!(messages, references, &log_fixed) {
      let variable = msg.variable.get(pos);
      child_variable.push(variable.is_some());
      if let Some(var) = variable {
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
    if let Some(scales) = &mut scales {
      scales.push(messages, alphabet, pos, log_norm, &child_variable, !is_shared);
    }
  }

  for (index, (&state, SharedProfile { dis, log_norm })) in izip!(&canonical, shared).enumerate() {
    seq_dis.log_lh += LogLh::new(weights[index] * log_norm);
    seq_dis.fixed.insert(state, kept_counts[index], dis);
  }
  if let Some(scales) = scales {
    seq_dis.scales = scales.finish();
  }
  Ok(seq_dis)
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct ScaleInputs<'a> {
  pub(crate) sequence: &'a Seq,
  pub(crate) child_sequences: &'a [&'a Seq],
}

struct SiteScaleBuilder<'a> {
  inputs: ScaleInputs<'a>,
  groups: SmallVec<[(AsciiChar, f64); INLINE_STATES]>,
  sites: Vec<(usize, f64)>,
}

impl<'a> SiteScaleBuilder<'a> {
  fn new(
    inputs: ScaleInputs<'a>,
    messages: &[&SparseSeqDistribution],
    canonical: &[AsciiChar],
    shared: &[SharedProfile],
    capacity: usize,
  ) -> Self {
    let groups = izip!(canonical, shared)
      .map(|(&state, profile)| {
        let children: f64 = messages.iter().map(|msg| msg.scales.group(state)).sum();
        (state, children + profile.log_norm)
      })
      .collect();
    Self {
      inputs,
      groups,
      sites: Vec::with_capacity(capacity),
    }
  }

  fn push(
    &mut self,
    messages: &[&SparseSeqDistribution],
    alphabet: &Alphabet,
    pos: usize,
    log_norm: f64,
    child_variable: &[bool],
    own_variable: bool,
  ) {
    let child_scales: f64 = izip!(messages, self.inputs.child_sequences, child_variable)
      .map(|(msg, child_sequence, &variable)| {
        let state = child_sequence[pos];
        if variable || !alphabet.is_canonical(state) {
          0.0
        } else {
          msg.scales.group(state)
        }
      })
      .sum();
    let state = self.inputs.sequence[pos];
    let base = if own_variable || !alphabet.is_canonical(state) {
      0.0
    } else {
      self.group(state)
    };
    self.sites.push((pos, log_norm + child_scales - base));
  }

  fn group(&self, state: AsciiChar) -> f64 {
    self
      .groups
      .iter()
      .find(|(group_state, _)| *group_state == state)
      .map_or(0.0, |(_, scale)| *scale)
  }

  fn finish(self) -> SiteScales {
    SiteScales::new(self.groups, self.sites)
  }
}

pub(crate) fn propagate_raw(exp_qt: &Array2<f64>, seq_dis: &SparseSeqDistribution) -> SparseSeqDistribution {
  let n_states = exp_qt.nrows();
  let mut message = SparseSeqDistribution {
    variable: VarPosMap::with_capacity(seq_dis.variable.len()),
    fixed: FixedProfiles::default(),
    log_lh: seq_dis.log_lh,
    scales: seq_dis.scales.groups_only(),
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
