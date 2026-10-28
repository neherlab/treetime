use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::partition::marginal::posterior_states::most_likely_state;
use crate::partition::marginal::shared::normalize::divide_out_value;
use crate::partition::marginal::sparse::candidate_sites::CandidateSites;
use crate::partition::marginal::sparse::site_grouping::{GroupedSites, ProfileSource};
use crate::partition::storage::site_groups::SiteGroups;
use crate::partition::storage::sparse::{FixedProfiles, INLINE_STATES, SiteProfile, SparseSeqDistribution};
use crate::partition::storage::var_pos_map::VarPosMap;
use crate::seq::composition::Composition;
use crate::{make_internal_error, make_internal_report};
use eyre::Report;
use itertools::{Itertools, izip};
use ndarray::linalg::general_mat_vec_mul;
use ndarray::{Array1, Array2, ArrayViewMut1, aview_mut1, aview1};
use smallvec::{SmallVec, smallvec};
use treetime_primitives::{AsciiChar, LogLh};
use treetime_utils::array::ndarray::is_max_above;
use treetime_utils::array::softmax_with_log_norm::softmax_with_log_norm_inplace;
use treetime_utils::cast::ToF64Exact;

const EPS: f64 = 1e-4;

pub(crate) fn combine_messages(
  composition: &Composition,
  messages: &[&SparseSeqDistribution],
  sites: &CandidateSites,
  groups: Vec<GroupedSites>,
  alphabet: &Alphabet,
  gtr_weight: Option<&Array1<f64>>,
) -> Result<SparseSeqDistribution, Report> {
  let mut product = MessageProduct::new(composition, messages, sites.len(), alphabet, gtr_weight)?;
  for (pos, state, references) in sites.iter() {
    product.add_site(pos, state, references);
  }
  let group_profiles = groups
    .into_iter()
    .map(|group| product.add_group(group))
    .filter_ok(|(_, profiles)| !profiles.is_empty())
    .collect::<Result<Vec<_>, Report>>()?;
  product.seq_dis.groups = SiteGroups::from_groups(group_profiles);
  Ok(product.finish())
}

struct MessageProduct<'a> {
  messages: &'a [&'a SparseSeqDistribution],
  alphabet: &'a Alphabet,
  canonical: SmallVec<[AsciiChar; INLINE_STATES]>,
  initial_log: SiteProfile,
  log_fixed: SmallVec<[FixedProfiles; 2]>,
  kept_counts: SmallVec<[usize; INLINE_STATES]>,
  weights: SmallVec<[f64; INLINE_STATES]>,
  seq_dis: SparseSeqDistribution,
}

impl<'a> MessageProduct<'a> {
  fn new(
    composition: &Composition,
    messages: &'a [&'a SparseSeqDistribution],
    n_sites: usize,
    alphabet: &'a Alphabet,
    gtr_weight: Option<&Array1<f64>>,
  ) -> Result<Self, Report> {
    let canonical: SmallVec<[AsciiChar; INLINE_STATES]> = alphabet.canonical().collect();
    let kept_counts = canonical
      .iter()
      .map(|&state| {
        composition
          .get(state)
          .ok_or_else(|| make_internal_report!("Canonical state {state} has no count"))
      })
      .collect::<Result<SmallVec<[usize; INLINE_STATES]>, Report>>()?;
    let weights = kept_counts.iter().map(|count| count.to_f64_exact()).collect();
    let initial_log = gtr_weight.map_or_else(
      || smallvec![0.0; alphabet.n_canonical()],
      |w| w.iter().map(|p| p.ln()).collect(),
    );
    let log_fixed = messages
      .iter()
      .map(|msg| msg.fixed.map_dis(|profile| profile.iter().map(|p| p.ln()).collect()))
      .collect();
    let seq_dis = SparseSeqDistribution {
      variable: VarPosMap::with_capacity(n_sites),
      groups: SiteGroups::default(),
      fixed: FixedProfiles::default(),
      log_lh: messages.iter().map(|m| m.log_lh).sum(),
    };
    Ok(Self {
      messages,
      alphabet,
      canonical,
      initial_log,
      log_fixed,
      kept_counts,
      weights,
      seq_dis,
    })
  }

  fn add_site(&mut self, pos: usize, state: AsciiChar, references: &[Option<AsciiChar>]) {
    let alphabet = self.alphabet;
    let mut all_states_equal = true;
    let mut log_vec = self.initial_log.clone();
    for (msg, reference, log_fixed) in izip!(self.messages, references, &self.log_fixed) {
      if let Some(var) = msg.variable.get(pos) {
        add_log(&mut log_vec, var.dis.iter().map(|p| p.ln()));
        if var.state != state {
          all_states_equal = false;
        }
      } else if let Some(ref_state) = *reference {
        if alphabet.is_canonical(ref_state) {
          add_shared_log(&mut log_vec, msg, log_fixed, pos, ref_state);
        }
        if ref_state != state {
          all_states_equal = false;
        }
      } else {
        add_shared_log(&mut log_vec, msg, log_fixed, pos, state);
      }
    }

    let log_norm = softmax_with_log_norm_inplace(aview_mut1(&mut log_vec));
    let dis = log_vec;
    self.seq_dis.log_lh += LogLh::new(log_norm);
    let index = self.canonical.iter().position(|&c| c == state);
    if let Some(index) = index {
      self.weights[index] -= 1.0;
    }

    if !is_site_resolved(&dis, EPS) || !all_states_equal {
      let label = if let Some(index) = index {
        self.kept_counts[index] = self.kept_counts[index].saturating_sub(1);
        state
      } else {
        alphabet.char(most_likely_state(&aview1(&dis)).unwrap_or(0))
      };
      self.seq_dis.variable.insert(pos, label, &dis);
    }
  }

  fn add_group(
    &mut self,
    GroupedSites {
      sources,
      ranges,
      counts,
    }: GroupedSites,
  ) -> Result<(Vec<(usize, usize)>, FixedProfiles), Report> {
    let mut profiles = FixedProfiles::default();
    for (state, count) in counts {
      let mut log_vec = self.initial_log.clone();
      let mut shared = true;
      let mut all_states_equal = true;
      for (msg, source, log_fixed) in izip!(self.messages, &sources, &self.log_fixed) {
        match source {
          ProfileSource::Missing => {
            shared = false;
            all_states_equal = false;
          },
          ProfileSource::Group(group) => {
            if let Some(entry) = msg.groups.profiles(*group).get(state) {
              shared = false;
              add_log(&mut log_vec, entry.dis.iter().map(|p| p.ln()));
            } else {
              add_log(&mut log_vec, log_fixed[state].dis.iter().copied());
            }
          },
          ProfileSource::Fixed => add_log(&mut log_vec, log_fixed[state].dis.iter().copied()),
        }
      }
      if shared {
        continue;
      }
      let index = self
        .canonical
        .iter()
        .position(|&c| c == state)
        .ok_or_else(|| make_internal_report!("Grouped state {state} is not canonical"))?;
      let log_norm = softmax_with_log_norm_inplace(aview_mut1(&mut log_vec));
      let multiplicity = count.to_f64_exact();
      self.seq_dis.log_lh += LogLh::new(multiplicity * log_norm);
      self.weights[index] -= multiplicity;
      if !is_site_resolved(&log_vec, EPS) || !all_states_equal {
        self.kept_counts[index] = self.kept_counts[index].saturating_sub(count);
        profiles.insert(state, count, log_vec);
      }
    }
    Ok((ranges, profiles))
  }

  fn finish(mut self) -> SparseSeqDistribution {
    for (index, &state) in self.canonical.iter().enumerate() {
      let mut log_vec = self.initial_log.clone();
      for log_fixed in &self.log_fixed {
        add_log(&mut log_vec, log_fixed[state].dis.iter().copied());
      }
      let log_norm = softmax_with_log_norm_inplace(aview_mut1(&mut log_vec));
      self.seq_dis.log_lh += LogLh::new(self.weights[index] * log_norm);
      self.seq_dis.fixed.insert(state, self.kept_counts[index], log_vec);
    }
    self.seq_dis
  }
}

pub(crate) fn propagate_raw(exp_qt: &Array2<f64>, seq_dis: &SparseSeqDistribution) -> SparseSeqDistribution {
  let n_states = exp_qt.nrows();
  let mut message = SparseSeqDistribution {
    variable: VarPosMap::with_capacity(seq_dis.variable.len()),
    groups: seq_dis.groups.map_profiles(|dis| propagate_site(exp_qt, dis)),
    fixed: FixedProfiles::default(),
    log_lh: seq_dis.log_lh,
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

pub(crate) fn propagate_raw_per_site(
  gtr: &GTR,
  site_rates: &Array1<f64>,
  branch_length: f64,
  transpose: bool,
  seq_dis: &SparseSeqDistribution,
) -> Result<SparseSeqDistribution, Report> {
  if !seq_dis.groups.is_empty() {
    return make_internal_error!("Grouped sites cannot be propagated with per-site rates");
  }
  let default_exp_qt = if transpose {
    gtr.expQt(branch_length).t().to_owned()
  } else {
    gtr.expQt(branch_length)
  };

  let mut message = SparseSeqDistribution {
    variable: VarPosMap::with_capacity(seq_dis.variable.len()),
    groups: SiteGroups::default(),
    fixed: FixedProfiles::default(),
    log_lh: seq_dis.log_lh,
  };

  for (&pos, var) in &seq_dis.variable {
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

  Ok(message)
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

fn add_shared_log(
  log_vec: &mut SiteProfile,
  message: &SparseSeqDistribution,
  log_fixed: &FixedProfiles,
  pos: usize,
  state: AsciiChar,
) {
  if let Some(entry) = message.groups.entry(pos, state) {
    add_log(log_vec, entry.dis.iter().map(|p| p.ln()));
  } else {
    add_log(log_vec, log_fixed[state].dis.iter().copied());
  }
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
    let uniform = 1.0 / dis.len().to_f64_exact();
    dis.fill(uniform);
    f64::NEG_INFINITY
  }
}
