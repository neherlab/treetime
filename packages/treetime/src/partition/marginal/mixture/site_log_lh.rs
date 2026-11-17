use crate::alphabet::alphabet::Alphabet;
use crate::partition::storage::sparse::{SparseEdgeBackward, SparseNodeState};
use crate::seq::composition::Composition;
use eyre::Report;
use itertools::{Itertools, izip};
use ndarray::{Array1, Array2, ArrayView1};
use smallvec::{SmallVec, smallvec};
use std::collections::BTreeMap;
use treetime_graph::key_map::EdgeMap;
use treetime_primitives::{AsciiChar, LogLh};
use treetime_utils::cast::ToF64Exact;
use treetime_utils::make_internal_error;

pub(crate) const INLINE_CATEGORIES: usize = 8;

pub(crate) type CategoryValues = SmallVec<[f64; INLINE_CATEGORIES]>;

#[derive(Clone, Debug)]
pub(crate) struct DenseMixture {
  pub(crate) site_log_lh: Array2<f64>,
  pub(crate) posterior: Array2<f64>,
  pub(crate) log_lh: LogLh,
}

pub(crate) fn dense_mixture(site_log_lh: &[Array1<f64>], weights: &Array1<f64>) -> Result<DenseMixture, Report> {
  let Some(length) = site_log_lh.first().map(Array1::len) else {
    return make_internal_error!("A rate mixture needs at least one category");
  };
  if site_log_lh.len() != weights.len() || site_log_lh.iter().any(|values| values.len() != length) {
    return make_internal_error!(
      "Rate mixture of {} categories got {} weights and per-site values of unequal length",
      site_log_lh.len(),
      weights.len()
    );
  }
  let site_log_lh = Array2::from_shape_fn((length, weights.len()), |(site, k)| site_log_lh[k][site]);
  let log_weights = weights.mapv(f64::ln);
  let mut posterior = Array2::zeros(site_log_lh.raw_dim());
  let mut log_lh = 0.0;
  for (values, mut row) in izip!(site_log_lh.rows(), posterior.rows_mut()) {
    let (site_total, site_posterior) = mix_site(&log_weights.view(), values.iter().copied());
    log_lh += site_total;
    row.assign(&ArrayView1::from(site_posterior.as_slice()));
  }
  Ok(DenseMixture {
    site_log_lh,
    posterior,
    log_lh: LogLh::new(log_lh),
  })
}

#[derive(Clone, Debug)]
pub(crate) struct SparseSiteLogLh {
  groups: SmallVec<[(AsciiChar, f64); 4]>,
  sites: Vec<(usize, AsciiChar, f64)>,
}

impl SparseSiteLogLh {
  pub(crate) fn from_pass(root: &SparseNodeState, backward: &EdgeMap<SparseEdgeBackward>, alphabet: &Alphabet) -> Self {
    let scales = &root.profile.scales;
    let sites = backward
      .values()
      .map(|edge| edge.msg_to_parent.scales.sites())
      .chain([scales.sites()])
      .map(|sites| sites.iter().copied())
      .kmerge_by(|a, b| a.0 < b.0)
      .coalesce(|(pos_a, a), (pos_b, b)| {
        if pos_a == pos_b {
          Ok((pos_a, a + b))
        } else {
          Err(((pos_a, a), (pos_b, b)))
        }
      })
      .map(|(pos, offset)| {
        let state = root.sequence[pos];
        let base = if root.profile.variable.contains_key(pos) || !alphabet.is_canonical(state) {
          0.0
        } else {
          scales.group(state)
        };
        (pos, state, offset + base)
      })
      .collect();
    let groups = alphabet.canonical().map(|state| (state, scales.group(state))).collect();
    Self { groups, sites }
  }

  fn group(&self, state: AsciiChar) -> f64 {
    self
      .groups
      .iter()
      .find(|(group_state, _)| *group_state == state)
      .map_or(0.0, |(_, value)| *value)
  }
}

#[derive(Clone, Debug)]
pub(crate) struct SparseMixture {
  pub(crate) log_lh: LogLh,
  prior: CategoryValues,
  groups: SmallVec<[MixedSite; 4]>,
  patterns: Vec<MixedSite>,
  sites: Vec<(usize, usize)>,
}

#[derive(Clone, Debug)]
pub(crate) struct MixedSite {
  pub(crate) state: AsciiChar,
  pub(crate) site_log_lh: CategoryValues,
  pub(crate) posterior: CategoryValues,
}

impl SparseMixture {
  pub(crate) fn new(
    categories: &[SparseSiteLogLh],
    weights: &Array1<f64>,
    composition: &Composition,
    alphabet: &Alphabet,
  ) -> Result<Self, Report> {
    if categories.len() != weights.len() || categories.is_empty() {
      return make_internal_error!(
        "Rate mixture of {} categories got {} weights",
        categories.len(),
        weights.len()
      );
    }
    let log_weights = weights.mapv(f64::ln);
    let mut listed: Vec<(usize, AsciiChar, SmallVec<[Option<f64>; INLINE_CATEGORIES]>)> = vec![];
    let entries = categories
      .iter()
      .enumerate()
      .map(|(k, category)| {
        category
          .sites
          .iter()
          .map(move |&(pos, state, value)| (pos, k, state, value))
      })
      .kmerge_by(|a, b| (a.0, a.1) < (b.0, b.1));
    for (pos, k, state, value) in entries {
      if listed.last().is_none_or(|&(last, ..)| last != pos) {
        listed.push((pos, state, smallvec![None; categories.len()]));
      }
      if let Some((_, _, values)) = listed.last_mut() {
        values[k] = Some(value);
      }
    }
    let mut log_lh = 0.0;
    let mut pattern_index = BTreeMap::<(AsciiChar, SmallVec<[u64; INLINE_CATEGORIES]>), usize>::new();
    let mut patterns = vec![];
    let mut pattern_totals = vec![];
    let sites = listed
      .into_iter()
      .map(|(pos, state, values)| {
        let site_log_lh: CategoryValues = izip!(categories, values)
          .map(|(category, value)| value.unwrap_or_else(|| category.group(state)))
          .collect();
        let key = (state, site_log_lh.iter().map(|value| value.to_bits()).collect());
        let pattern = *pattern_index.entry(key).or_insert_with(|| {
          let (total, posterior) = mix_site(&log_weights.view(), site_log_lh.iter().copied());
          pattern_totals.push(total);
          patterns.push(MixedSite {
            state,
            site_log_lh,
            posterior,
          });
          patterns.len() - 1
        });
        log_lh += pattern_totals[pattern];
        (pos, pattern)
      })
      .collect_vec();
    let groups = alphabet
      .canonical()
      .map(|state| {
        let listed_count = sites
          .iter()
          .filter(|&&(_, pattern)| patterns[pattern].state == state)
          .count();
        let count = composition.get(state).unwrap_or(0).saturating_sub(listed_count);
        let site_log_lh: CategoryValues = categories.iter().map(|category| category.group(state)).collect();
        let (total, posterior) = mix_site(&log_weights.view(), site_log_lh.iter().copied());
        if count > 0 {
          log_lh += total * count.to_f64_exact();
        }
        MixedSite {
          state,
          site_log_lh,
          posterior,
        }
      })
      .collect();
    Ok(Self {
      log_lh: LogLh::new(log_lh),
      prior: weights.iter().copied().collect(),
      groups,
      patterns,
      sites,
    })
  }

  pub(crate) fn listed(&self) -> impl Iterator<Item = (usize, usize, &MixedSite)> {
    self
      .sites
      .iter()
      .map(|&(pos, pattern)| (pos, pattern, &self.patterns[pattern]))
  }

  pub(crate) fn pattern_count(&self) -> usize {
    self.patterns.len()
  }

  pub(crate) fn pattern(&self, pattern: usize) -> &MixedSite {
    &self.patterns[pattern]
  }

  pub(crate) fn group(&self, state: AsciiChar) -> Option<&MixedSite> {
    self.groups.iter().find(|group| group.state == state)
  }

  pub(crate) fn site(&self, pos: usize, state: AsciiChar) -> Option<&MixedSite> {
    self
      .sites
      .binary_search_by_key(&pos, |(site, _)| *site)
      .ok()
      .map(|index| &self.patterns[self.sites[index].1])
      .or_else(|| self.group(state))
  }

  pub(crate) fn posterior(&self, pos: usize, state: AsciiChar) -> &[f64] {
    self.site(pos, state).map_or(&self.prior, |site| &site.posterior)
  }

  pub(crate) fn unlisted_posterior(&self, state: AsciiChar) -> &[f64] {
    self.group(state).map_or(&self.prior, |site| &site.posterior)
  }
}

pub(crate) fn mix_site(
  log_weights: &ArrayView1<'_, f64>,
  site_log_lh: impl Iterator<Item = f64>,
) -> (f64, CategoryValues) {
  let terms: CategoryValues = izip!(log_weights, site_log_lh).map(|(&w, l)| w + l).collect();
  let max = terms.iter().copied().fold(f64::NEG_INFINITY, f64::max);
  if !max.is_finite() {
    let n = terms.len().to_f64_exact();
    return (max, terms.iter().map(|_| 1.0 / n).collect());
  }
  let total = max + terms.iter().map(|term| (term - max).exp()).sum::<f64>().ln();
  (total, terms.iter().map(|term| (term - total).exp()).collect())
}
