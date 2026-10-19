use crate::make_internal_report;
use crate::optimize::indel::IndelTerm;
use crate::optimize::likelihood::evaluate_with_indels_log_lh_only;
use crate::partition::optimize::contribution::OptimizationContribution;
use eyre::Report;
use ndarray::Array1;
use ordered_float::OrderedFloat;

const GRID_SEARCH_POINTS: usize = 100;

const MIN_SEARCH_UPPER_BOUND: f64 = 0.5;

pub(super) fn start_branch_length(
  branch_length: f64,
  indels: IndelTerm,
  zero_length: ZeroLength,
  one_mutation: f64,
) -> f64 {
  if branch_length > 0.0 {
    branch_length
  } else if indels.requires_positive_length() {
    indels
      .length_estimate()
      .map_or(one_mutation, |length| length.max(one_mutation))
  } else if zero_length == ZeroLength::ExcludedBySites {
    one_mutation
  } else {
    branch_length
  }
}

pub(super) fn min_branch_length_for_indels(indels: IndelTerm, one_mutation: f64) -> f64 {
  if indels.requires_positive_length() {
    one_mutation * 0.01
  } else {
    0.0
  }
}

pub(crate) fn is_zero_branch_optimal(contributions: &[OptimizationContribution], zero_length: ZeroLength) -> bool {
  if zero_length != ZeroLength::Admissible || shared_certified_decay_rate(contributions).is_none() {
    return false;
  }

  let derivative: f64 = contributions
    .iter()
    .map(|contrib| contrib.zero_branch_length_derivative())
    .sum();

  if !derivative.is_finite() {
    return false;
  }

  derivative < 0.0
}

pub(super) fn reconcile_zero_boundary(
  candidate: f64,
  start_branch_length: f64,
  contributions: &[OptimizationContribution],
  indels: IndelTerm,
  zero_length: ZeroLength,
  one_mutation: f64,
) -> Result<f64, Report> {
  let zero_beats_candidate = candidate > 0.0 && is_zero_better_than(contributions, indels, zero_length, candidate)?;

  let zero_candidate_uncertified =
    candidate == 0.0 && zero_length == ZeroLength::Admissible && shared_certified_decay_rate(contributions).is_none();

  let zero_candidate_excluded = candidate == 0.0 && zero_length == ZeroLength::ExcludedBySites;

  if zero_beats_candidate || zero_candidate_uncertified || zero_candidate_excluded {
    grid_search_branch_length(start_branch_length, contributions, indels, zero_length, one_mutation)
  } else {
    Ok(candidate)
  }
}

pub(super) fn grid_search_branch_length(
  branch_length: f64,
  contributions: &[OptimizationContribution],
  indels: IndelTerm,
  zero_length: ZeroLength,
  one_mutation: f64,
) -> Result<f64, Report> {
  let branch_lengths = grid_search_branch_lengths(branch_length, one_mutation)?;

  let best_positive = branch_lengths
    .iter()
    .copied()
    .try_fold(
      None,
      |best, branch_length| -> Result<Option<(f64, OrderedFloat<f64>)>, Report> {
        let log_lh = OrderedFloat(evaluate_with_indels_log_lh_only(contributions, indels, branch_length)?.value());
        Ok(match best {
          Some((_, best_log_lh)) if best_log_lh >= log_lh => best,
          _ => Some((branch_length, log_lh)),
        })
      },
    )?
    .map(|(branch_length, _)| branch_length)
    .ok_or_else(|| {
      make_internal_report!("grid_search_branch_length: empty grid (GRID_SEARCH_POINTS = {GRID_SEARCH_POINTS})")
    })?;

  let zero_is_better = is_zero_better_than(contributions, indels, zero_length, best_positive)?;
  Ok(if zero_is_better { 0.0 } else { best_positive })
}

pub(super) fn grid_search_branch_lengths(branch_length: f64, one_mutation: f64) -> Result<Array1<f64>, Report> {
  let lower = 0.1 * one_mutation;
  let upper = search_upper_bound(branch_length, one_mutation);
  Array1::geomspace(lower, upper, GRID_SEARCH_POINTS).ok_or_else(|| {
    make_internal_report!(
      "grid_search_branch_lengths: geomspace requires strictly positive same-sign bounds, got lower={lower}, upper={upper} (branch_length={branch_length}, one_mutation={one_mutation})"
    )
  })
}

pub(super) fn search_upper_bound(branch_length: f64, one_mutation: f64) -> f64 {
  f64::max(1.5 * branch_length + one_mutation, MIN_SEARCH_UPPER_BOUND)
}

pub(super) fn is_zero_better_than(
  contributions: &[OptimizationContribution],
  indels: IndelTerm,
  zero_length: ZeroLength,
  branch_length: f64,
) -> Result<bool, Report> {
  if zero_length != ZeroLength::Admissible {
    return Ok(false);
  }
  let log_lh_zero = evaluate_with_indels_log_lh_only(contributions, indels, 0.0)?;
  let log_lh_other = evaluate_with_indels_log_lh_only(contributions, indels, branch_length)?;
  Ok(log_lh_zero > log_lh_other)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ZeroLength {
  Admissible,
  ExcludedBySites,
  ExcludedByIndels,
}

impl ZeroLength {
  pub(crate) fn of_edge(contributions: &[OptimizationContribution], indels: IndelTerm) -> Self {
    if !contributions
      .iter()
      .all(OptimizationContribution::all_sites_valid_at_zero)
    {
      Self::ExcludedBySites
    } else if indels.requires_positive_length() {
      Self::ExcludedByIndels
    } else {
      Self::Admissible
    }
  }
}

fn shared_certified_decay_rate(contributions: &[OptimizationContribution]) -> Option<f64> {
  let (first, rest) = contributions.split_first()?;
  let decay_rate = first.certified_decay_rate()?;
  rest
    .iter()
    .all(|contrib| contrib.certified_decay_rate() == Some(decay_rate))
    .then_some(decay_rate)
}
