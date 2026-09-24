use crate::make_internal_report;
use crate::optimize::likelihood::evaluate_with_indels_log_lh_only;
use crate::partition::optimize::contribution::OptimizationContribution;
use eyre::Report;
use ndarray::Array1;
use ordered_float::OrderedFloat;

const GRID_SEARCH_POINTS: usize = 100;

pub(super) const GRID_SEARCH_MIN_UPPER: f64 = 0.5;

pub(super) fn min_branch_length_for_indels(indel_count: usize, one_mutation: f64) -> f64 {
  if indel_count > 0 { one_mutation * 0.01 } else { 0.0 }
}

pub(crate) fn is_zero_branch_optimal(contributions: &[OptimizationContribution]) -> bool {
  if !contributions
    .iter()
    .all(|contrib| contrib.has_unimodal_branch_likelihood())
  {
    return false;
  }

  if !contributions.iter().all(|contrib| contrib.all_sites_valid_at_zero()) {
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
  branch_length_for_extent: f64,
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  one_mutation: f64,
) -> Result<f64, Report> {
  let all_valid_at_zero = contributions.iter().all(|c| c.all_sites_valid_at_zero());

  let positive_might_lose_to_zero =
    candidate > 0.0 && is_zero_better_than_grid_best(contributions, indel_count, indel_rate, candidate)?;

  let zero_might_lose_to_positive = candidate == 0.0
    && indel_count == 0
    && all_valid_at_zero
    && !contributions.iter().all(|c| c.has_unimodal_branch_likelihood());

  let zero_invalid_passthrough = candidate == 0.0 && !all_valid_at_zero;

  if positive_might_lose_to_zero || zero_might_lose_to_positive || zero_invalid_passthrough {
    grid_search_inner(
      branch_length_for_extent,
      contributions,
      indel_count,
      indel_rate,
      one_mutation,
    )
  } else {
    Ok(candidate)
  }
}

pub(super) fn grid_search_inner(
  branch_length: f64,
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  one_mutation: f64,
) -> Result<f64, Report> {
  let branch_lengths = grid_search_branch_lengths(branch_length, one_mutation)?;

  let best_positive = branch_lengths
    .iter()
    .copied()
    .try_fold(
      None,
      |best, branch_length| -> Result<Option<(f64, OrderedFloat<f64>)>, Report> {
        let log_lh = OrderedFloat(
          evaluate_with_indels_log_lh_only(contributions, indel_count, indel_rate, branch_length)?.value(),
        );
        Ok(match best {
          Some((_, best_log_lh)) if best_log_lh >= log_lh => best,
          _ => Some((branch_length, log_lh)),
        })
      },
    )?
    .map(|(branch_length, _)| branch_length)
    .ok_or_else(|| {
      make_internal_report!("grid_search_inner: empty grid (GRID_SEARCH_POINTS = {GRID_SEARCH_POINTS})")
    })?;

  let zero_is_better = is_zero_better_than_grid_best(contributions, indel_count, indel_rate, best_positive)?;
  Ok(if zero_is_better { 0.0 } else { best_positive })
}

pub(super) fn grid_search_branch_lengths(branch_length: f64, one_mutation: f64) -> Result<Array1<f64>, Report> {
  let lower = 0.1 * one_mutation;
  let upper = f64::max(1.5 * branch_length + one_mutation, GRID_SEARCH_MIN_UPPER);
  Array1::geomspace(lower, upper, GRID_SEARCH_POINTS).ok_or_else(|| {
    make_internal_report!(
      "grid_search_branch_lengths: geomspace requires strictly positive same-sign bounds, got lower={lower}, upper={upper} (branch_length={branch_length}, one_mutation={one_mutation})"
    )
  })
}

pub(super) fn is_zero_better_than_grid_best(
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  best_positive: f64,
) -> Result<bool, Report> {
  if indel_count > 0 || !contributions.iter().all(|c| c.all_sites_valid_at_zero()) {
    return Ok(false);
  }
  let log_lh_zero = evaluate_with_indels_log_lh_only(contributions, indel_count, indel_rate, 0.0)?;
  let log_lh_best = evaluate_with_indels_log_lh_only(contributions, indel_count, indel_rate, best_positive)?;
  Ok(log_lh_zero > log_lh_best)
}
