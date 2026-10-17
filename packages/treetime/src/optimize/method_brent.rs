use crate::brent::minimize_brent;
use crate::optimize::indel::IndelTerm;
use crate::optimize::likelihood::evaluate_with_indels_log_lh_only;
use crate::optimize::zero_boundary::GRID_SEARCH_MIN_UPPER;
use crate::partition::optimize::contribution::OptimizationContribution;
use eyre::{Report, WrapErr};

const BRENT_MAX_ITER: u64 = 50;

pub(super) fn brent_sqrt_inner(
  branch_length: f64,
  contributions: &[OptimizationContribution],
  indels: IndelTerm,
  min_branch_length: f64,
  one_mutation: f64,
) -> Result<f64, Report> {
  let BranchLengthBracket { lower, upper } = brent_bracket(branch_length, min_branch_length, one_mutation);
  let minimum = minimize_brent(lower.sqrt(), upper.sqrt(), None, BRENT_MAX_ITER, |s| {
    Ok(-evaluate_with_indels_log_lh_only(contributions, indels, s * s)?.value())
  })
  .wrap_err_with(|| {
    format!("When optimizing the branch length {branch_length} by Brent's method on its square root")
  })?;
  Ok(minimum.x * minimum.x)
}

pub(super) fn brent_bracket(branch_length: f64, min_branch_length: f64, one_mutation: f64) -> BranchLengthBracket {
  BranchLengthBracket {
    lower: min_branch_length.max(1e-12),
    upper: f64::max(1.5 * branch_length + one_mutation, GRID_SEARCH_MIN_UPPER),
  }
}

pub(super) struct BranchLengthBracket {
  pub lower: f64,
  pub upper: f64,
}
