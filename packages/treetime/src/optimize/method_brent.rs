use crate::optimize::likelihood::evaluate_with_indels_log_lh_only;
use crate::optimize::zero_boundary::GRID_SEARCH_MIN_UPPER;
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::{make_internal_report, make_report};
use argmin::core::{CostFunction, Error, Executor};
use argmin::solver::brent::BrentOpt;
use eyre::Report;

const BRENT_MAX_ITER: u64 = 50;

pub(super) fn brent_inner(
  branch_length: f64,
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  min_branch_length: f64,
  one_mutation: f64,
) -> Result<f64, Report> {
  let (lower, upper) = brent_bracket(branch_length, min_branch_length, one_mutation);

  let cost_fn = BranchLengthCostFn {
    contributions,
    indel_count,
    indel_rate,
    to_t: |t: f64| t,
  };

  let solver = BrentOpt::new(lower, upper);
  let res = Executor::new(&cost_fn, solver)
    .configure(|cfg| cfg.max_iters(BRENT_MAX_ITER))
    .run()
    .map_err(|e| make_report!("brent_inner: argmin BrentOpt failed at branch_length={branch_length}: {e}"))?;

  res.state().best_param.ok_or_else(|| {
    make_internal_report!("brent_inner: solver succeeded but reported no best_param at branch_length={branch_length}")
  })
}

pub(super) fn brent_sqrt_inner(
  branch_length: f64,
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  min_branch_length: f64,
  one_mutation: f64,
) -> Result<f64, Report> {
  let (lower, upper) = brent_bracket(branch_length, min_branch_length, one_mutation);

  let cost_fn = BranchLengthCostFn {
    contributions,
    indel_count,
    indel_rate,
    to_t: |s: f64| s * s,
  };

  let solver = BrentOpt::new(lower.sqrt(), upper.sqrt());
  let res = Executor::new(&cost_fn, solver)
    .configure(|cfg| cfg.max_iters(BRENT_MAX_ITER))
    .run()
    .map_err(|e| make_report!("brent_sqrt_inner: argmin BrentOpt failed at branch_length={branch_length}: {e}"))?;

  let s = res.state().best_param.ok_or_else(|| {
    make_internal_report!(
      "brent_sqrt_inner: solver succeeded but reported no best_param at branch_length={branch_length}"
    )
  })?;
  Ok(s * s)
}

pub(super) fn brent_log_inner(
  branch_length: f64,
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  min_branch_length: f64,
  one_mutation: f64,
) -> Result<f64, Report> {
  let (lower, upper) = brent_bracket(branch_length, min_branch_length, one_mutation);

  let cost_fn = BranchLengthCostFn {
    contributions,
    indel_count,
    indel_rate,
    to_t: f64::exp,
  };

  let solver = BrentOpt::new(lower.ln(), upper.ln());
  let res = Executor::new(&cost_fn, solver)
    .configure(|cfg| cfg.max_iters(BRENT_MAX_ITER))
    .run()
    .map_err(|e| make_report!("brent_log_inner: argmin BrentOpt failed at branch_length={branch_length}: {e}"))?;

  let u = res.state().best_param.ok_or_else(|| {
    make_internal_report!(
      "brent_log_inner: solver succeeded but reported no best_param at branch_length={branch_length}"
    )
  })?;
  Ok(u.exp())
}

pub(super) fn brent_bracket(branch_length: f64, min_branch_length: f64, one_mutation: f64) -> (f64, f64) {
  let lower = min_branch_length.max(1e-12);
  let upper = f64::max(1.5 * branch_length + one_mutation, GRID_SEARCH_MIN_UPPER);
  (lower, upper)
}

impl<F: Fn(f64) -> f64> CostFunction for &BranchLengthCostFn<'_, F> {
  type Param = f64;
  type Output = f64;

  fn cost(&self, p: &Self::Param) -> Result<Self::Output, Error> {
    let t = (self.to_t)(*p);
    let log_lh = evaluate_with_indels_log_lh_only(self.contributions, self.indel_count, self.indel_rate, t)
      .map_err(|error| Error::msg(error.to_string()))?;
    Ok(-log_lh)
  }
}

struct BranchLengthCostFn<'a, F: Fn(f64) -> f64> {
  contributions: &'a [OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  to_t: F,
}
