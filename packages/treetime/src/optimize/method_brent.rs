use crate::optimize::likelihood::evaluate_with_indels_log_lh_only;
use crate::optimize::zero_boundary::GRID_SEARCH_MIN_UPPER;
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::{make_internal_report, make_report};
use argmin::core::{CostFunction, Error, Executor};
use argmin::solver::brent::BrentOpt;
use eyre::Report;

const BRENT_MAX_ITER: u64 = 50;

pub(super) fn brent_sqrt_inner(
  branch_length: f64,
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  min_branch_length: f64,
  one_mutation: f64,
) -> Result<f64, Report> {
  let BranchLengthBracket { lower, upper } = brent_bracket(branch_length, min_branch_length, one_mutation);

  let cost_fn = SqrtBranchLengthCostFn {
    contributions,
    indel_count,
    indel_rate,
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

impl CostFunction for &SqrtBranchLengthCostFn<'_> {
  type Param = f64;
  type Output = f64;

  fn cost(&self, s: &Self::Param) -> Result<Self::Output, Error> {
    let t = s * s;
    let log_lh = evaluate_with_indels_log_lh_only(self.contributions, self.indel_count, self.indel_rate, t)
      .map_err(|error| Error::msg(error.to_string()))?;
    Ok(-log_lh)
  }
}

struct SqrtBranchLengthCostFn<'a> {
  contributions: &'a [OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
}
