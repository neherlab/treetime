use crate::optimize::indel::IndelTerm;
use crate::optimize::likelihood::evaluate_with_indels_log_lh_only;
use crate::optimize::solver_outcome::{Convergence, EdgeSolution};
use crate::optimize::zero_boundary::search_upper_bound;
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::{make_internal_report, make_report};
use argmin::core::{CostFunction, Error, Executor, State, TerminationReason};
use argmin::solver::brent::BrentOpt;
use eyre::{Report, WrapErr};

const BRENT_MAX_ITER: u64 = 50;

pub(super) fn brent_sqrt_inner(
  branch_length: f64,
  contributions: &[OptimizationContribution],
  indels: IndelTerm,
  min_branch_length: f64,
  one_mutation: f64,
) -> Result<EdgeSolution, Report> {
  let bracket = brent_bracket(branch_length, min_branch_length, one_mutation);
  let cost_fn = SqrtBranchLengthCostFn { contributions, indels };
  brent_sqrt_minimize(&cost_fn, &bracket, BRENT_MAX_ITER)
    .wrap_err_with(|| format!("Brent branch-length optimization failed at branch_length={branch_length}"))
}

pub(super) fn brent_sqrt_minimize<F>(
  cost_fn: F,
  bracket: &BranchLengthBracket,
  max_iters: u64,
) -> Result<EdgeSolution, Report>
where
  F: CostFunction<Param = f64, Output = f64>,
{
  let solver = BrentOpt::new(bracket.lower.sqrt(), bracket.upper.sqrt());
  let res = Executor::new(cost_fn, solver)
    .configure(|cfg| cfg.max_iters(max_iters))
    .run()
    .map_err(|e| make_report!("argmin BrentOpt failed: {e}"))?;

  let state = res.state();
  let convergence = match state.get_termination_reason() {
    Some(TerminationReason::SolverConverged) => Convergence::Converged,
    Some(TerminationReason::MaxItersReached) => Convergence::IterationLimit,
    reason => {
      return Err(make_internal_report!(
        "argmin BrentOpt stopped for an unexpected reason: {reason:?}"
      ));
    },
  };
  let s = state
    .best_param
    .ok_or_else(|| make_internal_report!("argmin BrentOpt finished without a best parameter"))?;
  Ok(EdgeSolution {
    branch_length: s * s,
    convergence,
  })
}

pub(super) fn brent_bracket(branch_length: f64, min_branch_length: f64, one_mutation: f64) -> BranchLengthBracket {
  BranchLengthBracket {
    lower: min_branch_length.max(1e-12),
    upper: search_upper_bound(branch_length, one_mutation),
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
    let log_lh = evaluate_with_indels_log_lh_only(self.contributions, self.indels, t)
      .map_err(|error| Error::msg(error.to_string()))?;
    Ok(-log_lh)
  }
}

struct SqrtBranchLengthCostFn<'a> {
  contributions: &'a [OptimizationContribution],
  indels: IndelTerm,
}
