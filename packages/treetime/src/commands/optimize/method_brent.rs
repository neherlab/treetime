use crate::commands::optimize::optimize_unified::{
  GRID_SEARCH_MIN_UPPER, OptimizationContribution, evaluate_with_indels_log_lh_only,
};
use argmin::core::{CostFunction, Error, Executor};
use argmin::solver::brent::BrentOpt;

/// Maximum iterations for Brent's method.
pub(crate) const BRENT_MAX_ITER: u64 = 50;

/// Brent's method inner loop (derivative-free, bracket-based).
///
/// Wraps the combined (substitution + indel) log-likelihood in an
/// `argmin::CostFunction` and runs `BrentOpt` to find the minimum of
/// the negated log-likelihood within a bracket.
///
/// The bracket lower bound is `min_branch_length` (or a small epsilon when
/// no indels are present), ensuring coverage of optima near zero. The upper
/// bound matches the grid search range. For concave objectives (JC69/F81
/// substitution + Poisson indel), the bracket is guaranteed to contain the
/// unique maximum: at the lower bound the Poisson derivative dominates
/// (positive, pushing right), and at the upper bound the substitution
/// derivative dominates (negative, pushing left).
pub(crate) fn brent_inner(
  branch_length: f64,
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  min_branch_length: f64,
  one_mutation: f64,
) -> f64 {
  // Lower bound: smallest t where the objective is well-defined.
  // When indels are present, min_branch_length > 0 (prevents Poisson
  // singularity at t=0). When no indels, use a small epsilon.
  let lower = min_branch_length.max(1e-12);
  let upper = f64::max(1.5 * branch_length + one_mutation, GRID_SEARCH_MIN_UPPER);

  let cost_fn = BranchLengthCostFunction {
    contributions,
    indel_count,
    indel_rate,
  };

  let solver = BrentOpt::new(lower, upper);
  let result = Executor::new(&cost_fn, solver)
    .configure(|cfg| cfg.max_iters(BRENT_MAX_ITER))
    .run();

  match result {
    Ok(res) => res.state().best_param.unwrap_or(branch_length),
    Err(_) => branch_length,
  }
}

/// Cost function for Brent optimization of branch length.
///
/// Negates the combined (substitution + indel) log-likelihood because
/// `BrentOpt` minimizes, but we want to maximize.
struct BranchLengthCostFunction<'a> {
  contributions: &'a [OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
}

impl CostFunction for &BranchLengthCostFunction<'_> {
  type Param = f64;
  type Output = f64;

  fn cost(&self, t: &Self::Param) -> Result<Self::Output, Error> {
    let log_lh = evaluate_with_indels_log_lh_only(self.contributions, self.indel_count, self.indel_rate, *t);
    Ok(-log_lh)
  }
}
