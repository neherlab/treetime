use crate::commands::optimize::optimize_unified::{
  GRID_SEARCH_MIN_UPPER, OptimizationContribution, evaluate_with_indels_log_lh_only,
};
use argmin::core::{CostFunction, Error, Executor};
use argmin::solver::brent::BrentOpt;

/// Maximum iterations for Brent's method.
pub(crate) const BRENT_MAX_ITER: u64 = 50;

/// Brent bracket endpoints in $t$-space shared by all three parameterizations.
///
/// The lower bound is `min_branch_length.max(1e-12)`: the `1e-12` floor avoids
/// $\ln(0)$ in the log-space cost function when the caller passes
/// `min_branch_length = 0` (no-indel case). The upper bound matches the grid
/// search range.
pub(crate) fn brent_bracket(branch_length: f64, min_branch_length: f64, one_mutation: f64) -> (f64, f64) {
  let lower = min_branch_length.max(1e-12);
  let upper = f64::max(1.5 * branch_length + one_mutation, GRID_SEARCH_MIN_UPPER);
  (lower, upper)
}

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
  let (lower, upper) = brent_bracket(branch_length, min_branch_length, one_mutation);

  let cost_fn = BranchLengthCostFn {
    contributions,
    indel_count,
    indel_rate,
    to_t: |t: f64| t,
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

/// Brent's method in $\sqrt{t}$ space.
///
/// Reparameterizes as $s = \sqrt{t}$: the cost function evaluates
/// $-\ell(s^2)$ and bracket endpoints transform as $\sqrt{\text{lower}}$,
/// $\sqrt{\text{upper}}$. The result is squared back to $t$-space.
///
/// This matches v0 exactly (`scipy.optimize.minimize_scalar(method='brent')`
/// with `t**2` inside `_neg_prob`). The $\sqrt{t}$ transform smooths the
/// objective near $t = 0$, giving Brent's parabolic interpolation a better
/// fit than raw $t$-space. Tolerance $\epsilon_s$ in $s$-space maps to
/// $t$-space precision $\approx 2 s^* \epsilon_s$, tighter near zero.
pub(crate) fn brent_sqrt_inner(
  branch_length: f64,
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  min_branch_length: f64,
  one_mutation: f64,
) -> f64 {
  let (lower, upper) = brent_bracket(branch_length, min_branch_length, one_mutation);

  let cost_fn = BranchLengthCostFn {
    contributions,
    indel_count,
    indel_rate,
    to_t: |s: f64| s * s,
  };

  let solver = BrentOpt::new(lower.sqrt(), upper.sqrt());
  let result = Executor::new(&cost_fn, solver)
    .configure(|cfg| cfg.max_iters(BRENT_MAX_ITER))
    .run();

  match result {
    Ok(res) => {
      let s = res.state().best_param.unwrap_or_else(|| branch_length.sqrt());
      s * s
    },
    Err(_) => branch_length,
  }
}

/// Brent's method in $\ln(t)$ space.
///
/// Reparameterizes as $u = \ln(t)$: the cost function evaluates
/// $-\ell(e^u)$ and bracket endpoints transform as $\ln(\text{lower})$,
/// $\ln(\text{upper})$. The result is exponentiated back to $t$-space.
///
/// The lower bound must be strictly positive to avoid $\ln(0)$; it is
/// clamped to $\max(\text{min\_branch\_length}, 10^{-12})$ (same as other
/// Brent variants). Tolerance $\epsilon_u$ in $u$-space is a natural
/// relative tolerance ($dt/t \approx du$).
pub(crate) fn brent_log_inner(
  branch_length: f64,
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  min_branch_length: f64,
  one_mutation: f64,
) -> f64 {
  let (lower, upper) = brent_bracket(branch_length, min_branch_length, one_mutation);

  let cost_fn = BranchLengthCostFn {
    contributions,
    indel_count,
    indel_rate,
    to_t: f64::exp,
  };

  let solver = BrentOpt::new(lower.ln(), upper.ln());
  let result = Executor::new(&cost_fn, solver)
    .configure(|cfg| cfg.max_iters(BRENT_MAX_ITER))
    .run();

  match result {
    Ok(res) => {
      let u = res.state().best_param.unwrap_or_else(|| branch_length.max(lower).ln());
      u.exp()
    },
    Err(_) => branch_length,
  }
}

/// Cost function for Brent optimization of branch length under a parameter
/// transform $t = \text{to\_t}(p)$.
///
/// Negates the combined (substitution + indel) log-likelihood because
/// `BrentOpt` minimizes, but we want to maximize. Specializing the transform:
/// $t$-space uses the identity, $\sqrt{t}$-space uses $p \mapsto p^2$,
/// $\ln(t)$-space uses $p \mapsto e^p$.
struct BranchLengthCostFn<'a, F: Fn(f64) -> f64> {
  contributions: &'a [OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  to_t: F,
}

impl<F: Fn(f64) -> f64> CostFunction for &BranchLengthCostFn<'_, F> {
  type Param = f64;
  type Output = f64;

  fn cost(&self, p: &Self::Param) -> Result<Self::Output, Error> {
    let t = (self.to_t)(*p);
    let log_lh = evaluate_with_indels_log_lh_only(self.contributions, self.indel_count, self.indel_rate, t);
    Ok(-log_lh)
  }
}
