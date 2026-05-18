use crate::make_internal_report;
use crate::optimize::likelihood::evaluate_with_indels_log_lh_only;
use crate::partition::optimization_contribution::OptimizationContribution;
use eyre::Report;
use ndarray::Array1;
use ordered_float::OrderedFloat;

/// Number of points in the grid search fallback.
const GRID_SEARCH_POINTS: usize = 100;

/// Minimum upper bound (subs/site) for the grid search fallback.
/// Ensures the grid covers the biologically plausible range regardless
/// of the current branch length estimate. Branches longer than 0.5
/// subs/site are rare in real phylogenies; the proportional upper
/// bound `1.5 * branch_length + one_mutation` extends beyond this
/// when the current estimate is already large.
pub(crate) const GRID_SEARCH_MIN_UPPER: f64 = 0.5;

/// Lower bound on the per-edge branch length used by Newton/Brent dispatch.
///
/// On indel-bearing edges the Poisson log-likelihood derivative diverges at
/// $t = 0$ (`-k/t^2`), so the optimizer must stay strictly positive. The
/// `0.01 * one_mutation` floor is well below the smallest meaningful branch
/// length ($1/L \sim 10^{-4}$ for $L = 10^4$) yet far enough from zero to
/// keep the indel Hessian finite. On no-indel edges zero is admissible, so
/// the floor is zero.
pub(crate) fn min_branch_length_for_indels(indel_count: usize, one_mutation: f64) -> f64 {
  if indel_count > 0 { one_mutation * 0.01 } else { 0.0 }
}

/// Log-spaced grid of candidate branch lengths for the grid search fallback.
///
/// The grid lower bound is `0.1 * one_mutation` (smallest meaningful branch length).
/// The upper bound is `max(1.5 * branch_length + one_mutation, GRID_SEARCH_MIN_UPPER)`,
/// ensuring coverage of the biologically plausible range (up to 0.5 subs/site) even
/// when the current branch length estimate is zero or very small.
///
/// Log spacing provides uniform resolution per decade across the 3-4 orders of
/// magnitude that branch lengths span in real phylogenies. A linear grid wastes
/// resolution near zero while under-sampling the long-branch tail.
///
/// Errors when `one_mutation <= 0` produces non-positive grid bounds (which would
/// make `geomspace` undefined).
pub(crate) fn grid_search_branch_lengths(branch_length: f64, one_mutation: f64) -> Result<Array1<f64>, Report> {
  let lower = 0.1 * one_mutation;
  let upper = f64::max(1.5 * branch_length + one_mutation, GRID_SEARCH_MIN_UPPER);
  Array1::geomspace(lower, upper, GRID_SEARCH_POINTS).ok_or_else(|| {
    make_internal_report!(
      "grid_search_branch_lengths: geomspace requires strictly positive same-sign bounds, got lower={lower}, upper={upper} (branch_length={branch_length}, one_mutation={one_mutation})"
    )
  })
}

/// Grid search fallback for non-concave regions.
pub(crate) fn grid_search_inner(
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
    .max_by_key(|&bl| {
      let log_lh = evaluate_with_indels_log_lh_only(contributions, indel_count, indel_rate, bl);
      OrderedFloat(log_lh)
    })
    .ok_or_else(|| {
      make_internal_report!("grid_search_inner: empty grid (GRID_SEARCH_POINTS = {GRID_SEARCH_POINTS})")
    })?;

  let zero_is_better = is_zero_better_than_grid_best(contributions, indel_count, indel_rate, best_positive);
  Ok(if zero_is_better { 0.0 } else { best_positive })
}

/// Whether zero branch length beats the best positive grid point.
///
/// For non-unimodal models that bypass the derivative shortcut, the grid
/// search must compare zero against the best positive candidate. When indels
/// are present ($k > 0$), the Poisson log-likelihood at $t = 0$ is $-\infty$,
/// so zero is never optimal. The `indel_count == 0` guard skips the comparison
/// entirely, matching the Newton path logic.
///
/// When `indel_count == 0`, both sides of the comparison use the indel-aware
/// evaluator for consistency with the grid evaluation, though the Poisson
/// contribution is zero in that case.
pub(crate) fn is_zero_better_than_grid_best(
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  best_positive: f64,
) -> bool {
  indel_count == 0 && contributions.iter().all(|c| c.all_sites_valid_at_zero()) && {
    let log_lh_zero = evaluate_with_indels_log_lh_only(contributions, indel_count, indel_rate, 0.0);
    let log_lh_best = evaluate_with_indels_log_lh_only(contributions, indel_count, indel_rate, best_positive);
    log_lh_zero > log_lh_best
  }
}

/// Check if zero branch length is optimal across all contributions.
///
/// The scientific criterion for zero-length optimality is the sign of the
/// log-likelihood derivative at t=0. For independent sites, the total
/// derivative is:
///
///   d/dt log L(0) = sum_i (sum_c k_{ic} lambda_c) / (sum_c k_{ic})
///
/// If this sum is negative, the log-likelihood decreases as t moves away
/// from zero, meaning zero is a local maximum.
///
/// This shortcut is restricted to models with proven unimodal branch-length
/// likelihood: JC69, F81, and binary models (Dinh & Matsen 2017, Corollary 3.1).
/// For these models, at most one stationary point exists on (0, infinity), so a
/// negative derivative at t=0 guarantees zero is the global maximum.
///
/// For models with multiple distinct nonzero eigenvalues (K80, HKY85, TN93,
/// general GTR), the likelihood can have two local maxima. A negative
/// derivative at t=0 only proves zero is a local maximum. The shortcut
/// returns false and the caller falls through to Newton/grid search.
///
/// The derivative-sign approach replaces an earlier implementation that
/// multiplied raw per-site likelihoods into a product and compared against
/// a fixed threshold (0.01). That approach suffered from two defects:
/// - Underflow: many small site likelihoods multiplied to f64 zero
/// - Scale dependence: the same local likelihood shape passed or failed
///   the threshold depending on alignment length
///
/// Before evaluating the derivative, each site's likelihood at t=0 must be
/// verified as positive and finite. If any site is degenerate (L_i(0) <= 0
/// or non-finite), the derivative formula divides by zero or produces
/// overflow. In that case, this function returns false and the caller falls
/// back to full Newton/grid optimization, which evaluates the likelihood
/// surface at positive branch lengths where the issue does not arise.
pub fn is_zero_branch_optimal(contributions: &[OptimizationContribution]) -> bool {
  // The derivative-sign shortcut is valid only when every partition uses a
  // model with proven unimodal branch-length likelihood. If any partition
  // uses a model that can be multimodal (K80, HKY85, TN93, general GTR),
  // skip the shortcut and let Newton/grid search evaluate the full surface.
  if !contributions
    .iter()
    .all(|contrib| contrib.has_unimodal_branch_likelihood())
  {
    return false;
  }

  // Verify every site's likelihood at t=0 is positive and finite.
  // If any site is degenerate, the derivative is undefined.
  if !contributions.iter().all(|contrib| contrib.all_sites_valid_at_zero()) {
    return false;
  }

  // Compute total derivative of log-likelihood at t=0.
  let derivative: f64 = contributions
    .iter()
    .map(|contrib| contrib.zero_branch_length_derivative())
    .sum();

  // Guard against overflow from near-canceling denominators. For general
  // GTR models with distinct eigenvalues, coefficients can nearly cancel
  // in L_i(0) while the numerator stays large, producing an infinite
  // derivative ratio. Finiteness of the sum catches this.
  if !derivative.is_finite() {
    return false;
  }

  derivative < 0.0
}

/// Reconcile a per-edge optimizer result against the boundary $t = 0$.
///
/// The pre-dispatch `is_zero_branch_optimal` shortcut catches the boundary
/// only for models with proven unimodal branch-length likelihood (JC69,
/// F81, binary; Dinh and Matsen 2017, Corollary 3.1). For every other
/// model the inner solver runs, and its result can be wrong in two
/// distinct ways:
///
/// 1. **Positive but worse than zero.** Methods that cannot evaluate at
///    $t = 0$ (NewtonLog, Brent variants) place a strictly positive
///    lower bound on the bracket. Their output can be a tiny positive
///    value like $10^{-12}$ when the true ML branch length is zero.
///    The cheap two-point check `is_zero_better_than_grid_best` flags
///    this case directly.
///
/// 2. **Exactly zero by step clamping.** `newton_inner` and
///    `newton_sqrt_inner` clamp the Newton step against `[-1.0, bl]`
///    (resp. `[sqrt_step_lower_bound(s), s]`); when the unclamped
///    step exceeds `bl`, the new branch length is exactly $0$.
///    `newton_sqrt_inner` in particular is observed to clamp to $0$
///    from $t_0 \approx 0.6$ on the Dinh and Matsen 2017 K80
///    $\kappa = 3$ counterexample
///    (`test_dispatch_zero_boundary_newton_sqrt_inner_can_clamp_to_zero_on_dinh_matsen_k80`),
///    even though that surface has a better positive mode near
///    $t \approx 0.2$. A `candidate > 0.0`-only gate would miss this.
///
/// Both failure modes route to `grid_search_inner` when at least one
/// partition is non-unimodal, all sites are valid at zero, and there are
/// no indels. The grid search enumerates 100 log-spaced positive
/// candidates and internally compares its best against zero, so it
/// returns $0$ when zero is genuinely the global maximum and the better
/// positive mode otherwise.
///
/// `branch_length_for_extent` is the input branch length used to size
/// the grid (the optimizer's output is unsuitable because it may be
/// clamped to $0$ or to a tiny floor like $10^{-12}$).
///
/// # Gate conditions
///
/// - `indel_count == 0`: the Poisson log-likelihood at $t = 0$ is
///   $-\infty$ when $k > 0$, so zero is never optimal in the presence
///   of indels.
/// - `all_sites_valid_at_zero`: zero is in the well-defined evaluation
///   domain. A degenerate site ($L_i(0) \leq 0$ or non-finite) makes
///   both sides of the comparison meaningless.
/// - Exact-zero case only: at least one partition is non-unimodal.
///   On unimodal models the pre-dispatch shortcut would have fired
///   (if the boundary is the global max) or the inner solver would
///   have found a positive mode (if it isn't). An exact-zero return
///   from the inner solver on an all-unimodal set of partitions is
///   either correct (and grid search would confirm) or unreachable,
///   so skipping the grid scan avoids wasted work in the common case.
///
/// 3. **Exactly zero on a partition where zero is invalid.** When any
///    site has $L_i(0) \leq 0$ or non-finite, the run_optimize_mixed
///    caller already bumps a zero input to `one_mutation` to keep the
///    optimizer in the well-defined domain. The optimizer can still
///    return $0$ via step clamping; that result re-enters the invalid
///    zero-evaluation domain and must be replaced. Route to grid search
///    on this case as well; `grid_search_inner` returns its best
///    positive candidate because `is_zero_better_than_grid_best` itself
///    gates on `all_sites_valid_at_zero`.
///
/// For Newton and NewtonSqrt on monotone-decreasing surfaces (e.g.
/// identical sequences) the helper is a no-op: the optimizer already
/// reached zero through its own grid-search fallback, and the
/// exact-zero re-verification returns zero.
pub(crate) fn reconcile_zero_boundary(
  candidate: f64,
  branch_length_for_extent: f64,
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  one_mutation: f64,
) -> Result<f64, Report> {
  let all_valid_at_zero = contributions.iter().all(|c| c.all_sites_valid_at_zero());

  let positive_might_lose_to_zero =
    candidate > 0.0 && is_zero_better_than_grid_best(contributions, indel_count, indel_rate, candidate);

  let zero_might_lose_to_positive = candidate == 0.0
    && indel_count == 0
    && all_valid_at_zero
    && !contributions.iter().all(|c| c.has_unimodal_branch_likelihood());

  // Optimizer returned 0.0 on a partition where zero is in the undefined
  // evaluation domain (some site has L_i(0) <= 0 or non-finite). The caller
  // bumped the input away from zero for exactly this reason; preserving the
  // optimizer's clamped 0.0 would re-enter the invalid domain on the next
  // marginal update.
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
