use crate::commands::optimize::args::BranchOptMethod;
use crate::commands::optimize::method_brent::{brent_inner, brent_log_inner, brent_sqrt_inner};
use crate::commands::optimize::method_newton::{
  NEWTON_ABS_TOL, NEWTON_REL_TOL, newton_inner, newton_log_inner, newton_sqrt_inner,
};
use crate::commands::optimize::optimize_dense;
use crate::commands::optimize::optimize_dense_eval::{evaluate_dense_contribution, evaluate_dense_contribution_impl};
use crate::commands::optimize::optimize_indel::{estimate_indel_rate, poisson_indel_log_lh};
use crate::commands::optimize::optimize_sparse;
use crate::commands::optimize::optimize_sparse_eval::{
  evaluate_sparse_contribution, evaluate_sparse_contribution_impl,
};
use crate::commands::optimize::partition_ops::PartitionOptimizeOps;
use crate::gtr::gtr::GTR;
use crate::representation::partition::marginal_dense::PartitionMarginalDense;
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::{make_error, make_internal_report};
use eyre::Report;
use itertools::Either;
use ndarray::{Array1, ArrayView1};
use ordered_float::OrderedFloat;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;

/// Number of points in the grid search fallback.
const GRID_SEARCH_POINTS: usize = 100;

/// Minimum upper bound (subs/site) for the grid search fallback.
/// Ensures the grid covers the biologically plausible range regardless
/// of the current branch length estimate. Branches longer than 0.5
/// subs/site are rare in real phylogenies; the proportional upper
/// bound `1.5 * branch_length + one_mutation` extends beyond this
/// when the current estimate is already large.
pub(crate) const GRID_SEARCH_MIN_UPPER: f64 = 0.5;

/// Newton convergence tolerance in $t$-space.
///
/// A step $|\Delta t| \leq \eta t$ corresponds to relative tolerance $\eta$
/// in branch-length space. The absolute floor prevents degeneration at $t = 0$.
pub(crate) fn newton_tolerance_t(t: f64) -> f64 {
  f64::max(NEWTON_REL_TOL * t, NEWTON_ABS_TOL)
}

/// Newton convergence tolerance in $\sqrt{t}$-space.
///
/// Since $\Delta t / t \approx 2 \Delta s / s$ for small steps, matching the
/// target relative tolerance $\eta$ in branch-length space requires
/// $|\Delta s| \leq \frac{\eta}{2} s$. The factor 0.5 corrects the
/// coordinate-space scaling.
pub(crate) fn newton_tolerance_sqrt(s: f64) -> f64 {
  f64::max(0.5 * NEWTON_REL_TOL * s, NEWTON_ABS_TOL)
}

/// Newton convergence tolerance in $\ln(t)$-space.
///
/// Since $\Delta t / t \approx \Delta u$ for small steps, a constant tolerance
/// $\eta$ in $u$-space directly corresponds to relative tolerance $\eta$ in
/// branch-length space, regardless of $|u|$. No scaling by the current value.
pub(crate) fn newton_tolerance_log() -> f64 {
  // ln(1 + η) ≈ η for small η; use the exact form for correctness
  NEWTON_REL_TOL.ln_1p()
}

/// Metrics computed during branch length optimization
#[derive(Clone, Debug, Default)]
pub struct OptimizationMetrics {
  /// Log likelihood value
  pub log_lh: f64,
  /// First derivative (gradient) of log likelihood with respect to branch length
  pub derivative: f64,
  /// Second derivative (hessian) of log likelihood with respect to branch length
  pub second_derivative: f64,
}

impl OptimizationMetrics {
  pub fn new(log_lh: f64, derivative: f64, second_derivative: f64) -> Self {
    Self {
      log_lh,
      derivative,
      second_derivative,
    }
  }

  /// Add another set of metrics to this one
  pub fn add(&mut self, other: &OptimizationMetrics) {
    self.log_lh += other.log_lh;
    self.derivative += other.derivative;
    self.second_derivative += other.second_derivative;
  }
}

#[allow(clippy::large_enum_variant)]
pub enum OptimizationContribution {
  Dense(optimize_dense::PartitionContribution),
  Sparse(optimize_sparse::PartitionContribution),
}

impl OptimizationContribution {
  /// Create optimization contribution from a dense partition
  ///
  /// Extracts the message data from the partition edge and computes the coefficient
  /// matrix using the dense optimization approach.
  pub fn from_dense(edge_key: GraphEdgeKey, partition: &PartitionMarginalDense) -> Self {
    let edge_partition = &partition.edges[&edge_key];
    let contribution = optimize_dense::get_coefficients(
      &edge_partition.msg_to_parent,
      &edge_partition.msg_to_child,
      &partition.gtr,
    );
    OptimizationContribution::Dense(contribution)
  }

  /// Create optimization contribution from a sparse partition
  ///
  /// Extracts variable positions and mutations from the sparse partition and
  /// computes site contributions for optimization.
  pub fn from_sparse(edge_key: GraphEdgeKey, partition: &PartitionMarginalSparse) -> Result<Self, Report> {
    let contribution = optimize_sparse::get_coefficients(edge_key, partition)?;
    Ok(OptimizationContribution::Sparse(contribution))
  }

  /// Evaluate this contribution's metrics for a given branch length
  ///
  /// Returns OptimizationMetrics containing likelihood, log_likelihood, derivative,
  /// and second_derivative for the given branch length. Both dense and sparse
  /// partitions calculate meaningful likelihood values.
  pub fn evaluate(&self, branch_length: f64) -> OptimizationMetrics {
    match self {
      OptimizationContribution::Dense(contribution) => evaluate_dense_contribution(contribution, branch_length),
      OptimizationContribution::Sparse(contribution) => evaluate_sparse_contribution(contribution, branch_length),
    }
  }

  /// Iterate over this contribution's site patterns as (multiplicity, coefficients).
  ///
  /// Dense and sparse representations differ only in iteration: dense yields
  /// `(1.0, row)` for each alignment position, sparse yields
  /// `(multiplicity, coefficients)` for each compressed site pattern. Callers
  /// that only need per-site coefficient vectors and their weights can consume
  /// both representations through this iterator without dispatching on the
  /// enum variant. The signature matches `evaluate_site_contributions` so the
  /// same iteration contract is used everywhere.
  pub fn sites(&self) -> impl Iterator<Item = (f64, ArrayView1<'_, f64>)> {
    match self {
      OptimizationContribution::Dense(contribution) => {
        Either::Left(contribution.coefficients.outer_iter().map(|row| (1.0, row)))
      },
      OptimizationContribution::Sparse(contribution) => Either::Right(
        contribution
          .site_contributions
          .iter()
          .map(|sc| (sc.multiplicity, sc.coefficients.view())),
      ),
    }
  }

  /// The GTR model backing this contribution.
  ///
  /// Dense and sparse variants each carry their own `GTR` field; this accessor
  /// provides a single interface so code that only reads model-level properties
  /// (eigenvalues, `unimodal_branch_likelihood`) does not need to dispatch on
  /// the variant.
  pub fn gtr(&self) -> &GTR {
    match self {
      OptimizationContribution::Dense(contribution) => &contribution.gtr,
      OptimizationContribution::Sparse(contribution) => &contribution.gtr,
    }
  }

  /// Check whether every site's likelihood at t=0 is positive and finite.
  ///
  /// The per-site likelihood at zero branch length is L_i(0) = sum_c k_{ic},
  /// the sum of eigenvalue-space coefficients. This must be positive (physical:
  /// likelihood is a probability measure) and finite (numerical: representable
  /// in f64) for the zero-branch derivative to be well-defined. If any site
  /// fails, the shortcut declines to decide and falls back to full optimization.
  pub fn all_sites_valid_at_zero(&self) -> bool {
    self.sites().all(|(_, coefficients)| {
      let site_lh = coefficients.sum();
      site_lh > 0.0 && site_lh.is_finite()
    })
  }

  /// Whether this contribution's model has proven unimodal branch-length likelihood.
  ///
  /// When true, a negative derivative at t=0 guarantees zero is the global
  /// maximum on [0, infinity). When false, the derivative test only proves
  /// zero is a local maximum - a better positive-length maximum may exist.
  pub fn has_unimodal_branch_likelihood(&self) -> bool {
    self.gtr().unimodal_branch_likelihood
  }

  /// Derivative of log-likelihood with respect to branch length at t=0.
  ///
  /// For the eigendecomposition-based likelihood L_i(t) = sum_c k_{ic} exp(lambda_c t),
  /// the per-site derivative of log L_i at t=0 is:
  ///
  ///   d/dt log L_i(0) = (sum_c k_{ic} lambda_c) / (sum_c k_{ic})
  ///
  /// This is a weighted average of eigenvalues, weighted by the coefficients.
  /// The total derivative sums over all sites (dense: positions, sparse:
  /// multiplicity-weighted site patterns).
  ///
  /// Caller must verify `all_sites_valid_at_zero()` first. If any site has
  /// L_i(0) <= 0 or non-finite, the division produces NaN/Inf.
  pub fn zero_branch_length_derivative(&self) -> f64 {
    debug_assert!(
      self.all_sites_valid_at_zero(),
      "zero_branch_length_derivative called without verifying all_sites_valid_at_zero"
    );
    let eigvals = &self.gtr().eigvals;
    self
      .sites()
      .map(|(multiplicity, coefficients)| multiplicity * coefficients.dot(eigvals) / coefficients.sum())
      .sum()
  }
}

pub fn evaluate_mixed(contributions: &[OptimizationContribution], branch_length: f64) -> OptimizationMetrics {
  evaluate_mixed_impl(contributions, branch_length, true)
}

pub fn evaluate_mixed_log_lh_only(contributions: &[OptimizationContribution], branch_length: f64) -> f64 {
  evaluate_mixed_impl(contributions, branch_length, false).log_lh
}

fn evaluate_mixed_impl(
  contributions: &[OptimizationContribution],
  branch_length: f64,
  compute_derivatives: bool,
) -> OptimizationMetrics {
  let mut total_metrics = OptimizationMetrics::default();
  for contribution in contributions {
    let metrics = match contribution {
      OptimizationContribution::Dense(c) => evaluate_dense_contribution_impl(c, branch_length, compute_derivatives),
      OptimizationContribution::Sparse(c) => evaluate_sparse_contribution_impl(c, branch_length, compute_derivatives),
    };
    total_metrics.add(&metrics);
  }
  total_metrics
}

/// Evaluate substitution + indel contributions for a given branch length.
pub(crate) fn evaluate_with_indels(
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  branch_length: f64,
) -> OptimizationMetrics {
  let mut metrics = evaluate_mixed(contributions, branch_length);
  metrics.add(&poisson_indel_log_lh(indel_count, indel_rate, branch_length));
  metrics
}

/// Log-likelihood only (no derivatives) for substitution + indel contributions.
pub(crate) fn evaluate_with_indels_log_lh_only(
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  branch_length: f64,
) -> f64 {
  let sub_lh = evaluate_mixed_log_lh_only(contributions, branch_length);
  let indel_lh = poisson_indel_log_lh(indel_count, indel_rate, branch_length).log_lh;
  sub_lh + indel_lh
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

/// Unified optimization function for mixed partition types.
///
/// Main optimization loop that works with both sparse and dense partitions simultaneously.
/// For each edge, it collects contributions from all partitions and optimizes the branch
/// length using the selected method.
pub fn run_optimize_mixed<N, E, P>(
  graph: &Graph<N, E, ()>,
  partitions: &[Arc<RwLock<P>>],
  method: BranchOptMethod,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  P: PartitionOptimizeOps + ?Sized,
{
  let total_length = total_sequence_length(partitions);
  if total_length == 0 {
    return make_error!("Total sequence length across all partitions is zero; cannot optimize branch lengths");
  }

  let indel_rate = estimate_indel_rate(graph, partitions);
  run_optimize_mixed_inner(graph, partitions, method, indel_rate, false)
}

#[cfg(test)]
pub(crate) fn run_optimize_mixed_with_indel_rate<N, E, P>(
  graph: &Graph<N, E, ()>,
  partitions: &[Arc<RwLock<P>>],
  method: BranchOptMethod,
  indel_rate: f64,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  P: PartitionOptimizeOps + ?Sized,
{
  run_optimize_mixed_inner(graph, partitions, method, indel_rate, false)
}

pub(crate) fn run_optimize_mixed_inner<N, E, P>(
  graph: &Graph<N, E, ()>,
  partitions: &[Arc<RwLock<P>>],
  method: BranchOptMethod,
  indel_rate: f64,
  no_indels: bool,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  P: PartitionOptimizeOps + ?Sized,
{
  let total_length = total_sequence_length(partitions);

  if total_length == 0 {
    return make_error!("Total sequence length across all partitions is zero; cannot optimize branch lengths");
  }

  let one_mutation = 1.0 / total_length as f64;

  graph.get_edges().iter().try_for_each(|edge_ref| -> Result<(), Report> {
    let edge_key = edge_ref.read_arc().key();
    let mut edge = edge_ref.write_arc().payload().write_arc();
    let mut branch_length = edge.branch_length().unwrap_or(0.0);

    let contributions: Vec<OptimizationContribution> = partitions
      .iter()
      .map(|partition| partition.read_arc().create_edge_contribution(edge_key))
      .collect::<Result<_, _>>()?;

    let indel_count: usize = if no_indels {
      0
    } else {
      partitions
        .iter()
        .map(|partition| partition.read_arc().edge_indel_count(edge_key))
        .sum()
    };

    // The Poisson log-likelihood derivative diverges at t=0 when k > 0, producing
    // inf/NaN in Newton's method. Use a non-zero starting point for indel-bearing edges.
    if branch_length == 0.0 && indel_count > 0 {
      branch_length = if indel_rate > 0.0 {
        (indel_count as f64 / indel_rate).max(one_mutation)
      } else {
        one_mutation
      };
    }

    // When branch length is zero and any site has non-positive likelihood at t=0
    // (mismatched certain states), evaluating ln(site_lh) or dividing by site_lh
    // produces -inf/inf. Bump to a small positive value so the evaluator operates
    // in the well-defined domain.
    if branch_length == 0.0 && !contributions.iter().all(|c| c.all_sites_valid_at_zero()) {
      branch_length = one_mutation;
    }

    // When indels are present on this edge, the Poisson derivative at t=0 is +infinity,
    // so zero branch length is never optimal. Only check the substitution-based criterion
    // when there are no indels.
    if indel_count == 0 && is_zero_branch_optimal(&contributions) {
      edge.set_branch_length(Some(0.0));
      return Ok(());
    }

    // Lower bound for Newton/Brent steps on indel-bearing edges. The Poisson derivative
    // diverges at t=0, so we must prevent the optimizer from landing exactly at zero.
    let min_branch_length = min_branch_length_for_indels(indel_count, one_mutation);

    let new_branch_length = match method {
      BranchOptMethod::Brent => brent_inner(
        branch_length,
        &contributions,
        indel_count,
        indel_rate,
        min_branch_length,
        one_mutation,
      ),
      BranchOptMethod::BrentSqrt => brent_sqrt_inner(
        branch_length,
        &contributions,
        indel_count,
        indel_rate,
        min_branch_length,
        one_mutation,
      ),
      BranchOptMethod::BrentLog => brent_log_inner(
        branch_length,
        &contributions,
        indel_count,
        indel_rate,
        min_branch_length,
        one_mutation,
      ),
      BranchOptMethod::Newton => {
        let metrics = evaluate_with_indels(&contributions, indel_count, indel_rate, branch_length);
        newton_inner(
          branch_length,
          &metrics,
          &contributions,
          indel_count,
          indel_rate,
          min_branch_length,
          one_mutation,
        )
      },
      BranchOptMethod::NewtonSqrt => {
        let metrics = evaluate_with_indels(&contributions, indel_count, indel_rate, branch_length);
        newton_sqrt_inner(
          branch_length,
          &metrics,
          &contributions,
          indel_count,
          indel_rate,
          min_branch_length,
          one_mutation,
        )
      },
      BranchOptMethod::NewtonLog => {
        // ln(t) requires t > 0; bump zero branch lengths to one_mutation
        let bl = if branch_length == 0.0 {
          one_mutation
        } else {
          branch_length
        };
        let metrics = evaluate_with_indels(&contributions, indel_count, indel_rate, bl);
        newton_log_inner(
          bl,
          &metrics,
          &contributions,
          indel_count,
          indel_rate,
          min_branch_length,
          one_mutation,
        )
      },
    }?;

    // Post-optimization boundary reconciliation. See `reconcile_zero_boundary`
    // for the full rationale. The input `branch_length` is used to size the
    // verification grid (the optimizer's output may be clamped to zero or to a
    // tiny floor like $10^{-12}$ and is unsuitable as an extent).
    let new_branch_length = reconcile_zero_boundary(
      new_branch_length,
      branch_length,
      &contributions,
      indel_count,
      indel_rate,
      one_mutation,
    )?;

    edge.set_branch_length(Some(new_branch_length));
    Ok(())
  })
}

/// Initial estimation of branch lengths for mixed partitions.
///
/// Computes per-edge substitution count over canonical (non-ambiguous,
/// non-deletion) positions via `edge_subs().len()` for both sparse and
/// dense partitions.
///
/// The denominator is the per-edge effective alignment length rather than the raw
/// sequence length, so gap-heavy edges get correctly scaled rates.
///
/// For edges with indels but no substitutions, the Poisson maximum likelihood estimate (MLE) $\hat{t} = k / \mu$
/// provides a non-zero initial estimate so Newton optimization starts from a
/// reasonable point (the indel derivative diverges at $t = 0$).
///
/// When `overwrite_valid` is false, edges that already have a finite branch
/// length are skipped, preserving calibrated input values while filling in
/// edges with missing (`None`) or invalid (`NaN`) branch lengths.
pub fn initial_guess_mixed<N, E, P>(
  graph: &Graph<N, E, ()>,
  partitions: &[Arc<RwLock<P>>],
  overwrite_valid: bool,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  P: PartitionOptimizeOps + ?Sized,
{
  let total_length: usize = partitions
    .iter()
    .map(|partition| partition.read_arc().sequence_length())
    .sum();

  if total_length == 0 {
    return make_error!("Total sequence length across all partitions is zero; cannot compute initial guess");
  }

  let one_mutation = 1.0 / total_length as f64;
  let indel_rate = estimate_indel_rate(graph, partitions);

  for edge_ref in graph.get_edges() {
    let edge_key = edge_ref.read_arc().key();
    let mut edge = edge_ref.write_arc().payload().write_arc();

    let indel_count: usize = partitions
      .iter()
      .map(|partition| partition.read_arc().edge_indel_count(edge_key))
      .sum();

    if !overwrite_valid {
      if let Some(bl) = edge.branch_length() {
        // A finite positive BL is always valid. A zero BL is valid only
        // if the edge has no indels: with indels present, the Poisson
        // derivative diverges at t=0 and estimate_indel_rate() needs
        // positive total BL to produce a nonzero rate.
        if bl.is_finite() && (bl > 0.0 || indel_count == 0) {
          continue;
        }
      }
    }

    let sub_count: usize = partitions
      .iter()
      .map(|partition| partition.read_arc().edge_subs(graph, edge_key).map(|subs| subs.len()))
      .sum::<Result<_, _>>()?;

    let effective_length: usize = partitions
      .iter()
      .map(|partition| partition.read_arc().edge_effective_length(graph, edge_key))
      .sum::<Result<_, _>>()?;

    let branch_length = if effective_length > 0 {
      let sub_estimate = sub_count as f64 / effective_length as f64;
      if sub_estimate == 0.0 && indel_count > 0 {
        if indel_rate > 0.0 {
          // Poisson MLE for indel-only branches: t = k / mu
          indel_count as f64 / indel_rate
        } else {
          // Bootstrap: rate unknown (e.g. all-zero input tree), seed with one_mutation
          // so that estimate_indel_rate produces a positive rate in run_optimize_mixed
          one_mutation
        }
      } else {
        sub_estimate
      }
    } else if indel_count > 0 {
      // All positions are gaps/ambiguous but indels are present
      one_mutation
    } else {
      one_mutation * 0.1
    };
    edge.set_branch_length(Some(branch_length));
  }

  Ok(())
}

fn total_sequence_length<P>(partitions: &[Arc<RwLock<P>>]) -> usize
where
  P: PartitionOptimizeOps + ?Sized,
{
  partitions
    .iter()
    .map(|partition| partition.read_arc().sequence_length())
    .sum()
}
