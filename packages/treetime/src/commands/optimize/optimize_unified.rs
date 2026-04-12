use crate::commands::optimize::optimize_dense;
use crate::commands::optimize::optimize_dense_eval::{evaluate_dense_contribution, evaluate_dense_contribution_impl};
use crate::commands::optimize::optimize_indel::{estimate_indel_rate, poisson_indel_log_lh};
use crate::commands::optimize::optimize_method::BranchOptMethod;
use crate::commands::optimize::optimize_sparse;
use crate::commands::optimize::optimize_sparse_eval::{
  evaluate_sparse_contribution, evaluate_sparse_contribution_impl,
};
use crate::commands::optimize::partition_ops::PartitionOptimizeOps;
use crate::make_error;
use crate::representation::partition::marginal_dense::PartitionMarginalDense;
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::representation::partition::traits::{ExactStateCache, GraphNodePathLookup};
use crate::representation::payload::ancestral::GraphAncestral;
use argmin::core::{CostFunction, Error, Executor};
use argmin::solver::brent::BrentOpt;
use eyre::Report;
use ndarray::{Array1, Axis};
use num::clamp;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::edge::{GraphEdgeKey, HasBranchLength};

/// Number of points in the grid search fallback.
const GRID_SEARCH_POINTS: usize = 100;

/// Minimum upper bound (subs/site) for the grid search fallback.
/// Ensures the grid covers the biologically plausible range regardless
/// of the current branch length estimate. Branches longer than 0.5
/// subs/site are rare in real phylogenies; the proportional upper
/// bound `1.5 * branch_length + one_mutation` extends beyond this
/// when the current estimate is already large.
const GRID_SEARCH_MIN_UPPER: f64 = 0.5;

/// Relative tolerance for Newton inner-loop convergence.
const NEWTON_REL_TOL: f64 = 0.001;

/// Absolute tolerance floor for Newton inner-loop convergence (subs/site).
/// Prevents the purely relative tolerance from degenerating to zero when
/// the current branch length is zero or very small. The value 1e-8 is well
/// below the smallest meaningful branch length (1/L ~ 1e-4 for L=10000)
/// so it does not mask real convergence.
const NEWTON_ABS_TOL: f64 = 1e-8;

/// Maximum iterations for Newton inner loop.
const NEWTON_MAX_ITER: usize = 10;

/// Maximum iterations for Brent's method.
const BRENT_MAX_ITER: u64 = 50;

/// Newton inner-loop convergence tolerance for branch length updates.
///
/// Combines relative tolerance (0.1% of current branch length) with an
/// absolute floor to avoid degenerate behavior at zero branch length.
pub(crate) fn newton_tolerance(branch_length: f64) -> f64 {
  f64::max(NEWTON_REL_TOL * branch_length, NEWTON_ABS_TOL)
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
  pub fn from_sparse(
    graph: &dyn GraphNodePathLookup,
    edge_key: GraphEdgeKey,
    partition: &PartitionMarginalSparse,
    cache: &mut ExactStateCache,
  ) -> Result<Self, Report> {
    let contribution = optimize_sparse::get_coefficients(graph, edge_key, partition, cache)?;
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

  /// Check whether every site's likelihood at t=0 is positive and finite.
  ///
  /// The per-site likelihood at zero branch length is L_i(0) = sum_c k_{ic},
  /// the sum of eigenvalue-space coefficients. This must be positive (physical:
  /// likelihood is a probability measure) and finite (numerical: representable
  /// in f64) for the zero-branch derivative to be well-defined. If any site
  /// fails, the shortcut declines to decide and falls back to full optimization.
  pub fn all_sites_valid_at_zero(&self) -> bool {
    match self {
      OptimizationContribution::Dense(contribution) => contribution
        .coefficients
        .sum_axis(Axis(1))
        .iter()
        .all(|&site_lh| site_lh > 0.0 && site_lh.is_finite()),
      OptimizationContribution::Sparse(contribution) => contribution.site_contributions.iter().all(|coeff| {
        let site_lh = coeff.coefficients.sum();
        site_lh > 0.0 && site_lh.is_finite()
      }),
    }
  }

  /// Whether this contribution's model has proven unimodal branch-length likelihood.
  ///
  /// When true, a negative derivative at t=0 guarantees zero is the global
  /// maximum on [0, infinity). When false, the derivative test only proves
  /// zero is a local maximum - a better positive-length maximum may exist.
  pub fn has_unimodal_branch_likelihood(&self) -> bool {
    match self {
      OptimizationContribution::Dense(contribution) => contribution.gtr.unimodal_branch_likelihood,
      OptimizationContribution::Sparse(contribution) => contribution.gtr.unimodal_branch_likelihood,
    }
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
    match self {
      OptimizationContribution::Dense(contribution) => {
        let gtr = &contribution.gtr;
        let coefficients = &contribution.coefficients;
        ((coefficients * &gtr.eigvals).sum_axis(Axis(1)) / coefficients.sum_axis(Axis(1))).sum()
      },
      OptimizationContribution::Sparse(contribution) => contribution
        .site_contributions
        .iter()
        .map(|coeff| {
          coeff.multiplicity * ((&coeff.coefficients * &contribution.gtr.eigvals).sum() / coeff.coefficients.sum())
        })
        .sum(),
    }
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
fn evaluate_with_indels(
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
fn evaluate_with_indels_log_lh_only(
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
pub(crate) fn grid_search_branch_lengths(branch_length: f64, one_mutation: f64) -> Array1<f64> {
  let lower = 0.1 * one_mutation;
  let upper = f64::max(1.5 * branch_length + one_mutation, GRID_SEARCH_MIN_UPPER);
  Array1::linspace(lower.ln(), upper.ln(), GRID_SEARCH_POINTS).mapv(f64::exp)
}

/// Newton-Raphson inner loop in $t$-space.
///
/// Takes an initial `branch_length` and `metrics` (already evaluated at that point),
/// runs up to `NEWTON_MAX_ITER` Newton steps with the convergence criterion
/// $|t_{\text{new}} - t_{\text{old}}| < \text{tol}(t)$, and returns the optimized
/// branch length. Falls through to grid search if the Hessian becomes non-negative.
fn newton_inner(
  branch_length: f64,
  metrics: &OptimizationMetrics,
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  min_branch_length: f64,
  one_mutation: f64,
) -> f64 {
  if metrics.second_derivative < 0.0 {
    let mut bl = branch_length;
    let mut new_bl = (bl - clamp(metrics.derivative / metrics.second_derivative, -1.0, bl)).max(min_branch_length);

    for _ in 0..NEWTON_MAX_ITER {
      if (new_bl - bl).abs() <= newton_tolerance(bl) {
        break;
      }
      let new_metrics = evaluate_with_indels(contributions, indel_count, indel_rate, new_bl);
      if new_metrics.second_derivative < 0.0 {
        bl = new_bl;
        new_bl = (bl - clamp(new_metrics.derivative / new_metrics.second_derivative, -1.0, bl)).max(min_branch_length);
      } else {
        break;
      }
    }
    new_bl
  } else {
    grid_search_inner(branch_length, contributions, indel_count, indel_rate, one_mutation)
  }
}

/// Newton-Raphson inner loop in $\sqrt{t}$ space.
///
/// Reparameterizes the optimization variable as $s = \sqrt{t}$ and applies the
/// chain rule to transform $t$-space derivatives into $s$-space:
///
///   $d\ell/ds = 2s \cdot d\ell/dt$
///   $d^2\ell/ds^2 = 4s^2 \cdot d^2\ell/dt^2 + 2 \cdot d\ell/dt$
///
/// This reduces the indel Hessian singularity from $O(1/t^2)$ to $O(1/t)$,
/// allowing Newton to make progress toward the combined (substitution + indel)
/// optimum instead of stalling at the indel-only MLE.
///
/// The convergence criterion is applied in $s$-space:
/// $|s_{\text{new}} - s_{\text{old}}| < \text{tol}(s)$.
fn newton_sqrt_inner(
  branch_length: f64,
  metrics: &OptimizationMetrics,
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  min_branch_length: f64,
  one_mutation: f64,
) -> f64 {
  let mut s = branch_length.sqrt();
  let min_s = min_branch_length.sqrt();

  // Transform initial metrics to s-space
  let (ds, d2s) = chain_rule_sqrt(s, metrics.derivative, metrics.second_derivative);

  if d2s >= 0.0 {
    return grid_search_inner(branch_length, contributions, indel_count, indel_rate, one_mutation);
  }

  let mut new_s = (s - clamp(ds / d2s, -1.0, s)).max(min_s);

  for _ in 0..NEWTON_MAX_ITER {
    if (new_s - s).abs() <= newton_tolerance(s) {
      break;
    }
    let t = new_s * new_s;
    let new_metrics = evaluate_with_indels(contributions, indel_count, indel_rate, t);
    let (ds_new, d2s_new) = chain_rule_sqrt(new_s, new_metrics.derivative, new_metrics.second_derivative);
    if d2s_new < 0.0 {
      s = new_s;
      new_s = (s - clamp(ds_new / d2s_new, -1.0, s)).max(min_s);
    } else {
      break;
    }
  }

  new_s * new_s
}

/// Transform $t$-space derivatives to $\sqrt{t}$-space via the chain rule.
///
/// Given $s = \sqrt{t}$:
///   $d\ell/ds = 2s \cdot d\ell/dt$
///   $d^2\ell/ds^2 = 4s^2 \cdot d^2\ell/dt^2 + 2 \cdot d\ell/dt$
pub(crate) fn chain_rule_sqrt(s: f64, dl_dt: f64, d2l_dt2: f64) -> (f64, f64) {
  let dl_ds = 2.0 * s * dl_dt;
  let d2l_ds2 = 4.0 * s * s * d2l_dt2 + 2.0 * dl_dt;
  (dl_ds, d2l_ds2)
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
fn brent_inner(
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

/// Grid search fallback for non-concave regions.
fn grid_search_inner(
  branch_length: f64,
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  one_mutation: f64,
) -> f64 {
  let branch_lengths = grid_search_branch_lengths(branch_length, one_mutation);

  let best_positive = branch_lengths
    .iter()
    .max_by_key(|&&bl| {
      let log_lh = evaluate_with_indels_log_lh_only(contributions, indel_count, indel_rate, bl);
      ordered_float::OrderedFloat(log_lh)
    })
    .copied()
    .unwrap();

  let zero_is_better = is_zero_better_than_grid_best(contributions, indel_count, indel_rate, best_positive);
  if zero_is_better { 0.0 } else { best_positive }
}

/// Unified optimization function for mixed partition types.
///
/// Main optimization loop that works with both sparse and dense partitions simultaneously.
/// For each edge, it collects contributions from all partitions and optimizes the branch
/// length using the selected method.
pub fn run_optimize_mixed<P>(
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<P>>],
  method: BranchOptMethod,
) -> Result<(), Report>
where
  P: PartitionOptimizeOps + ?Sized,
{
  let total_length: usize = partitions
    .iter()
    .map(|partition| partition.read_arc().sequence_length())
    .sum();

  if total_length == 0 {
    return make_error!("Total sequence length across all partitions is zero; cannot optimize branch lengths");
  }

  let one_mutation = 1.0 / total_length as f64;
  let indel_rate = estimate_indel_rate(graph, partitions);
  graph.get_edges().iter().try_for_each(|edge_ref| -> Result<(), Report> {
    let edge_key = edge_ref.read_arc().key();
    let mut edge = edge_ref.write_arc().payload().write_arc();
    let mut branch_length = edge.branch_length().unwrap_or(0.0);

    let mut contributions = Vec::with_capacity(partitions.len());
    for partition in partitions {
      let mut exact_state_cache = ExactStateCache::new();
      contributions.push(
        partition
          .read_arc()
          .create_edge_contribution(graph, edge_key, &mut exact_state_cache)?,
      );
    }

    let indel_count: usize = partitions
      .iter()
      .map(|partition| partition.read_arc().edge_indel_count(edge_key))
      .sum();

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
    let min_branch_length = if indel_count > 0 { one_mutation * 0.01 } else { 0.0 };

    let new_branch_length = match method {
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
      BranchOptMethod::Brent => brent_inner(
        branch_length,
        &contributions,
        indel_count,
        indel_rate,
        min_branch_length,
        one_mutation,
      ),
    };

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
pub fn initial_guess_mixed<P>(
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<P>>],
  overwrite_valid: bool,
) -> Result<(), Report>
where
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
      0.0
    };

    edge.set_branch_length(Some(branch_length));
  }

  Ok(())
}
