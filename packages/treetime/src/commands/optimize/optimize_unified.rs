use crate::commands::optimize::optimize_dense;
use crate::commands::optimize::optimize_dense_eval::{evaluate_dense_contribution, evaluate_dense_contribution_impl};
use crate::commands::optimize::optimize_indel::{estimate_indel_rate, poisson_indel_log_lh};
use crate::commands::optimize::optimize_sparse;
use crate::commands::optimize::optimize_sparse_eval::{
  evaluate_sparse_contribution, evaluate_sparse_contribution_impl,
};
use crate::commands::optimize::partition_ops::PartitionOptimizeOps;
use crate::representation::partition::marginal_dense::PartitionMarginalDense;
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::representation::payload::ancestral::GraphAncestral;
use eyre::Report;
use ndarray::Axis;
use num::clamp;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::edge::{GraphEdgeKey, HasBranchLength};

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
      OptimizationContribution::Sparse(contribution) => contribution.unimodal_branch_likelihood,
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
          coeff.multiplicity * ((&coeff.coefficients * &contribution.eigenvalues).sum() / coeff.coefficients.sum())
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

/// Unified optimization function for mixed partition types
///
/// Main optimization loop that works with both sparse and dense partitions simultaneously.
/// For each edge, it collects contributions from all partitions and optimizes the branch
/// length using either Newton's method or grid search.
pub fn run_optimize_mixed<P>(graph: &GraphAncestral, partitions: &[Arc<RwLock<P>>]) -> Result<(), Report>
where
  P: PartitionOptimizeOps + ?Sized,
{
  let total_length: usize = partitions
    .iter()
    .map(|partition| partition.read_arc().sequence_length())
    .sum();

  let one_mutation = 1.0 / total_length as f64;
  let indel_rate = estimate_indel_rate(graph, partitions);

  graph.get_edges().iter().try_for_each(|edge_ref| -> Result<(), Report> {
    let edge_key = edge_ref.read_arc().key();
    let mut edge = edge_ref.write_arc().payload().write_arc();
    let mut branch_length = edge.branch_length().unwrap_or(0.0);

    let contributions: Vec<OptimizationContribution> = partitions
      .iter()
      .map(|partition| partition.read_arc().create_edge_contribution(edge_key))
      .collect::<Result<_, _>>()?;

    let indel_count: usize = partitions
      .iter()
      .map(|partition| partition.read_arc().edge_indel_count(edge_key))
      .sum();

    // When indels are present on this edge, the Poisson derivative at t=0 is +infinity,
    // so zero branch length is never optimal. Only check the substitution-based criterion
    // when there are no indels.
    if indel_count == 0 && is_zero_branch_optimal(&contributions) {
      edge.set_branch_length(Some(0.0));
      return Ok(());
    }

    let metrics = evaluate_with_indels(&contributions, indel_count, indel_rate, branch_length);
    let mut new_branch_length;

    if metrics.second_derivative < 0.0 {
      // Newton's method to find the optimal branch length
      new_branch_length = branch_length - clamp(metrics.derivative / metrics.second_derivative, -1.0, branch_length);
      let max_iter = 10;
      let mut n_iter = 0;

      while (new_branch_length - branch_length).abs() > 0.001 * branch_length && n_iter < max_iter {
        let new_metrics = evaluate_with_indels(&contributions, indel_count, indel_rate, new_branch_length);
        if new_metrics.second_derivative < 0.0 {
          branch_length = new_branch_length;
          new_branch_length = branch_length
            - clamp(
              new_metrics.derivative / new_metrics.second_derivative,
              -1.0,
              branch_length,
            );
        } else {
          break;
        }
        n_iter += 1;
      }
    } else {
      // Grid search over positive branch lengths
      let branch_lengths = ndarray::Array1::linspace(0.1 * one_mutation, 1.5 * branch_length + one_mutation, 100);

      let best_positive = branch_lengths
        .iter()
        .max_by_key(|&&bl| {
          let log_lh = evaluate_with_indels_log_lh_only(&contributions, indel_count, indel_rate, bl);
          ordered_float::OrderedFloat(log_lh)
        })
        .copied()
        .unwrap();

      // Evaluate t=0 as a separate candidate: for non-unimodal models that
      // bypass the derivative shortcut, zero may still be optimal. Guard with
      // site validity to avoid ln(0) for degenerate coefficient rows.
      let zero_is_better = contributions.iter().all(|c| c.all_sites_valid_at_zero()) && {
        let log_lh_zero = evaluate_mixed_log_lh_only(&contributions, 0.0);
        let log_lh_best = evaluate_mixed_log_lh_only(&contributions, best_positive);
        log_lh_zero > log_lh_best
      };
      new_branch_length = if zero_is_better { 0.0 } else { best_positive };
    }

    edge.set_branch_length(Some(new_branch_length));
    Ok(())
  })
}

/// Initial estimation of branch lengths for mixed partitions.
///
/// Computes per-edge substitution count over canonical (non-gap, non-ambiguous)
/// positions via `edge_subs().len()` for both sparse and dense partitions.
///
/// The denominator is the per-edge effective alignment length rather than the raw
/// sequence length, so gap-heavy edges get correctly scaled rates.
///
/// For edges with indels but no substitutions, the Poisson MLE $\hat{t} = k / \mu$
/// provides a non-zero initial estimate so Newton optimization starts from a
/// reasonable point (the indel derivative diverges at $t = 0$).
pub fn initial_guess_mixed<P>(graph: &GraphAncestral, partitions: &[Arc<RwLock<P>>]) -> Result<(), Report>
where
  P: PartitionOptimizeOps + ?Sized,
{
  let indel_rate = estimate_indel_rate(graph, partitions);

  for edge_ref in graph.get_edges() {
    let edge_key = edge_ref.read_arc().key();
    let mut edge = edge_ref.write_arc().payload().write_arc();

    let sub_count: usize = partitions
      .iter()
      .map(|partition| partition.read_arc().edge_subs(graph, edge_key).map(|subs| subs.len()))
      .sum::<Result<_, _>>()?;

    let effective_length: usize = partitions
      .iter()
      .map(|partition| partition.read_arc().edge_effective_length(graph, edge_key))
      .sum::<Result<_, _>>()?;

    let indel_count: usize = partitions
      .iter()
      .map(|partition| partition.read_arc().edge_indel_count(edge_key))
      .sum();

    let branch_length = if effective_length > 0 {
      let sub_estimate = sub_count as f64 / effective_length as f64;
      if sub_estimate == 0.0 && indel_count > 0 && indel_rate > 0.0 {
        // Poisson MLE for indel-only branches: t = k / mu
        indel_count as f64 / indel_rate
      } else {
        sub_estimate
      }
    } else {
      0.0
    };

    edge.set_branch_length(Some(branch_length));
  }

  Ok(())
}
