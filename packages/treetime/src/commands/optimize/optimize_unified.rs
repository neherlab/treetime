use crate::commands::optimize::optimize_dense;
use crate::commands::optimize::optimize_dense_eval::{evaluate_dense_contribution, evaluate_dense_contribution_impl};
use crate::commands::optimize::optimize_sparse;
use crate::commands::optimize::optimize_sparse_eval::{
  evaluate_sparse_contribution, evaluate_sparse_contribution_impl,
};
use crate::representation::partition::marginal_dense::PartitionMarginalDense;
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::representation::payload::ancestral::GraphAncestral;
use eyre::Report;
use itertools::{chain, izip};
use ndarray::Axis;
use ndarray_stats::QuantileExt;
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

  /// Check if zero branch length is optimal for this contribution
  ///
  /// Returns the likelihood value when branch length is zero, which can be used
  /// to determine if optimization should converge to zero length.
  pub fn zero_branch_length_likelihood(&self) -> f64 {
    match self {
      OptimizationContribution::Dense(contribution) => contribution.coefficients.sum_axis(Axis(1)).product(),
      OptimizationContribution::Sparse(contribution) => contribution
        .site_contributions
        .iter()
        .map(|coeff| coeff.coefficients.sum().powf(coeff.multiplicity))
        .product(),
    }
  }

  /// Get zero branch length derivative for this contribution
  ///
  /// Returns the derivative of log likelihood with respect to branch length
  /// when branch length is zero. Negative values indicate zero is optimal.
  pub fn zero_branch_length_derivative(&self) -> f64 {
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

/// Check if zero branch length is optimal across all contributions
///
/// Analyzes the combined likelihood and derivative at zero branch length to
/// determine if the optimization should converge to zero.
pub fn is_zero_branch_optimal(contributions: &[OptimizationContribution]) -> bool {
  // Calculate combined zero branch length likelihood
  let zero_branch_length_lh: f64 = contributions
    .iter()
    .map(|contrib| contrib.zero_branch_length_likelihood())
    .product();

  if zero_branch_length_lh > 0.01 {
    // Check if derivative is negative (indicating zero is optimal)
    let zero_branch_length_derivative: f64 = contributions
      .iter()
      .map(|contrib| contrib.zero_branch_length_derivative())
      .sum();

    if zero_branch_length_derivative < 0.0 {
      return true;
    }
  }

  false
}

/// Unified optimization function for mixed partition types
///
/// Main optimization loop that works with both sparse and dense partitions simultaneously.
/// For each edge, it collects contributions from all partitions and optimizes the branch
/// length using either Newton's method or grid search.
pub fn run_optimize_mixed(
  graph: &GraphAncestral,
  dense_partitions: &[Arc<RwLock<PartitionMarginalDense>>],
  sparse_partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
) -> Result<(), Report> {
  let dense_lengths = dense_partitions.iter().map(|p| p.read_arc().get_sequence_length());
  let sparse_lengths = sparse_partitions.iter().map(|p| p.read_arc().get_sequence_length());
  let total_length: usize = chain!(dense_lengths, sparse_lengths).sum();

  let one_mutation = 1.0 / total_length as f64;

  graph.get_edges().iter().try_for_each(|edge_ref| -> Result<(), Report> {
    let edge_key = edge_ref.read_arc().key();
    let mut edge = edge_ref.write_arc().payload().write_arc();
    let mut branch_length = edge.branch_length().unwrap_or(0.0);

    // Collect contributions from all partitions
    let contributions = {
      let mut contributions = Vec::with_capacity(dense_partitions.len() + sparse_partitions.len());

      // Add dense partition contributions
      for partition in dense_partitions {
        let partition = partition.read_arc();
        contributions.push(OptimizationContribution::from_dense(edge_key, &partition));
      }

      // Add sparse partition contributions
      for partition in sparse_partitions {
        let partition = partition.read_arc();
        contributions.push(OptimizationContribution::from_sparse(edge_key, &partition)?);
      }

      contributions
    };

    // Check if zero branch length is optimal
    if is_zero_branch_optimal(&contributions) {
      edge.set_branch_length(Some(0.0));
      return Ok(());
    }

    // Otherwise, optimize the branch length
    let metrics = evaluate_mixed(&contributions, branch_length);
    let mut new_branch_length;

    if metrics.second_derivative < 0.0 {
      // Newton's method to find the optimal branch length
      new_branch_length = branch_length - clamp(metrics.derivative / metrics.second_derivative, -1.0, branch_length);
      let max_iter = 10;
      let mut n_iter = 0;

      while (new_branch_length - branch_length).abs() > 0.001 * branch_length && n_iter < max_iter {
        let new_metrics = evaluate_mixed(&contributions, new_branch_length);
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
      // Evaluate on a vector of branch lengths to find the maximum
      let branch_lengths = ndarray::Array1::linspace(0.1 * one_mutation, 1.5 * branch_length + one_mutation, 100);

      let best_branch_length = branch_lengths
        .iter()
        .max_by_key(|&&bl| {
          let metrics = evaluate_mixed(&contributions, bl);
          ordered_float::OrderedFloat(metrics.log_lh)
        })
        .copied()
        .unwrap();
      new_branch_length = best_branch_length;
    }

    edge.set_branch_length(Some(new_branch_length));
    Ok(())
  })
}

/// Initial estimation of branch lengths for mixed partitions
///
/// Provides an initial estimate of branch lengths based on the number of observed
/// differences across all partitions. This gives the optimization a reasonable
/// starting point.
pub fn initial_guess_mixed(
  graph: &GraphAncestral,
  dense_partitions: &[Arc<RwLock<PartitionMarginalDense>>],
  sparse_partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
) {
  for edge_ref in graph.get_edges() {
    let edge_key = edge_ref.read_arc().key();
    let mut edge = edge_ref.write_arc().payload().write_arc();
    let mut differences: usize = 0;

    // Count differences from dense partitions
    for partition in dense_partitions {
      let partition = partition.read_arc();
      let edge = &partition.edges[&edge_key];
      let rows_parent = edge.msg_to_parent.dis.rows();
      let rows_child = edge.msg_to_child.dis.rows();
      for (prof_parent, prof_child) in izip!(rows_parent, rows_child) {
        if prof_parent[prof_child.argmax().unwrap()] < 0.5 {
          differences += 1;
        }
      }
    }

    // Count differences from sparse partitions
    for partition in sparse_partitions {
      let partition = partition.read_arc();
      let edge_partition = &partition.edges[&edge_key];
      differences += edge_partition.subs.len();
    }

    let dense_lengths = dense_partitions.iter().map(|p| p.read_arc().get_sequence_length());
    let sparse_lengths = sparse_partitions.iter().map(|p| p.read_arc().get_sequence_length());
    let total_length: usize = chain!(dense_lengths, sparse_lengths).sum();

    let branch_length = (differences as f64) / (total_length as f64);

    edge.set_branch_length(Some(branch_length));
  }
}
