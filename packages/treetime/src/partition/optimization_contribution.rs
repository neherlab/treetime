use crate::gtr::gtr::GTR;
use crate::partition::marginal_dense::PartitionMarginalDense;
use crate::partition::marginal_sparse::PartitionMarginalSparse;
use crate::partition::optimize_dense;
use crate::partition::optimize_sparse;
use eyre::Report;
use itertools::Either;
use ndarray::ArrayView1;
use treetime_graph::edge::GraphEdgeKey;

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
    let edge_partition = &partition.data.edges[&edge_key];
    let contribution = optimize_dense::get_coefficients(
      &edge_partition.msg_to_parent,
      &edge_partition.msg_to_child,
      &partition.data.gtr,
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
  pub fn gtr(&self) -> &GTR {
    match self {
      OptimizationContribution::Dense(contribution) => &contribution.gtr,
      OptimizationContribution::Sparse(contribution) => &contribution.gtr,
    }
  }

  /// Check whether every site's likelihood at t=0 is positive and finite.
  pub fn all_sites_valid_at_zero(&self) -> bool {
    self.sites().all(|(_, coefficients)| {
      let site_lh = coefficients.sum();
      site_lh > 0.0 && site_lh.is_finite()
    })
  }

  /// Whether this contribution's model has proven unimodal branch-length likelihood.
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
  /// Caller must verify `all_sites_valid_at_zero()` first.
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
