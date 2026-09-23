use crate::gtr::gtr::GTR;
use crate::partition::optimize;
use crate::partition::storage::dense::{DenseEdgeBackward, DenseEdgeForward};
use crate::partition::storage::sparse::{SparseEdgeBackward, SparseEdgeForward, SparseEdgeObs};
use eyre::Report;
use itertools::Either;
use ndarray::ArrayView1;

pub enum OptimizationContribution {
  Dense(optimize::dense::PartitionContribution),
  Sparse(optimize::sparse::PartitionContribution),
}

impl OptimizationContribution {
  pub(crate) fn from_dense(gtr: &GTR, backward: &DenseEdgeBackward, forward: &DenseEdgeForward) -> Self {
    let contribution = optimize::dense::get_coefficients(&backward.msg_to_parent, &forward.msg_to_child, gtr);
    OptimizationContribution::Dense(contribution)
  }

  pub(crate) fn from_sparse(
    gtr: &GTR,
    backward: &SparseEdgeBackward,
    forward: &SparseEdgeForward,
    edge_obs: &SparseEdgeObs,
  ) -> Result<Self, Report> {
    let contribution = optimize::sparse::get_coefficients(gtr, backward, forward, edge_obs)?;
    Ok(OptimizationContribution::Sparse(contribution))
  }

  fn sites(&self) -> impl Iterator<Item = (f64, ArrayView1<'_, f64>)> {
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

  fn gtr(&self) -> &GTR {
    match self {
      OptimizationContribution::Dense(contribution) => &contribution.gtr,
      OptimizationContribution::Sparse(contribution) => &contribution.gtr,
    }
  }

  pub(crate) fn all_sites_valid_at_zero(&self) -> bool {
    self.sites().all(|(_, coefficients)| {
      let site_lh = coefficients.sum();
      site_lh > 0.0 && site_lh.is_finite()
    })
  }

  pub(crate) fn has_unimodal_branch_likelihood(&self) -> bool {
    self.gtr().unimodal_branch_likelihood
  }

  pub(crate) fn zero_branch_length_derivative(&self) -> f64 {
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
