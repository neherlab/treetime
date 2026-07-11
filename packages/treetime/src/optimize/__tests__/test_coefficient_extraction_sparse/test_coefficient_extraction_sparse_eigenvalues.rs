#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::sparse_eval::evaluate_sparse_contribution;
  use crate::partition::optimize_sparse::{PartitionContribution, SiteContribution};
  use ndarray::array;

  #[test]
  fn test_eigenvalues_affect_branch_length_evaluation() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let site = SiteContribution {
      multiplicity: 1.0,
      coefficients: array![0.5, 0.2, 0.2, 0.1],
    };

    let contribution = PartitionContribution {
      site_contributions: vec![site],
      gtr,
    };

    // At different branch lengths, exp(λt) changes for non-zero eigenvalues
    let metrics_short = evaluate_sparse_contribution(&contribution, 0.01).expect("valid branch length");
    let metrics_long = evaluate_sparse_contribution(&contribution, 1.0).expect("valid branch length");

    // Log-LH should differ at different branch lengths
    assert!(
      (metrics_short.log_lh - metrics_long.log_lh).abs() > 1e-6,
      "log-LH should differ at different branch lengths"
    );
  }

  #[test]
  fn test_jc69_eigenvalues_structure() {
    // JC69 has one zero eigenvalue and three equal negative eigenvalues
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    // Find the zero eigenvalue (should be one)
    let zero_count = gtr.eigvals.iter().filter(|&&ev| ev.abs() < 1e-10).count();
    assert_eq!(zero_count, 1, "JC69 should have exactly one zero eigenvalue");

    // All eigenvalues should be <= 0
    for &ev in &gtr.eigvals {
      assert!(ev <= 1e-10, "JC69 eigenvalues should be non-positive");
    }
  }
}
