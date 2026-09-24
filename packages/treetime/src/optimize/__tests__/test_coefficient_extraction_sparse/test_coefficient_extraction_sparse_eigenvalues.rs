#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::sparse_eval::evaluate_sparse_contribution;
  use crate::partition::optimize::sparse::{PartitionContribution, SiteContribution};
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

    let metrics_short = evaluate_sparse_contribution(&contribution, 0.01, true).expect("valid branch length");
    let metrics_long = evaluate_sparse_contribution(&contribution, 1.0, true).expect("valid branch length");

    assert!(
      (metrics_short.log_lh.value() - metrics_long.log_lh.value()).abs() > 1e-6,
      "log-LH should differ at different branch lengths"
    );
  }

  #[test]
  fn test_jc69_eigenvalues_structure() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let zero_count = gtr.eigvals.iter().filter(|&&ev| ev.abs() < 1e-10).count();
    assert_eq!(1, zero_count, "JC69 should have exactly one zero eigenvalue");

    for &ev in &gtr.eigvals {
      assert!(ev <= 1e-10, "JC69 eigenvalues should be non-positive");
    }
  }
}
