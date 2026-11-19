#[cfg(test)]
mod tests {
  use super::super::test_marginal_stability_support::tests::{
    assert_dense_profile_stable, assert_sparse_profile_stable,
  };
  use crate::test_utils::{run_dense_marginal_with_newick, run_sparse_marginal_with_newick};
  use eyre::Report;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::rare_states_observed(     [0.97, 0.01, 0.01, 0.01],         "(A:0.1,B:0.2)root;", ">A\nCCCC\n>B\nGGGG\n", -10.0)]
  #[case::dominant_state_observed(  [0.97, 0.01, 0.01, 0.01],         "(A:0.1,B:0.2)root;", ">A\nAAAA\n>B\nAAAA\n", 0.0)]
  #[case::extremely_skewed_pi(      [0.9997, 0.0001, 0.0001, 0.0001], "(A:0.1,B:0.1)root;", ">A\nACGT\n>B\nACGT\n", 0.0)]
  #[trace]
  fn test_marginal_stability_near_zero_pi_dense_profiles_normalized(
    #[case] pi: [f64; 4],
    #[case] newick: &str,
    #[case] aln: &str,
    #[case] log_lh_bound: f64,
  ) -> Result<(), Report> {
    let gtr = helpers::nuc_gtr(pi)?;
    let (log_lh, partitions) = run_dense_marginal_with_newick(newick, aln, &gtr)?;

    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh < log_lh_bound, "Log-likelihood should be below {log_lh_bound}: {log_lh}");

    for node_data in partitions.node_states.values() {
      assert_dense_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }

  #[test]
  fn test_marginal_stability_near_zero_pi_sparse_rare_states_profiles_normalized() -> Result<(), Report> {
    let gtr = helpers::nuc_gtr([0.97, 0.01, 0.01, 0.01])?;

    let newick = "(A:0.1,B:0.2)root;";
    let aln = ">A\nCCCC\n>B\nGGGG\n";
    let (log_lh, partitions) = run_sparse_marginal_with_newick(newick, aln, &gtr)?;

    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(
      log_lh < -10.0,
      "Log-likelihood should be strongly negative for rare states: {log_lh}"
    );

    for node_data in partitions.node_states.values() {
      assert_sparse_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::{Alphabet, AlphabetName};
    use crate::gtr::gtr::{GTR, Mu};
    use eyre::Report;
    use ndarray::arr1;

    pub(super) fn nuc_gtr(pi: [f64; 4]) -> Result<GTR, Report> {
      let alphabet = Alphabet::new(AlphabetName::Nuc)?;
      GTR::builder()
        .n_states(alphabet.n_canonical())
        .mu(Mu::FromMatrix(1.0))
        .pi(arr1(&pi))
        .build()
    }
  }
}
