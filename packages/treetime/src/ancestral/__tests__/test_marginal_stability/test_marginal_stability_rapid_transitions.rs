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
  #[case::high_mutation_rate(          (10.0,  [0.25, 0.25, 0.25, 0.25]), ("(A:0.1,B:0.1)root;",                   ">A\nACGT\n>B\nTGCA\n"))]
  #[case::very_high_mutation_rate(     (100.0, [0.25, 0.25, 0.25, 0.25]), ("(A:0.01,B:0.01)root;",                 ">A\nAAAA\n>B\nTTTT\n"))]
  #[case::high_mutation_nonuniform_pi( (10.0,  [0.4, 0.1, 0.2, 0.3]),     ("((A:0.05,B:0.05)AB:0.1,C:0.15)root;", ">A\nACGT\n>B\nTGCA\n>C\nGGGG\n"))]
  #[case::combined_extreme_parameters( (50.0,  [0.9, 0.03, 0.04, 0.03]),  ("(A:1e-8,B:1e-8)root;",                 ">A\nCCCC\n>B\nGGGG\n"))]
  #[trace]
  fn test_marginal_stability_rapid_transitions_dense_profiles_normalized(
    #[case] (mu, pi): (f64, [f64; 4]),
    #[case] (newick, aln): (&str, &str),
  ) -> Result<(), Report> {
    let gtr = helpers::nuc_gtr(mu, pi)?;
    let (log_lh, partitions) = run_dense_marginal_with_newick(newick, aln, &gtr)?;

    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");

    for node_data in partitions.node_states.values() {
      assert_dense_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::high_mutation_rate(          (10.0, [0.25, 0.25, 0.25, 0.25]), ("(A:0.1,B:0.1)root;",  ">A\nACGT\n>B\nTGCA\n"))]
  #[case::combined_extreme_parameters( (50.0, [0.9, 0.03, 0.04, 0.03]),  ("(A:1e-8,B:1e-8)root;", ">A\nCCCC\n>B\nGGGG\n"))]
  #[trace]
  fn test_marginal_stability_rapid_transitions_sparse_profiles_normalized(
    #[case] (mu, pi): (f64, [f64; 4]),
    #[case] (newick, aln): (&str, &str),
  ) -> Result<(), Report> {
    let gtr = helpers::nuc_gtr(mu, pi)?;
    let (log_lh, partitions) = run_sparse_marginal_with_newick(newick, aln, &gtr)?;

    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");

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

    pub(super) fn nuc_gtr(mu: f64, pi: [f64; 4]) -> Result<GTR, Report> {
      let alphabet = Alphabet::new(AlphabetName::Nuc)?;
      GTR::builder()
        .n_states(alphabet.n_canonical())
        .mu(Mu::FromMatrix(mu))
        .pi(arr1(&pi))
        .build()
    }
  }
}
