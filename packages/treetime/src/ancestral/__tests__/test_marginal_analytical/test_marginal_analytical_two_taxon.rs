#[cfg(test)]
mod tests {
  use super::super::test_marginal_analytical_support::tests::{analytical_two_taxon_likelihood, state_index};
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::{GTR, Mu};
  use crate::pretty_assert_ulps_eq;
  use crate::test_utils::run_dense_marginal_with_newick;
  use eyre::Report;
  use ndarray::array;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::same_state(           'A', 'A', 0.1,  0.2)]
  #[case::different_states(     'A', 'T', 0.1,  0.2)]
  #[case::asymmetric_branches(  'G', 'G', 0.01, 1.0)]
  #[trace]
  fn test_marginal_analytical_two_taxon_jc69_log_lh_matches_closed_form(
    #[case] state_a: char,
    #[case] state_b: char,
    #[case] t1: f64,
    #[case] t2: f64,
  ) -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;

    let expected_log_lh = analytical_two_taxon_likelihood(&gtr, state_index(state_a), state_index(state_b), t1, t2).ln();

    let newick = format!("(A:{t1},B:{t2})root;");
    let aln = format!(">A\n{state_a}\n>B\n{state_b}\n");
    let (actual_log_lh, _) = run_dense_marginal_with_newick(&newick, &aln, &gtr)?;

    pretty_assert_ulps_eq!(expected_log_lh, actual_log_lh, epsilon = 1e-7);
    Ok(())
  }

  #[test]
  fn test_marginal_analytical_two_taxon_nonuniform_pi_log_lh_matches_closed_form() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let n_states = alphabet.n_canonical();
    let gtr = GTR::builder()
      .n_states(n_states)
      .mu(Mu::FromMatrix(1.0))
      .pi(array![0.4, 0.1, 0.2, 0.3])
      .build()?;

    let t1 = 0.15;
    let t2 = 0.25;

    let expected_lh = analytical_two_taxon_likelihood(&gtr, state_index('A'), state_index('C'), t1, t2);
    let expected_log_lh = expected_lh.ln();

    let newick = "(A:0.15,B:0.25)root;";
    let aln = ">A\nA\n>B\nC\n";
    let (actual_log_lh, _) = run_dense_marginal_with_newick(newick, aln, &gtr)?;

    pretty_assert_ulps_eq!(expected_log_lh, actual_log_lh, epsilon = 1e-7);
    Ok(())
  }

  #[test]
  fn test_marginal_analytical_two_taxon_log_lh_sums_over_positions() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let t1 = 0.1;
    let t2 = 0.2;

    let positions = [
      (state_index('A'), state_index('T')),
      (state_index('C'), state_index('C')),
      (state_index('G'), state_index('A')),
    ];

    let mut expected_log_lh = 0.0;
    for (obs_a, obs_b) in positions {
      let lh = analytical_two_taxon_likelihood(&gtr, obs_a, obs_b, t1, t2);
      expected_log_lh += lh.ln();
    }

    let newick = "(A:0.1,B:0.2)root;";
    let aln = ">A\nACG\n>B\nTCA\n";
    let (actual_log_lh, _) = run_dense_marginal_with_newick(newick, aln, &gtr)?;

    pretty_assert_ulps_eq!(expected_log_lh, actual_log_lh, epsilon = 1e-7);
    Ok(())
  }
}
