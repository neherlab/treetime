#[cfg(test)]
mod tests {
  use super::super::test_marginal_analytical_support::tests::{analytical_three_taxon_likelihood, state_index};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::test_utils::run_dense_marginal_with_newick;
  use eyre::Report;
  use rstest::rstest;
  use treetime_utils::pretty_assert_ulps_eq;

  #[rstest]
  fn test_marginal_analytical_three_taxon_log_lh_matches_closed_form_for_all_states(
    #[values('A', 'C', 'G', 'T')] state_a: char,
    #[values('A', 'C', 'G', 'T')] state_b: char,
    #[values('A', 'C', 'G', 'T')] state_c: char,
  ) -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let expected_log_lh = analytical_three_taxon_likelihood(
      &gtr,
      state_index(state_a),
      state_index(state_b),
      state_index(state_c),
      0.1,
      0.2,
      0.15,
      0.25,
    )
    .ln();

    let newick = "((A:0.1,B:0.2)AB:0.15,C:0.25)root;";
    let aln = format!(">A\n{state_a}\n>B\n{state_b}\n>C\n{state_c}\n");
    let (actual_log_lh, _) = run_dense_marginal_with_newick(newick, &aln, &gtr)?;

    pretty_assert_ulps_eq!(expected_log_lh, actual_log_lh, epsilon = 1e-7);
    Ok(())
  }
}
