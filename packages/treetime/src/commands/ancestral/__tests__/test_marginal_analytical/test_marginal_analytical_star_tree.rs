#[cfg(test)]
mod tests {
  use super::super::test_marginal_analytical_support::tests::{
    analytical_star_tree_likelihood, run_dense_marginal_get_log_lh, state_index,
  };
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::{GTR, GTRParams};
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use ndarray::array;

  /// Star tree (4 leaves directly from root) with all identical states under JC69.
  #[test]
  fn test_star_tree_analytical_jc69_all_same() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let t = 0.1;
    let observations = [0, 0, 0, 0];

    let expected_lh = analytical_star_tree_likelihood(&gtr, &observations, t);
    let expected_log_lh = expected_lh.ln();

    let newick = "(A:0.1,B:0.1,C:0.1,D:0.1)root;";
    let aln = ">A\nA\n>B\nA\n>C\nA\n>D\nA\n";
    let actual_log_lh = run_dense_marginal_get_log_lh(newick, aln, gtr)?;

    pretty_assert_ulps_eq!(expected_log_lh, actual_log_lh, epsilon = 1e-7);
    Ok(())
  }

  /// Star tree with all four distinct nucleotide states under JC69.
  #[test]
  fn test_star_tree_analytical_jc69_mixed_states() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let t = 0.2;
    let observations = [state_index('A'), state_index('C'), state_index('G'), state_index('T')];

    let expected_lh = analytical_star_tree_likelihood(&gtr, &observations, t);
    let expected_log_lh = expected_lh.ln();

    let newick = "(A:0.2,B:0.2,C:0.2,D:0.2)root;";
    let aln = ">A\nA\n>B\nC\n>C\nG\n>D\nT\n";
    let actual_log_lh = run_dense_marginal_get_log_lh(newick, aln, gtr)?;

    pretty_assert_ulps_eq!(expected_log_lh, actual_log_lh, epsilon = 1e-7);
    Ok(())
  }

  /// Star tree with non-uniform equilibrium frequencies (GTR model).
  #[test]
  fn test_star_tree_analytical_nonuniform_pi() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;
    let gtr = GTR::new(GTRParams {
      alphabet,
      mu: 1.0,
      W: None,
      pi: array![0.1, 0.2, 0.3, 0.4],
    })?;

    let t = 0.15;
    let observations = [state_index('T'), state_index('T'), state_index('G'), state_index('C')];

    let expected_lh = analytical_star_tree_likelihood(&gtr, &observations, t);
    let expected_log_lh = expected_lh.ln();

    let newick = "(A:0.15,B:0.15,C:0.15,D:0.15)root;";
    let aln = ">A\nT\n>B\nT\n>C\nG\n>D\nC\n";
    let actual_log_lh = run_dense_marginal_get_log_lh(newick, aln, gtr)?;

    pretty_assert_ulps_eq!(expected_log_lh, actual_log_lh, epsilon = 1e-7);
    Ok(())
  }
}
