#[cfg(test)]
mod tests {
  use super::super::test_marginal_analytical_support::tests::{
    analytical_two_taxon_likelihood, run_dense_marginal_get_log_lh, state_index,
  };
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::{GTR, GTRParams};
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use ndarray::array;

  /// Two-taxon tree with identical leaf states under JC69 (Jukes-Cantor 1969).
  ///
  /// Tree: `(A:0.1, B:0.2)root;` both leaves observe state `A`.
  ///
  /// Analytical: `L = sum_s pi[s] * P(A | s, 0.1) * P(A | s, 0.2)`.
  #[test]
  fn test_two_taxon_analytical_jc69_same_state() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let t1 = 0.1;
    let t2 = 0.2;

    let expected_lh = analytical_two_taxon_likelihood(&gtr, 0, 0, t1, t2);
    let expected_log_lh = expected_lh.ln();

    let newick = "(A:0.1,B:0.2)root;";
    let aln = ">A\nA\n>B\nA\n";
    let actual_log_lh = run_dense_marginal_get_log_lh(newick, aln, gtr)?;

    pretty_assert_ulps_eq!(expected_log_lh, actual_log_lh, epsilon = 1e-7);
    Ok(())
  }

  /// Two-taxon tree with different leaf states under JC69.
  ///
  /// Tree: `(A:0.1, B:0.2)root;` leaf `A` observes `A`, leaf `B` observes `T`.
  #[test]
  fn test_two_taxon_analytical_jc69_different_states() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let t1 = 0.1;
    let t2 = 0.2;

    let expected_lh = analytical_two_taxon_likelihood(&gtr, state_index('A'), state_index('T'), t1, t2);
    let expected_log_lh = expected_lh.ln();

    let newick = "(A:0.1,B:0.2)root;";
    let aln = ">A\nA\n>B\nT\n";
    let actual_log_lh = run_dense_marginal_get_log_lh(newick, aln, gtr)?;

    pretty_assert_ulps_eq!(expected_log_lh, actual_log_lh, epsilon = 1e-7);
    Ok(())
  }

  /// Two-taxon tree with non-uniform equilibrium frequencies (GTR model).
  #[test]
  fn test_two_taxon_analytical_nonuniform_pi() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;
    let n_states = alphabet.n_canonical();
    let gtr = GTR::new(GTRParams {
      n_states,
      mu: 1.0,
      W: None,
      pi: array![0.4, 0.1, 0.2, 0.3],
    })?;

    let t1 = 0.15;
    let t2 = 0.25;

    let expected_lh = analytical_two_taxon_likelihood(&gtr, state_index('A'), state_index('C'), t1, t2);
    let expected_log_lh = expected_lh.ln();

    let newick = "(A:0.15,B:0.25)root;";
    let aln = ">A\nA\n>B\nC\n";
    let actual_log_lh = run_dense_marginal_get_log_lh(newick, aln, gtr)?;

    pretty_assert_ulps_eq!(expected_log_lh, actual_log_lh, epsilon = 1e-7);
    Ok(())
  }

  /// Two-taxon tree with a 3-position alignment under JC69.
  #[test]
  fn test_two_taxon_analytical_multiple_positions() -> Result<(), Report> {
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
    let actual_log_lh = run_dense_marginal_get_log_lh(newick, aln, gtr)?;

    pretty_assert_ulps_eq!(expected_log_lh, actual_log_lh, epsilon = 1e-7);
    Ok(())
  }

  /// Two-taxon tree with highly asymmetric branch lengths under JC69.
  #[test]
  fn test_two_taxon_analytical_asymmetric_branches() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let t1 = 0.01;
    let t2 = 1.0;

    let expected_lh = analytical_two_taxon_likelihood(&gtr, state_index('G'), state_index('G'), t1, t2);
    let expected_log_lh = expected_lh.ln();

    let newick = "(A:0.01,B:1.0)root;";
    let aln = ">A\nG\n>B\nG\n";
    let actual_log_lh = run_dense_marginal_get_log_lh(newick, aln, gtr)?;

    pretty_assert_ulps_eq!(expected_log_lh, actual_log_lh, epsilon = 1e-7);
    Ok(())
  }
}
