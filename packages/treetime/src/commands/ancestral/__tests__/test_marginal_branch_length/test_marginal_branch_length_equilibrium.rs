#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::{GTR, GTRParams};
  use crate::pretty_assert_ulps_eq;
  use crate::test_utils::{run_dense_marginal_with_newick, run_sparse_marginal_with_newick};
  use eyre::Report;
  use ndarray::array;

  /// At very long branch lengths, Felsenstein site likelihood converges to the product of equilibrium frequencies (dense partition).
  #[test]
  fn test_equilibrium_convergence_dense() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let t = 100.0;
    let expected_equilibrium_log_lh = (0.25 * 0.25_f64).ln();

    let newick = format!("(A:{t},B:{t})root;");
    let aln = ">A\nA\n>B\nT\n";
    let actual_log_lh = run_dense_marginal_with_newick(&newick, aln, gtr)?;

    pretty_assert_ulps_eq!(expected_equilibrium_log_lh, actual_log_lh, epsilon = 1e-6);
    Ok(())
  }

  /// Equilibrium convergence test using sparse partition.
  #[test]
  fn test_equilibrium_convergence_sparse() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let t = 100.0;
    let expected_equilibrium_log_lh = (0.25 * 0.25_f64).ln();

    let newick = format!("(A:{t},B:{t})root;");
    let aln = ">A\nA\n>B\nT\n";
    let actual_log_lh = run_sparse_marginal_with_newick(&newick, aln, gtr)?;

    pretty_assert_ulps_eq!(expected_equilibrium_log_lh, actual_log_lh, epsilon = 1e-6);
    Ok(())
  }

  /// Equilibrium convergence with non-uniform equilibrium frequencies (GTR model).
  #[test]
  fn test_equilibrium_convergence_nonuniform_pi() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let n_states = alphabet.n_canonical();
    let gtr = GTR::new(GTRParams {
      n_states,
      mu: 1.0,
      W: None,
      pi: array![0.4, 0.1, 0.2, 0.3],
    })?;

    let t = 100.0;
    let expected_equilibrium_log_lh = (0.4 * 0.2_f64).ln();

    let newick = format!("(A:{t},B:{t})root;");
    let aln = ">A\nA\n>B\nG\n";
    let actual_log_lh = run_dense_marginal_with_newick(&newick, aln, gtr)?;

    pretty_assert_ulps_eq!(expected_equilibrium_log_lh, actual_log_lh, epsilon = 1e-6);
    Ok(())
  }

  /// Multi-position equilibrium convergence under JC69.
  #[test]
  fn test_equilibrium_convergence_multiple_positions() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let t = 100.0;
    let expected_equilibrium_log_lh = 3.0 * (0.25 * 0.25_f64).ln();

    let newick = format!("(A:{t},B:{t})root;");
    let aln = ">A\nACG\n>B\nTCA\n";
    let actual_log_lh = run_dense_marginal_with_newick(&newick, aln, gtr)?;

    pretty_assert_ulps_eq!(expected_equilibrium_log_lh, actual_log_lh, epsilon = 1e-6);
    Ok(())
  }

  /// Equilibrium convergence on a 4-leaf star tree (polytomy) under JC69.
  #[test]
  fn test_equilibrium_star_tree() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let t = 100.0;
    let expected_equilibrium_log_lh = (0.25_f64.powi(4)).ln();

    let newick = format!("(A:{t},B:{t},C:{t},D:{t})root;");
    let aln = ">A\nA\n>B\nC\n>C\nG\n>D\nT\n";
    let actual_log_lh = run_dense_marginal_with_newick(&newick, aln, gtr)?;

    pretty_assert_ulps_eq!(expected_equilibrium_log_lh, actual_log_lh, epsilon = 1e-6);
    Ok(())
  }

  /// Numerical stability at extremely short branch lengths.
  #[test]
  fn test_branch_length_sensitivity_near_zero() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let t = 1e-8;

    let newick = format!("(A:{t},B:{t})root;");
    let aln = ">A\nA\n>B\nA\n";
    let log_lh = run_dense_marginal_with_newick(&newick, aln, gtr)?;

    assert!(log_lh.is_finite(), "Log-likelihood should be finite at t={t}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive at t={t}");

    let approx_expected = 0.25_f64.ln();
    assert!(
      (log_lh - approx_expected).abs() < 1e-2,
      "At very short branches, log_lh should be close to ln(pi[A]). Expected ~{approx_expected}, got {log_lh}"
    );

    Ok(())
  }

  /// Dense and sparse partitions produce identical Felsenstein log-likelihoods across branch lengths.
  #[test]
  fn test_dense_sparse_consistency_across_branch_lengths() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let branch_lengths = [0.01, 0.1, 0.5, 1.0, 5.0, 20.0];

    let violation = branch_lengths
      .iter()
      .map(|t| {
        let newick = format!("(A:{t},B:{t})root;");
        let aln = ">A\nACGTACGT\n>B\nTGCATGCA\n";
        let dense_log_lh = run_dense_marginal_with_newick(&newick, aln, gtr.clone())?;
        let sparse_log_lh = run_sparse_marginal_with_newick(&newick, aln, gtr.clone())?;
        Ok::<_, Report>((*t, dense_log_lh, sparse_log_lh))
      })
      .collect::<Result<Vec<_>, _>>()?
      .into_iter()
      .find(|(_, dense_log_lh, sparse_log_lh)| (dense_log_lh - sparse_log_lh).abs() >= 1e-10);

    if let Some((t, dense_log_lh, sparse_log_lh)) = violation {
      panic!("Dense and sparse should match. At t={t}: dense={dense_log_lh}, sparse={sparse_log_lh}");
    }

    Ok(())
  }
}
